#!/usr/bin/env python3
"""
Two jobs, both run by Cirro before Nextflow launches:

1. Decode the SPATIA pipeline-params YAML that Cirro embeds as an inline
   base64 data URI, and write it to this dataset's own config/ folder as a
   real file.

   Why this exists (found 2026-09-28 by comparing a run's actual debug log
   against process-input.json): process-form.json's "config_file" field uses
   "format": "file" (direct browser upload). Cirro renders a "format": "file"
   form value into process-input.json as an inline base64 data URI --
     data:application/x-yaml;name=crc_tma_full_pipeline.params.yaml;base64,....
   -- per JSON Schema's file-format convention, NOT as an S3 path. Nextflow's
   Channel.fromPath() has no concept of a data: URI: it tried to treat the
   whole string, base64 blob included, as a literal filesystem path
   ("File name too long"), then as a literal S3 object key ("400 Bad
   Request"). Confirmed directly against the run's own nextflow debug log
   and the actual params.json Cirro generated -- not inferred.

2. UPDATED 2026-09-30 (v2) -- build the samplesheet.csv that run_pipeline.py
   needs, from Cirro's own dataset metadata, instead of guessing a static
   S3 path for it in process-input.json.

   Why: every earlier attempt to point process-input.json's "samplesheet"
   key at a fixed "$.dataset.s3|/samplesheet.csv" (or later, a hardcoded
   literal S3 URI) 404'd. Root cause, found by reading the actual installed
   cirro SDK source (cirro/helpers/preprocess_dataset.py in
   CirroBio/Cirro-SDK-Python) rather than guessing again: Cirro stores the
   dataset's sample metadata and file listing at
   "<dataset_root>/config/samplesheet.csv" and
   "<dataset_root>/config/files.csv" respectively -- one level deeper than
   we assumed ("config/" was missing) -- and PreprocessDataset.from_running()
   already loads both of these automatically into ds.samplesheet and
   ds.files for the CURRENT running dataset (dataset_root = the
   PW_S3_DATASET env var), with no hardcoded dataset ID needed. This is
   also exactly the documented, intended pattern: "One of the most common
   tasks performed in the preprocess script is to construct a sample
   sheet ... listing all of the files from the input dataset."
   (docs.cirro.bio/pipelines/preprocess-script/).

   ds.files carries Cirro's own resolved, authoritative S3 path per sample
   (columns: sample, file, dataset, sampleIndex, ...) -- no more guessing
   raw-image paths either. ds.samplesheet carries the metadata uploaded via
   Cirro's native "Upload Samplesheet" feature.

   CONFIRMED 2026-10-13, by directly inspecting a real combined 4-dataset
   run's config/samplesheet.csv: ds.samplesheet NEVER carries a "file"
   column, even when the uploaded CSV had one -- it only ever has the
   sample-level columns Cirro itself assigns (sample, group, ...). An
   earlier version of build_samplesheet() hard-required "file" in
   ds.samplesheet and matched on folder+filename; that assumption was
   wrong for every multi-dataset run, not just a one-off, so it's been
   replaced with a plain "sample" match against ds.files (see that
   function's docstring for the current approach and its safety check).
   We merge the two on "sample", write the result to this dataset's own
   config/ folder, and hand Nextflow the final path via ds.add_param() --
   so process-input.json no longer needs a "samplesheet" key at all.
"""
import base64

import boto3
import pandas as pd

from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.models.s3_path import S3Path


def decode_data_uri(data_uri: str) -> bytes:
    """Decode a `data:<mediatype>;base64,<payload>` URI to raw bytes."""
    header, sep, payload = data_uri.partition(",")
    if not sep:
        raise ValueError(f"Not a data URI (no comma found): {data_uri[:80]!r}...")
    if ";base64" not in header:
        raise ValueError(f"Expected a base64 data URI, got header: {header!r}")
    return base64.b64decode(payload)


def write_to_s3(data: bytes, s3_uri: str) -> None:
    path = S3Path(s3_uri)
    if not path.valid:
        raise ValueError(f"Not a valid s3:// URI: {s3_uri!r}")
    boto3.client("s3").put_object(Bucket=path.bucket, Key=path.key, Body=data)


def build_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    """
    Merge Cirro's auto-resolved file listing (ds.files) with the
    sample-level metadata from the uploaded samplesheet(s) (ds.samplesheet),
    matched on "sample".

    REVISED 2026-10-13: the previous version matched on folder+filename and
    hard-required ds.samplesheet to carry a "file" column. Confirmed by
    direct inspection of a real combined 4-dataset run (cc786886, 2026-10-13)
    that ds.samplesheet NEVER carries "file" -- it only ever has the
    sample-level columns Cirro itself assigns (sample, group, ...).
    File-path info lives in ds.files instead, which already has it
    natively, per physical file, alongside a "dataset" column. The old
    guard assumed a mechanism ds.samplesheet doesn't have, so it was wrong
    for every multi-dataset run, not just this one.

    Matching is on "sample" alone now, which is safe as long as ds.files'
    "sample" values are unique. That held for the 2026-10-13 run: each of
    the 8 known cross-dataset filename collisions (same bare sample name
    physically present in two TMAs) appeared exactly once in ds.files,
    apparently because Cirro's combined "files" listing only surfaces
    files that have matching assigned sample metadata in their OWN
    dataset's uploaded samplesheet -- which excludes the stray duplicate
    once its row was removed from the wrong per-dataset samplesheet
    (2026-10-12 fix). If a future dataset combination ever has a genuine
    duplicate "sample" in ds.files, this raises explicitly below instead
    of silently mis-assigning one metadata row's group to the wrong
    physical file.
    """
    files = ds.files.copy()
    dup_files_samples = (
        files["sample"][files["sample"].duplicated(keep=False)].unique().tolist()
    )
    if dup_files_samples:
        raise ValueError(
            "build_samplesheet() found duplicate 'sample' values in ds.files, "
            "which makes matching sample metadata to files ambiguous: "
            f"{dup_files_samples}. This usually means two selected datasets both "
            "contain a physical file with the same bare sample name -- remove or "
            "rename the stray duplicate in Cirro before re-running."
        )

    meta = ds.samplesheet.copy()
    dup_meta_samples = (
        meta["sample"][meta["sample"].duplicated(keep=False)].unique().tolist()
    )
    if dup_meta_samples:
        raise ValueError(
            "build_samplesheet() found duplicate 'sample' values in the uploaded "
            f"samplesheet metadata (ds.samplesheet): {dup_meta_samples}. Each "
            "sample name must be unique across every selected dataset's uploaded "
            "samplesheet."
        )
    # Drop any "file" column the uploaded samplesheet happens to carry --
    # ds.files' own "file" (the authoritative S3 path) is what's kept.
    meta = meta.drop(columns=[c for c in ("file",) if c in meta.columns])

    merged = pd.merge(files, meta, on="sample", how="left", validate="one_to_one")

    n_total = len(merged)
    if "group" in merged.columns:
        unmatched = merged["group"].isna()
    else:
        unmatched = pd.Series([False] * n_total, index=merged.index)
    n_unmatched = int(unmatched.sum())
    if n_unmatched:
        ds.logger.info(
            f"WARNING: {n_unmatched} of {n_total} files (matched by 'sample') have no "
            "row in any uploaded samplesheet -- EXCLUDING them from this run rather than "
            "processing them with an unknown group. If this count looks too high, check "
            "that EVERY selected dataset has its own uploaded samplesheet covering its "
            "samples."
        )
        for _, r in merged[unmatched].iterrows():
            ds.logger.info(f"  excluded (no samplesheet match): {r.get('file', '?')}")
        merged = merged[~unmatched].reset_index(drop=True)

    return merged


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()

    config_bytes = decode_data_uri(ds.params["config_file_upload"])
    write_to_s3(config_bytes, ds.params["config"])
    ds.logger.info(f"Wrote decoded config file to {ds.params['config']}")

    # Raw upload param only exists to carry the data URI into this script --
    # remove it so Nextflow's params never contain the huge blob
    # (mirrors sopa's own cleanup of its spatial_data/image_file params).
    ds.remove_param("config_file_upload", force=True)

    samplesheet_df = build_samplesheet(ds)
    samplesheet_s3_path = f"{ds.dataset_root}/config/pipeline_samplesheet.csv"
    write_to_s3(samplesheet_df.to_csv(index=False).encode(), samplesheet_s3_path)
    ds.logger.info(
        f"Wrote merged samplesheet ({len(samplesheet_df)} rows, "
        f"columns={list(samplesheet_df.columns)}) to {samplesheet_s3_path}"
    )
    ds.add_param("samplesheet", samplesheet_s3_path, overwrite=True)
