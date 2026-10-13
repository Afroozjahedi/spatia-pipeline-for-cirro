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
   (columns: sample, file, sampleIndex, ...) -- no more guessing raw-image
   paths either. ds.samplesheet carries the metadata uploaded via Cirro's
   native "Upload Samplesheet" feature (columns: sample, group, ...,
   possibly also a redundant "file" column from the original upload, which
   we drop in favor of ds.files' resolved path). We merge the two on
   "sample", write the result to this dataset's own config/ folder, and
   hand Nextflow the final path via ds.add_param() -- so process-input.json
   no longer needs a "samplesheet" key at all.

   Confidence this fixes it: ~90% (based on reading the actual installed
   SDK source, not on inference from error text). The ~10% uncertainty is
   whether ds.samplesheet is reliably populated for this account's Cirro
   the way the SDK docs describe -- the next run will confirm or refute
   that directly.
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


def _folder_filename_key(path: str) -> str:
    """
    Last two path components ("folder/filename.tif"), used as a
    match key that survives differing S3 prefixes (ds.files' full S3 URI
    vs. an uploaded samplesheet's dataset-relative path) while still
    disambiguating files that share a bare filename across different
    SOURCE DATASETS.

    ADDED 2026-10-12, multi-dataset combined runs: when this pipeline runs
    against several selected source datasets at once (Afrouz's CRC
    TMA_A / TMA_B / DII_TMA_A / DII_TMA_B layout), ds.files is the union of
    every selected dataset's raw files. Two different TMAs can each have a
    core that happens to be named the same thing (e.g. both TMA_A and
    TMA_B contain a file literally called reg034_X01_Y01_Z09.tif) -- these
    are different physical tissue cores, not duplicates. Matching on bare
    "sample" alone (the old behavior) can't tell them apart: a left-merge
    would hand ONE metadata row's group to BOTH physical files, silently
    mislabeling the wrong one -- and, worse, re-creating the exact
    image_id collision that run_pipeline.py's per-dataset samplesheets
    were built to avoid (see samplesheets/ and the 8 cross-dataset
    collisions Afrouz identified and corrected 2026-10-12). Matching on
    folder+filename instead is safe because each uploaded per-dataset
    samplesheet's "file" column already encodes which folder it came from
    (e.g. "TMA_A/reg034_X01_Y01_Z09.tif"), and ds.files' own S3 path ends
    the same way.
    """
    parts = str(path).rstrip("/").split("/")
    return "/".join(parts[-2:]) if len(parts) >= 2 else (parts[-1] if parts else "")


def build_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    """
    Merge Cirro's auto-resolved file listing (ds.files) with the
    sample-level metadata from the uploaded samplesheet(s) (ds.samplesheet),
    matched by folder+filename (see _folder_filename_key) rather than bare
    "sample" -- required for correctness once this pipeline runs against
    multiple source datasets that can share bare sample names (see that
    function's docstring). ds.files' own "file" column (the authoritative
    S3 path) is what ends up in the output; the uploaded samplesheet's own
    "file" column is used only to build the match key, then dropped.

    Rows with no metadata match are EXCLUDED from the returned samplesheet
    (not kept with group=NaN) -- run_pipeline.py's _apply_samplesheet()
    builds image_experiment_group_map from every row in this output, so an
    unmatched row would otherwise propagate as a literal missing/NaN group
    instead of being left out of the run.
    """
    if "file" not in ds.samplesheet.columns:
        raise ValueError(
            "build_samplesheet() requires the uploaded samplesheet to include "
            "a 'file' column (dataset-relative path, e.g. 'TMA_A/reg034_X01_Y01_Z09.tif') "
            "-- it's needed to disambiguate files that share a bare sample name "
            "across different source datasets when this pipeline runs against "
            "multiple selected datasets at once."
        )

    meta = ds.samplesheet.copy()
    meta["_match_key"] = meta["file"].map(_folder_filename_key)
    # Drop "file" (superseded by ds.files' own authoritative S3 path) AND
    # "sample" (ds.files' "sample" column is kept as authoritative too --
    # without dropping this, pandas would suffix both sides' "sample"
    # column into sample_x/sample_y instead of leaving a clean "sample"
    # column, which would break run_pipeline.py's _resolve_samplesheet_columns()).
    meta = meta.drop(columns=[c for c in ("file", "sample") if c in meta.columns])

    files = ds.files.copy()
    files["_match_key"] = files["file"].map(_folder_filename_key)

    merged = pd.merge(files, meta, on="_match_key", how="left", validate="many_to_one")
    merged = merged.drop(columns=["_match_key"])

    n_total = len(merged)
    if "group" in merged.columns:
        unmatched = merged["group"].isna()
    else:
        unmatched = pd.Series([False] * n_total, index=merged.index)
    n_unmatched = int(unmatched.sum())
    if n_unmatched:
        ds.logger.info(
            f"WARNING: {n_unmatched} of {n_total} files (matched by folder/filename) have no "
            "row in any uploaded samplesheet -- EXCLUDING them from this run rather than "
            "processing them with an unknown group. If this count looks too high, check that "
            "EVERY selected dataset has its own uploaded samplesheet with a 'file' column in "
            "the 'folder/filename.tif' form."
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
