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


def build_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    """
    Merge Cirro's auto-resolved file listing (ds.files) with the
    sample-level metadata from the uploaded samplesheet (ds.samplesheet),
    on "sample". ds.files' "file" column is authoritative for the actual
    S3 path -- if ds.samplesheet also has a "file" column (left over from
    the original upload), it's dropped rather than merged, to avoid a
    file_x/file_y collision and to avoid overriding Cirro's own resolved
    path with whatever was typed into the upload.
    """
    meta = ds.samplesheet.drop(columns=[c for c in ("file",) if c in ds.samplesheet.columns])
    merged = pd.merge(ds.files, meta, on="sample", how="left", validate="many_to_one")
    missing_group = merged["group"].isna().sum() if "group" in merged.columns else 0
    if missing_group:
        ds.logger.info(
            f"WARNING: {missing_group} of {len(merged)} samples have no metadata match "
            "in the uploaded samplesheet (check the 'sample' values line up exactly)."
        )
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
