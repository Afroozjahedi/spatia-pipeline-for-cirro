#!/usr/bin/env python3
"""
Decode the SPATIA params YAML and samplesheet CSV that Cirro embeds as
inline base64 data URIs, and write them to this dataset's own config/
folder as real files.

Why this exists (found 2026-09-28 by comparing a run's actual debug log
against process-input.json): process-form.json's two fields use
"format": "file" (direct browser upload). Cirro renders a "format": "file"
form value into process-input.json as an inline base64 data URI --
  data:application/x-yaml;name=crc_tma_full_pipeline.params.yaml;base64,....
-- per JSON Schema's file-format convention, NOT as an S3 path. Nextflow's
Channel.fromPath() has no concept of a data: URI: it tried to treat the
whole string, base64 blob included, as a literal filesystem path
("File name too long"), then as a literal S3 object key ("400 Bad
Request"). Confirmed directly against the run's own nextflow debug log
and the actual params.json Cirro generated -- not inferred.

Cirro's own nf-core/sopa pipeline (CirroBio/Cirro-pipelines) hits the same
"file picked in the form" situation for phenocycler/ome_tif and solves it
with exactly this pattern: a preprocess.py step (which Cirro already runs
before Nextflow launches -- our logs show "RUNNING PREPROCESS" every run,
just "Skipping preprocess script" since we didn't have one) reads the raw
form value from ds.params, and writes the real content to a FIXED S3
location that process-input.json maps separately via the
"$.dataset.s3|/config/..." convention -- e.g. sopa writes its derived
samplesheet to ds.params["input"], which process-input.json maps to
"$.dataset.s3|/config/manifest.csv". main.nf then reads that fixed,
already-resolved S3 URI via Channel.fromPath(), which Nextflow can
actually stage.

This script does the same thing for our two uploads:
  config_file_upload / samplesheet_upload  (raw data URIs, from the form)
  -> decoded and written to ->
  config / samplesheet  (fixed S3 paths process-input.json already points
                          main.nf's params.config / params.samplesheet at)

Then the raw upload params are removed so Nextflow never sees the data URIs.
"""
import base64

import boto3

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


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()

    config_bytes = decode_data_uri(ds.params["config_file_upload"])
    samplesheet_bytes = decode_data_uri(ds.params["samplesheet_upload"])

    write_to_s3(config_bytes, ds.params["config"])
    write_to_s3(samplesheet_bytes, ds.params["samplesheet"])

    ds.logger.info(f"Wrote decoded config file to {ds.params['config']}")
    ds.logger.info(f"Wrote decoded samplesheet to {ds.params['samplesheet']}")

    # Raw uploads only exist to carry the data URIs into this script --
    # remove them so Nextflow's params never contain the huge blobs
    # (mirrors sopa's own cleanup of its spatial_data/image_file params).
    ds.remove_param("config_file_upload", force=True)
    ds.remove_param("samplesheet_upload", force=True)
