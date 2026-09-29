#!/usr/bin/env python3
"""
Decode the SPATIA pipeline-params YAML that Cirro embeds as an inline
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

Cirro's own nf-core/sopa pipeline (CirroBio/Cirro-pipelines) hits the same
"file picked in the form" situation for phenocycler/ome_tif and solves it
with exactly this pattern: a preprocess.py step (which Cirro already runs
before Nextflow launches -- our logs show "RUNNING PREPROCESS" every run,
just "Skipping preprocess script" since we didn't have one) reads the raw
form value from ds.params, and writes the real content to a FIXED S3
location that process-input.json maps separately via the
"$.dataset.s3|/config/..." convention. main.nf then reads that fixed,
already-resolved S3 URI via Channel.fromPath(), which Nextflow can
actually stage.

UPDATED 2026-09-29 -- samplesheet handling REMOVED from this script.
Per Cirro admin (Dima), the samplesheet no longer comes through the
run-launch form at all: it now lives on the dataset itself, uploaded via
Cirro's native "Upload Samplesheet" feature. That native upload is a real
file at a real S3 path from the start -- it never becomes a data URI, so
there is nothing for this script to decode for it. process-input.json
reads it directly via a fixed "$.dataset.s3|/samplesheet.csv" mapping.
This script now only exists for config_file, which is still a raw
"format": "file" browser upload and still needs the decode-and-relocate
treatment described above.
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
    write_to_s3(config_bytes, ds.params["config"])
    ds.logger.info(f"Wrote decoded config file to {ds.params['config']}")

    # Raw upload param only exists to carry the data URI into this script --
    # remove it so Nextflow's params never contain the huge blob
    # (mirrors sopa's own cleanup of its spatial_data/image_file params).
    ds.remove_param("config_file_upload", force=True)
