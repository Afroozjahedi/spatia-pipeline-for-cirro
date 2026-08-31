#!/usr/bin/env python3
"""
build_image_patient_map.py
============================
Generates the {image_id: patient_id} JSON file survival.py's
`analysis.survival.image_patient_map` config key can point at, for
TMA-style cohorts where a patient-level annotation file lists which core
numbers ("TMA spot / region", or whatever your column is actually called)
belong to each patient, and processed images are named with a matching
region-number pattern (e.g. "CLR_reg001_A" -> region 1).

Why this exists
-----------------
survival.py was rewritten to be dataset-agnostic — it no longer auto-derives
image->patient from a "TMA spot / region"-shaped column (that assumption
doesn't hold for non-TMA cohorts). For a real TMA cohort, you still need
that mapping; this script builds it once from your annotation file instead
of either hand-typing ~100+ entries or baking TMA-specific parsing back into
the generic survival module.

Not specific to CRC or to the exact column names used in this repo's example
config — column names and the region-number pattern are all CLI arguments.

Usage
-----
    python build_image_patient_map.py \\
        --annotation-file /path/to/patient_with_tls_class.csv \\
        --patient-id-col Patient \\
        --spot-col "TMA spot / region" \\
        --image-dir /path/to/matched_cells_dir \\
        --output image_patient_map.json

Then in your experiment YAML:
    analysis:
      survival:
        image_patient_map: "image_patient_map.json"
"""

import argparse
import json
import os
import re
import sys

import pandas as pd


def _spot_to_regions(spot_str, region_sep: str) -> list:
    """'1,2' -> [1, 2] (or whatever separator your annotation file uses)."""
    try:
        return [int(s.strip()) for s in str(spot_str).split(region_sep) if s.strip()]
    except (TypeError, ValueError):
        return []


def build_map(
    annotation_file: str,
    patient_id_col: str,
    spot_col: str,
    region_sep: str,
    image_dir: str,
    region_regex: str,
    file_suffix: str,
) -> dict:
    surv_df = pd.read_csv(annotation_file)
    for col in (patient_id_col, spot_col):
        if col not in surv_df.columns:
            raise ValueError(
                f"Column {col!r} not found in {annotation_file} — "
                f"available columns: {list(surv_df.columns)}"
            )

    image_ids = sorted(
        f.replace(file_suffix, "")
        for f in os.listdir(image_dir)
        if f.endswith(file_suffix) and not f.startswith("._")
    )
    if not image_ids:
        raise ValueError(f"No files ending in {file_suffix!r} found in {image_dir}")

    pattern = re.compile(region_regex)
    image_region = {}
    unmatched_images = []
    for image_id in image_ids:
        m = pattern.search(image_id)
        if m:
            image_region[image_id] = int(m.group(1))
        else:
            unmatched_images.append(image_id)

    image_patient_map = {}
    collisions = []
    patients_with_no_images = []

    for _, row in surv_df.dropna(subset=[patient_id_col]).iterrows():
        patient_id = str(row[patient_id_col])
        regions = _spot_to_regions(row[spot_col], region_sep)
        if not regions:
            continue
        matched_any = False
        for image_id, region in image_region.items():
            if region in regions:
                matched_any = True
                if image_id in image_patient_map and image_patient_map[image_id] != patient_id:
                    collisions.append((image_id, image_patient_map[image_id], patient_id))
                image_patient_map[image_id] = patient_id
        if not matched_any:
            patients_with_no_images.append(patient_id)

    mapped_images = set(image_patient_map.keys())
    all_regioned_images = set(image_region.keys())
    images_with_no_patient = sorted(all_regioned_images - mapped_images)

    print(f"[build_image_patient_map] Annotation file : {annotation_file}  ({len(surv_df)} patients)")
    print(f"[build_image_patient_map] Image dir        : {image_dir}  ({len(image_ids)} images)")
    print(f"[build_image_patient_map] Matched          : {len(image_patient_map)} image -> patient entries")

    if unmatched_images:
        print(f"[build_image_patient_map] ⚠️  {len(unmatched_images)} image(s) didn't match "
              f"--region-regex {region_regex!r} at all (not counted as unmatched-to-patient, "
              f"just unparseable): {unmatched_images[:10]}"
              f"{' ...' if len(unmatched_images) > 10 else ''}")
    if images_with_no_patient:
        print(f"[build_image_patient_map] ⚠️  {len(images_with_no_patient)} image(s) parsed a region "
              f"number but no patient's {spot_col!r} list claims that region: "
              f"{images_with_no_patient[:10]}{' ...' if len(images_with_no_patient) > 10 else ''}")
    if patients_with_no_images:
        print(f"[build_image_patient_map] ⚠️  {len(patients_with_no_images)} patient(s) have a "
              f"{spot_col!r} entry but zero matching images in --image-dir: "
              f"{patients_with_no_images[:10]}{' ...' if len(patients_with_no_images) > 10 else ''}")
    if collisions:
        print(f"[build_image_patient_map] ⚠️  {len(collisions)} image(s) matched region numbers "
              f"claimed by MORE THAN ONE patient (last one wins, likely a data problem — check "
              f"your annotation file for duplicate/overlapping region numbers): {collisions[:10]}")

    return image_patient_map


def main():
    p = argparse.ArgumentParser(
        description="Build an {image_id: patient_id} JSON map from a TMA-style "
                     "patient annotation file, for survival.py's image_patient_map."
    )
    p.add_argument("--annotation-file", required=True, help="Patient-level annotation CSV.")
    p.add_argument("--patient-id-col", default="Patient", help='Patient ID column. Default: "Patient"')
    p.add_argument("--spot-col", default="TMA spot / region",
                    help='Column listing this patient\'s core/region numbers. Default: "TMA spot / region"')
    p.add_argument("--region-sep", default=",", help='Separator within --spot-col. Default: ","')
    p.add_argument("--image-dir", required=True,
                    help="Directory of *_matched_with_boundaries.csv files (triads.py's input_dir).")
    p.add_argument("--file-suffix", default="_matched_with_boundaries.csv",
                    help="Filename suffix stripped to get each image_id. "
                         'Default: "_matched_with_boundaries.csv"')
    p.add_argument("--region-regex", default=r"reg(\d+)",
                    help=r'Regex (first capture group = region number) applied to each image_id. '
                         r'Default: r"reg(\d+)" — matches "..._reg001_...").')
    p.add_argument("--output", default="image_patient_map.json",
                    help="Output JSON path. Default: image_patient_map.json")
    args = p.parse_args()

    if not os.path.exists(args.annotation_file):
        print(f"[build_image_patient_map] ERROR: --annotation-file not found: {args.annotation_file}",
              file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(args.image_dir):
        print(f"[build_image_patient_map] ERROR: --image-dir not found: {args.image_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        image_patient_map = build_map(
            args.annotation_file, args.patient_id_col, args.spot_col, args.region_sep,
            args.image_dir, args.region_regex, args.file_suffix,
        )
    except ValueError as e:
        print(f"[build_image_patient_map] ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    if not image_patient_map:
        print("[build_image_patient_map] ⚠️  No image->patient mappings produced — nothing written.")
        sys.exit(1)

    with open(args.output, "w") as f:
        json.dump(image_patient_map, f, indent=2, sort_keys=True)

    print(f"\n[build_image_patient_map] ✓  Wrote {len(image_patient_map)} entries → {args.output}")
    print(f"[build_image_patient_map]    In your experiment YAML:")
    print(f"[build_image_patient_map]      analysis:")
    print(f"[build_image_patient_map]        survival:")
    print(f"[build_image_patient_map]          image_patient_map: \"{args.output}\"")


if __name__ == "__main__":
    main()
