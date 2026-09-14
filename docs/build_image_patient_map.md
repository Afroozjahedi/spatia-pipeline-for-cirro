# Building an image-to-patient map (TMA cohorts)

## What this is for

If your cohort has multiple images per patient (e.g. several TMA cores from the same patient), survival analysis needs to know which images belong to which patient. This script builds that mapping automatically from your annotation file, instead of typing out potentially 100+ entries by hand.

## Usage

```bash
python build_image_patient_map.py \
    --annotation-file /path/to/patient_with_tls_class.csv \
    --patient-id-col Patient \
    --spot-col "TMA spot / region" \
    --image-dir /path/to/matched_cells_dir \
    --output image_patient_map.json
```

Then point survival analysis at the result:

```yaml
analysis:
  survival:
    image_patient_map: "image_patient_map.json"
```

## Arguments

- `--annotation-file` (required) — your patient-level annotation spreadsheet
- `--patient-id-col` (default `"Patient"`)
- `--spot-col` (default `"TMA spot / region"`) — the column listing which core/region numbers belong to each patient
- `--region-sep` (default `","`) — how those numbers are separated within that column
- `--image-dir` (required) — your folder of per-image cell data files
- `--region-regex` — how to pull the region number out of each image's filename (default matches `"reg"` followed by digits); change this if your filenames use a different pattern
- `--output` (default `image_patient_map.json`)

## What you'll get

A JSON file mapping each image to a patient, plus a report of anything that didn't match cleanly: images that couldn't be matched to any patient, patients with no matching images, and any region number claimed by more than one patient (a sign of a data entry issue worth checking in your annotation file).

## Things to know

- This only matches on the core/region number, not the whole image filename — if your filenames also encode something else your spot list doesn't capture (e.g. the same region number reused across different slides), double-check the output JSON the first time you use this on a new dataset shape.
- Always look at the match report before trusting the output — a real mismatch (unmatched images, or a region claimed by two patients) usually means something in the annotation file needs fixing, not something wrong with the script.
