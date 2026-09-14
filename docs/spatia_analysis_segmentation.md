# Cell segmentation

**Pipeline position:** Step 2 — after ROI masking, before preprocessing.

## What this step does

Finds individual cell boundaries in your masked ROI images and produces one row per cell with its size, shape, position, and average marker intensity — the table the rest of the pipeline is built on. It also saves overlay images so you can visually check that the segmentation actually looks right before trusting the numbers.

## Configuring it

- `paths.masked_roi_dir` — input images (from the ROI masking step)
- `paths.channel_file` — your `channelnames.txt`, listing channels in the order they appear in your images
- `paths.segmentation_results_dir` — where results are written
- `segmentation.seg_method` — `"mesmer"` (default) or `"cellpose"`
- `segmentation.nuclei_channel` — which channel is your nuclear stain; `"auto"` picks it automatically for CODEX-style acquisitions
- `segmentation.membrane_channel_list` — which channel(s) mark the cell boundary/membrane
- `segmentation.compartment` — `"whole-cell"` (default) or `"nucleus"` only
- `segmentation.resize_factor`, `.size_cutoff` — fine-tuning knobs; the defaults work for most data
- `segmentation.generate_overlays` — whether to also save QC overlay images (recommended, on by default)
- `segmentation.channel_check` — verifies each image's channels are in the order you expect before segmenting it (recommended, on by default) — this is what catches an image from a different scan batch or panel revision before it gets silently mislabeled

## What you'll get

- One segmentation result and one CSV per masked ROI image, with per-cell size/shape/position plus average marker intensity
- Overlay PNGs so you can spot-check segmentation quality
- A run summary reporting how many images succeeded and how many failed, and why

## Things to know

- If an image's channels aren't in the order your `channel_file` expects, that image is skipped rather than segmented with the wrong labels, and the run summary lists every image this happened to. This is usually a sign that file came from a different scan batch, panel revision, or export than the rest — worth checking those specific files rather than assuming the check is wrong.
- A one-off failure on a single image (a corrupted file, an out-of-memory error) is logged and the run continues rather than stopping the whole batch — check the run summary's failed-file list afterward. If too many images end up missing, the next pipeline step will flag it when it can't find enough segmentation results to continue.
- If you've changed `resize_factor` away from `1` and a channel's average intensity looks like it's missing from the output CSV entirely (not even a blank value, just absent), that combination is worth a closer look — it's a known edge case rather than a random glitch.
