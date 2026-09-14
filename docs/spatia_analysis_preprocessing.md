# Preprocessing

**Pipeline position:** Step 3 — after segmentation, before cell typing (always runs).

## What this step does

Cleans up the raw per-cell measurements coming out of segmentation before anything scientific is done with them: filters out cells that are too small or too dim to be real (likely segmentation artifacts), normalizes marker intensities (z-score) so different images and batches are comparable, and automatically detects and removes background noise per marker. Produces both per-tissue/group files and one combined file across everything, plus visualizations for sanity-checking the filtering, and a QC summary report generated automatically at the end.

## Configuring it

- `experiment.groups`
- `paths.segmentation_results_dir` — input (from segmentation)
- `paths.output_dir`
- `preprocessing.nuclei_channel` (default `"auto"`) — which column is your nuclear stain; `"auto"` picks it for you
- `preprocessing.last_marker` — where your marker columns end and metadata columns begin; `"auto"` resolves it from your panel file, or specify the actual last marker name yourself
- `preprocessing.qc_filter.size_percentile` / `.dapi_percentile` (both default `1`) — how aggressively to filter out too-small or too-dim cells
- `preprocessing.noise.cut_off` / `.count_bin` — tuning knobs for automatic noise removal; the defaults work for most data

## What you'll get

- Per-tissue and combined processed data files (CSV + h5ad), ready for cell typing
- Marker visualization images to check filtering and normalization worked as expected
- A QC summary report (tissue/group breakdown, why cells were removed, marker-level QC)
- Optionally, a QuPath-importable export classifying every original cell as kept or excluded (and why), for visual review

## Things to know

- If `nuclei_channel` or `last_marker` don't clearly match a column in your data, the run tells you which column it picked instead — worth reading that message carefully the first time you run a new panel, since a wrong guess here silently shifts which columns are treated as markers vs. metadata.
- A single image failing to process (a bad file, an unexpected format) is logged and skipped rather than stopping the whole batch — check the run's error summary if your final cell count looks lower than expected.
- The QC summary report is generated in a way that can never fail the actual preprocessing run — if the report itself looks wrong or incomplete, that's independent of whether your processed data is fine.
