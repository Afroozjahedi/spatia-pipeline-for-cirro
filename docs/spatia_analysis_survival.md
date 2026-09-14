# Survival analysis

**Pipeline position:** Step 6b (optional) — runs after triads (and after functional, if you're using marker exposures).

## What this step does

Links your triad findings — and, optionally, functional marker levels — to patient outcomes: Kaplan-Meier survival curves with a log-rank test, and optionally a multivariate Cox regression. Works with any cohort's own outcome/annotation file; nothing about column names or group codes is hardcoded to one dataset.

## Image-to-patient mapping

By default each image is treated as one patient, which is fine for a cohort with one image per subject. If your cohort has multiple images per patient (e.g. several TMA cores from the same patient), set `image_patient_map` to tell the pipeline which images belong to which patient — `build_image_patient_map.py` can generate this file for you from a TMA annotation sheet.

## Configuring it

- `experiment.name`, `.groups`
- `paths.output_dir`, `.input_dir`
- `analysis.survival.enabled` — must be `true`
- `analysis.survival.patient_annotation_file` (required) — your external clinical/outcome data; this can't be derived from the imaging data itself
- `analysis.survival.patient_id_col` — which column in that file identifies each patient
- `analysis.survival.image_patient_map` — only needed for multi-image-per-patient cohorts (see above)
- `analysis.survival.group_col` / `.group_code_map` — to split survival curves by experiment group
- `analysis.survival.area_per_image_um2` — a flat per-image tissue-area constant used for triad density (see caveat below)
- `analysis.survival.split_by` — `"median"` or a specific numeric threshold, for splitting patients into high/low groups
- `analysis.survival.outcomes` — which outcomes to analyze (e.g. overall survival, disease-free survival), each with its own duration/event columns
- `analysis.survival.marker_exposures` — optional; lets you also test whether a patient's average marker level (e.g. within-triad exhaustion marker expression) predicts outcome, alongside triad density
- `analysis.survival.covariates` — set this to also run a multivariate Cox regression alongside the simpler high/low split analysis

## What you'll get

- `patient_cohort_summary.csv` — one row per patient with triad density, any marker exposures, and your annotation data merged together
- Kaplan-Meier plots and `survival_logrank_results.csv` for every outcome you configured
- Cox regression results, if you set `covariates`

## Things to know

- Unlike triad detection, this step always uses one flat tissue-area number per image rather than each image's own measured area — for a cohort where image size varies a lot, this can bias density-based results for patients whose images are unusually small or large.
- Cox regression only runs if you set `covariates`, and it skips itself (printing why) rather than fitting on too little data, if you don't have enough patients relative to the number of covariates requested.
- If you configure several marker exposures alongside multiple outcomes, keep in mind the number of statistical tests grows with each one, and this step doesn't currently correct for that the way the functional-marker step does — worth thinking about how many you run at once for the paper.
