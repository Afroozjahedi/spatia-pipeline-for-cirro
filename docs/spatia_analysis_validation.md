# `spatia/analysis/validation.py`

**Pipeline position:** Cross-cutting — called by `run_pipeline.py` after each step, not itself a pipeline stage.

## Purpose

Post-step output validation. Each validator checks that a step's expected outputs (1) exist on disk, (2) are non-empty, (3) can be opened without error (for h5ad files), and (4) contain expected content (obs columns present, cell counts > 0, plausible cell-type diversity). Designed to catch a step that "succeeded" (no exception) but produced garbage or incomplete output before the next step wastes time consuming it.

## Entry points

```python
from spatia.analysis.validation import validate_step, validate_all

ok, errors = validate_step("preprocessing", cfg)   # single step
all_passed, results = validate_all(cfg)             # every registered step
```

## Key functions

| Function | Role |
|---|---|
| `ValidationError` | Small class (`step`, `message`, `path`) with a readable `__str__`. Not an exception — collected into lists and returned, not raised. |
| `_file_exists_and_nonempty(path)` | Existence + 0-byte check. |
| `_open_h5ad(path)` | Tries `anndata.read_h5ad`; returns success flag, error message, and the loaded object (or `None`). |
| `_check_obs_columns(adata, required_cols, path)` | Returns which of `required_cols` are missing from `adata.obs`. |
| `validate_segmentation(cfg)` | Checks ≥1 `*_mesmer_result.csv` exists under `segmentation_results_dir`, that each is non-empty with the columns `preprocessing.py` reads (`label`, `x`, `y`, `area`), and — the important one — compares the count of masked ROI TIFFs found in `masked_roi_dir` against the count of CSVs actually produced. A shortfall means some image(s) silently failed segmentation or were skipped for a channel mismatch; since this validator only receives `cfg` (not `run_cell_segmentation`'s in-memory `errors` list), it re-derives the same signal from disk instead. |
| `validate_preprocessing(cfg)` | Checks ≥1 `*_combined_all_experiment_groups.h5ad` exists, opens cleanly, has a `experiment_group` column, has >0 cells; and that every group in `experiment.groups` appears somewhere **across all combined h5ad files together** (fixed 2026-09-09 to check at the whole-run level instead of per file — see Notes/risks). |
| `validate_cell_typing(cfg)` | Branches on mode: for `semi_automatic` with `cluster_labels_file` still `null`, checks only that the clustered h5ad exists (an intentional early stop, not a failure) and prints next-step guidance. Otherwise checks the cell-typed h5ad exists, has `cell_type`/`experiment_group` columns, ≥2 distinct cell types, and flags if >50% of cells are `"Unassigned"`. |
| `validate_triads(cfg)` | Checks `triad_summary.csv` and `experiment_group_comparison_counts.csv` exist and are non-empty, that `triad_summary.csv` has the expected columns and ≥1 row, and that at least one `*_triad_pairs.csv` exists. |
| `_validate_passthrough(cfg)` | Stub validator (always passes) — currently registered for `functional` and `survival`, which have no dedicated checks yet. |
| `validate_step(step, cfg)` | Dispatches to the right validator by name; raises `ValueError` for an unregistered step name. |
| `validate_all(cfg, steps=None)` | Runs every (or a specified subset of) registered validators, returns a combined pass/fail plus per-step results dict. |

## Registered steps

`segmentation`, `preprocessing`, `cell_typing`, `triads`, `functional` (passthrough), `survival` (passthrough).

## Dependencies

`anndata` (lazily imported inside `_open_h5ad`), `pandas` (lazily imported inside `validate_triads` and `validate_segmentation`), stdlib only otherwise.

## Notes / risks

- **Fixed 2026-09-09 — `validate_preprocessing()`'s per-file experiment_group check produced a false "FAILED" on a fully-successful run (confidence: high, reproduced against a real 138/138, 0-error run).** The check previously required *each individual* combined h5ad (one per tissue) to contain both configured `experiment.groups`. That's only true when a dataset scans every experiment_group as replicate stains of the same physical tissue. This CRC TMA dataset (and any dataset using `image_experiment_group_map` to assign whole slide folders to a single group — see `spatia_analysis_preprocessing.md` / the yaml's `experiment.image_experiment_group_map`) scans CLR and DII as entirely separate slide folders, so **every** tissue's combined h5ad legitimately contains only one group — the old check failed on every single tissue, every run, regardless of correctness, and `run_pipeline.py` reported "preprocessing FAILED" even when `run_preprocessing()` itself succeeded completely. Fixed by aggregating `experiment_group` values across *all* combined h5ad files first, then checking `experiment.groups` is a subset of that aggregate — correct for both the split-slide-folder case (this dataset) and the paired-replicate case (where the old per-file check would also have passed, since it's subsumed by the aggregate check). Audited the rest of the codebase (`functional.py`, `roi_masking.py`, `survival.py`, `triads.py`, the other `validate_*` functions in this file) for the same per-item-must-have-every-group anti-pattern — found no other occurrences (confidence: 90%, grep + targeted read, not an exhaustive line-by-line audit).
- **`segmentation` is registered (confidence: high).** `validate_segmentation()` catches a gap `segmentation.py` deliberately leaves open by design: `run_cell_segmentation()` collects per-image failures (generic `spacec` exceptions, and channel mismatches) but doesn't raise on the generic-failure case, so one bad image doesn't halt a whole batch. `validate_segmentation()` turns "some images silently produced no output" into a pipeline-halting failure, by comparing masked-ROI-TIFF count against produced-CSV count on disk — see `spatia_analysis_segmentation.md` for the full story. `tif_conversion` and `roi_masking` remain unregistered (see below).
- **`functional` and `survival` have no real validation — silent gap (confidence: high, directly observable).** Both are wired to `_validate_passthrough`, which always returns `(True, [])`. So a `functional` or `survival` step that produced an empty `functional_marker_summary.csv` or crashed partway through KM plotting (if the exception were ever swallowed upstream) would report as "validated OK." Given these are two of the pipeline's scientific-output steps, this is probably the single highest-value gap to close next — even a basic "does `functional_marker_summary.csv` exist and have >0 rows" check would catch most failure modes.
- **`tif_conversion` and `roi_masking` still aren't registered (confidence: high, directly observable — narrower gap than before now that `segmentation` is covered).** `_VALIDATORS` covers `segmentation` onward, but not the two steps upstream of it. Both of those already return an `errors`/mismatch-style structure from their own run functions (`tif_conversion.py`'s `run_tif_conversion` returns `{converted_files, errors, output_dir}`), so a validator for either would follow the same disk-based pattern `validate_segmentation()` just established — worth doing next if these two steps see real use.
- **Errors are collected, not raised — callers must remember to check the boolean (confidence: medium).** This is a reasonable design (lets `run_pipeline.py` decide whether to halt or warn-and-continue), but it does mean any code path that calls a validator and forgets to check the returned `bool` will silently treat a failed validation as a no-op. Worth confirming `run_pipeline.py` always checks and acts on the result — see that script's own doc.
