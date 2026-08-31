# SPATIA pipeline — code documentation index

Generated 2026-08-12, trimmed 2026-08-19 to the Nature Protocols manuscript's scope: the SPATIA pipeline itself plus the CRC TMA dataset it validates on. Each doc covers purpose, entry points, key functions, config keys, inputs/outputs, dependencies, and a notes/risks section flagging gaps worth your attention.

**Out of scope, removed from this folder (2026-08-19):** the UI (`ui_app.md`), the separate agentic orchestrator entry point (`run_agentic.md` — the underlying `agentic.py` logic is still covered by `spatia_analysis_agentic.md` below), and all 5 Matusiak-dataset docs (`prepare_matusiak_data.md`, `matusiak_make_figures.md`, `matusiak_radius_sweep.md`, `matusiak_triad_counts_fast.md`, `matusiak_triad_permutation.md`) — the Matusiak/CODEX pan-cancer dataset was explicitly ruled out of scope for this paper (it belongs to a separate pan-cancer triads paper). If that scope ever changes, these were git-untracked to begin with, so nothing was lost from version history — they'd need to be regenerated.

The 8 core pipeline step docs also have a **Flow** section — a Mermaid flowchart of input → key functions → output — right after the entry-point code block, for a faster read than the prose below it. Mermaid code fences render automatically on GitHub. In VS Code they render only with the **Markdown Preview Mermaid Support** extension installed (search that name in the Extensions panel); without it, the diagram shows as a plain code block. Open Preview with `Cmd+Shift+V` either way.

## How the pieces fit together

The repo has three layers:

1. **`spatia/analysis/*.py`** — the actual pipeline logic, one module per stage. Config-driven, no hardcoded paths. Nothing in here has a `__main__` — everything is called by an orchestrator.
2. **Orchestrator** (`run_pipeline.py`) — drives the pipeline modules from a YAML config. (`run_agentic.py` and `ui/app.py` also exist in the repo but are out of scope for this doc set — see the note at the top of this file.)
3. **CRC TMA prep scripts** (`prepare_crc_data.py`, `prepare_crc_for_celltyping.py`, `compare_expert_vs_auto.py`) — dataset-specific glue that feeds the paper's validation dataset into the pipeline and scores cell-typing accuracy against expert labels.

## Pipeline stages (`spatia/analysis/`)

Ordered by pipeline position — raw images in, statistics out:

| # | Step | Doc | Status |
|---|---|---|---|
| 0 | tif_conversion | [spatia_analysis_tif_conversion.md](spatia_analysis_tif_conversion.md) | off by default |
| 1 | roi_masking | [spatia_analysis_roi_masking.md](spatia_analysis_roi_masking.md) | off by default |
| 2 | segmentation | [spatia_analysis_segmentation.md](spatia_analysis_segmentation.md) | off by default |
| 3 | preprocessing | [spatia_analysis_preprocessing.md](spatia_analysis_preprocessing.md) | always on |
| 4 | cell_typing | [spatia_analysis_cell_typing.md](spatia_analysis_cell_typing.md) | always on |
| — | agentic (LLM layer inside cell_typing) | [spatia_analysis_agentic.md](spatia_analysis_agentic.md) | optional, separate entry point |
| 5 | triads | [spatia_analysis_triads.md](spatia_analysis_triads.md) | config-gated |
| 6a | functional | [spatia_analysis_functional.md](spatia_analysis_functional.md) | config-gated |
| 6b | survival | [spatia_analysis_survival.md](spatia_analysis_survival.md) | config-gated |
| — | validation (cross-cutting, runs after every step) | [spatia_analysis_validation.md](spatia_analysis_validation.md) | — |
| — | visualization (not wired into the pipeline yet) | [spatia_analysis_visualization.md](spatia_analysis_visualization.md) | standalone utility |

## Orchestrator

| File | Doc |
|---|---|
| `run_pipeline.py` | [run_pipeline.md](run_pipeline.md) |

## Dataset-specific scripts

**CRC TMA dataset** — the dataset this Nature Protocols paper validates on:

| File | Doc |
|---|---|
| `prepare_crc_data.py` | [prepare_crc_data.md](prepare_crc_data.md) |
| `prepare_crc_for_celltyping.py` | [prepare_crc_for_celltyping.md](prepare_crc_for_celltyping.md) |
| `compare_expert_vs_auto.py` | [compare_expert_vs_auto.md](compare_expert_vs_auto.md) — computes the accuracy metric behind the paper's current validation blocker |

**Dataset-agnostic** — usable with any dataset, not tied to CRC:

| File | Doc |
|---|---|
| `prepare_matched_cells.py` | [prepare_matched_cells.md](prepare_matched_cells.md) — canonical starting point for feeding a new dataset into `triads.py`, generalizes what `prepare_crc_data.py`/`prepare_matusiak_data.py` do per-dataset |
| `export_triads_for_qupath.py` | [export_triads_for_qupath.md](export_triads_for_qupath.md) — converts `triads.py`'s `*_triad_pairs.csv` into QuPath-importable TSVs for `07-2_triad_visualization.groovy` |
| `pool_tissues_for_celltyping.py` | [pool_tissues_for_celltyping.md](pool_tissues_for_celltyping.md) — concatenates `preprocessing.py`'s per-tissue h5ads into one cohort-level h5ad, since nothing else in the pipeline pools across tissues before `cell_typing.py` |
| `build_image_patient_map.py` | [build_image_patient_map.md](build_image_patient_map.md) — generates the `image_id -> patient_id` JSON `survival.py` needs for TMA-style cohorts (multiple cores per patient), from a patient annotation file's spot/region column |

## QuPath / Groovy scripts

Manual and semi-manual steps that sit alongside the Python pipeline — mostly upstream of `roi_masking.py`/`segmentation.py`, or downstream visualization ports. None of these are tracked in git (the repo's `.gitignore` intentionally excludes all Groovy/QuPath files, same as it excludes most of the Python docs).

| File | Doc | Note |
|---|---|---|
| `00_ROI_extract_mask_project.groovy` | [qupath_00_roi_extract_mask_project.md](qupath_00_roi_extract_mask_project.md) | Exists as an identical duplicate in two folders — see the script's own maintenance note |
| `qupath_scripts/00b_auto_tma_dearray.groovy` | [qupath_00b_auto_tma_dearray.md](qupath_00b_auto_tma_dearray.md) | Never run against real data — self-flagged in its own header |
| `05-4-1_cell-typing_segmentation-boundries_visualization.groovy` | [qupath_05-4-1_segmentation_boundaries_visualization.md](qupath_05-4-1_segmentation_boundaries_visualization.md) | Polygon (boundary-accurate) cell rendering |
| `05-4-2_cell-typing_visualization.groovy` | [qupath_05-4-2_celltyping_visualization.md](qupath_05-4-2_celltyping_visualization.md) | Source script `spatia/analysis/visualization.py` ports from |
| `07-2_triad_visualization.groovy` | [qupath_07-2_triad_visualization.md](qupath_07-2_triad_visualization.md) | Exists as an identical duplicate in two folders; never run against real data — self-flagged in its own header. Input TSV comes from `export_triads_for_qupath.py` |

## Cross-file findings worth your attention

These are the issues that came up more than once, or that matter most for the "agentic pipeline" goal stated in this project's instructions. Confidence scores reflect how directly each is observable in the code versus inferred.

**`survival.py` is dataset-agnostic; Cox PH is implemented (confidence: high).** No hardcoded annotation-file column names, image→patient mapping, or group-code map — all config-driven (`patient_id_col`, `image_patient_map`, `group_code_map`). `_cox_ph()` is config-gated behind `analysis.survival.covariates`. Optional per-patient functional-marker exposure covariates via `analysis.survival.marker_exposures`, re-aggregated from the same per-cell files `functional.py` reads. See `spatia_analysis_survival.md`; see git log for the rewrite history and verification detail.

**Validation coverage gap — partially closed (confidence: high).** `validation.py` now has real checks for `segmentation`, `preprocessing`, `cell_typing`, and `triads`. `functional` and `survival` — two of the pipeline's scientific-output steps — are still wired to a passthrough stub that always reports success. `tif_conversion` and `roi_masking` still aren't registered at all. `run_pipeline.py`'s "validate after every step, halt on failure" safety net now covers segmentation's own known failure mode (a per-image `spacec` crash that `segmentation.py` deliberately doesn't raise on — see `spatia_analysis_segmentation.md`), but `functional`/`survival`/`tif_conversion`/`roi_masking` are still open gaps.

**UI/orchestrator step-list drift (confidence: high, out of scope for this paper but still true in the code).** `ui/app.py` hand-copies `ALL_STEPS` as a 5-item list, missing `tif_conversion`/`roi_masking`/`segmentation` from `run_pipeline.py`'s canonical 8-item list. Its own doc was removed from this set since the UI isn't part of the Nature Protocols pipeline — flagging here so it isn't forgotten if the UI becomes relevant again.

**On "agentic pipeline" scope (confidence: medium).** Right now, "agentic" in this codebase means one thing: `agentic.py`'s LLM-assisted cluster labeling (still documented — [spatia_analysis_agentic.md](spatia_analysis_agentic.md)), invoked via the separate `run_agentic.py` entry point (its own doc was removed from this set as out of scope). Step *selection and sequencing* in `run_pipeline.py` is fully static and config-driven — there's no point where an agent decides which steps to run or reacts to intermediate results.

**No cross-tissue pooling before cell_typing — now closed by an optional new script (confidence: high, verified live).** `preprocessing.py` pools cells *within* a tissue (via `extract_tissue_identifier()`), writing one combined h5ad per `tissue_id` — never across tissues. `cell_typing.py`'s `input_file` is a single hardcoded path with no glob/loop, so a multi-tissue cohort (e.g. multiple timepoints or patients) typed by pointing `input_file` at one tissue's file gets GMM thresholds (and, in `semi_automatic` mode, Leiden clusters) fit independently per tissue — any batch drift between tissues gets baked into each tissue's own thresholds rather than calibrated out. `pool_tissues_for_celltyping.py` (dataset-agnostic, optional) now closes this: concatenates every `*_combined_all_experiment_groups.h5ad` in a directory into one cohort-level h5ad for `cell_typing.input_file`, with marker-panel-mismatch handling verified against a synthetic 3-tissue fixture. Whether joint typing is the *right* methodological choice for a given cohort (vs. per-tissue typing, or a batch-correction approach like scVI/scANVI) is still a call the user has to make — this script only removes the "there's no way to do it" blocker, not the judgment call.

**Duplicate, git-untracked Groovy file — confirmed intentional, not a gap (confidence: high, verified byte-for-byte).** `00_ROI_extract_mask_project.groovy` exists identically in both `pipeline with semi-automatic cell typing/` and `pipieline with automatic cell typing/` (note the second folder name is itself a typo). Confirmed deliberate by the user, not accidental drift. Since neither folder is under git, nothing enforces the two copies stay in sync going forward — a maintenance note inside the file itself (both copies) flags that every future edit needs to be applied to both. The file previously had a stray " 2" in its name (an accidental save/duplicate artifact, unrelated to the two-folder duplication) which didn't match how `00b_auto_tma_dearray.groovy` referred to it in its own comments — renamed to drop the " 2" so the actual filename and all cross-references agree.

## Verification

Function signatures, config keys, and file paths cited in each doc were checked against the source during writing (each doc was written directly from a full read of its corresponding file, not from memory or guessing). If the codebase changes, these docs will drift — there's no automated regeneration.
