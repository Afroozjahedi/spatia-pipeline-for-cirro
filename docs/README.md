# SPATIA pipeline — documentation index

This folder documents the SPATIA pipeline (the CRC TMA dataset it validates on, for the Nature Protocols paper) and the QuPath/manual steps that sit alongside it.

## Pipeline steps, in order

| # | Step | Doc | On by default? |
|---|---|---|---|
| 0 | Format conversion | [spatia_analysis_tif_conversion.md](spatia_analysis_tif_conversion.md) | Off — only needed for QPTIFF/OME-TIFF raw files |
| 1 | ROI masking | [spatia_analysis_roi_masking.md](spatia_analysis_roi_masking.md) | Off — only needed if you're masking ROIs from whole-slide images |
| 2 | Segmentation | [spatia_analysis_segmentation.md](spatia_analysis_segmentation.md) | Off — only needed if you're starting from raw/masked images rather than pre-segmented data |
| 3 | Preprocessing | [spatia_analysis_preprocessing.md](spatia_analysis_preprocessing.md) | Always on |
| 4 | Cell typing | [spatia_analysis_cell_typing.md](spatia_analysis_cell_typing.md) | Always on |
| 5 | Triad detection | [spatia_analysis_triads.md](spatia_analysis_triads.md) | Turned on via config |
| 6a | Functional marker analysis | [spatia_analysis_functional.md](spatia_analysis_functional.md) | Optional, turned on via config |
| 6b | Survival analysis | [spatia_analysis_survival.md](spatia_analysis_survival.md) | Optional, turned on via config |
| — | Image format support (CODEX/OME-TIFF/QPTIFF) | [spatia_analysis_channels.md](spatia_analysis_channels.md) | Used automatically by segmentation |
| — | Cell-type overlay plot | [spatia_analysis_visualization.md](spatia_analysis_visualization.md) | Standalone tool, run manually |

Run all of it with `run_pipeline.py` — see [run_pipeline.md](run_pipeline.md) for usage.

## Helper scripts

| File | Doc | When you need it |
|---|---|---|
| `build_image_patient_map.py` | [build_image_patient_map.md](build_image_patient_map.md) | Your cohort has multiple images per patient (e.g. several TMA cores) and you're running survival analysis |
| `prepare_matched_cells.py` | [prepare_matched_cells.md](prepare_matched_cells.md) | Feeding a new dataset's cell-typed data into triad detection for the first time |
| `pool_tissues_for_celltyping.py` | [pool_tissues_for_celltyping.md](pool_tissues_for_celltyping.md) | Your cohort has multiple tissues/cores and you want them cell-typed together rather than one at a time (recommended for TMA cohorts, so thresholds are calibrated consistently across cores) |
| `export_triads_for_qupath.py` | [export_triads_for_qupath.md](export_triads_for_qupath.md) | You want to visualize detected triads back in QuPath |

## QuPath / manual steps

These sit alongside the Python pipeline, mostly before ROI masking/segmentation or for visualization:

| File | Doc |
|---|---|
| `00_ROI_extract_mask_project.groovy` | [qupath_00_roi_extract_mask_project.md](qupath_00_roi_extract_mask_project.md) |
| `00b_auto_tma_dearray.groovy` | [qupath_00b_auto_tma_dearray.md](qupath_00b_auto_tma_dearray.md) |
| `05-4-1_...segmentation_boundaries_visualization.groovy` | [qupath_05-4-1_segmentation_boundaries_visualization.md](qupath_05-4-1_segmentation_boundaries_visualization.md) |
| `05-4-2_cell-typing_visualization.groovy` | [qupath_05-4-2_celltyping_visualization.md](qupath_05-4-2_celltyping_visualization.md) |
| `07-2_triad_visualization.groovy` | [qupath_07-2_triad_visualization.md](qupath_07-2_triad_visualization.md) |

## Things worth knowing across the whole pipeline

- **If a step fails partway through a large batch, check the run's error summary before assuming everything is broken** — most steps log and skip a single bad file rather than stopping the whole run, so a lower-than-expected cell/image count is the usual sign something needs a second look, not a hard crash.
- **Automatic pass/fail checking after each step doesn't yet cover every step** — segmentation, preprocessing, cell typing, and triads are checked automatically; functional, survival, format conversion, and ROI masking aren't yet, so it's worth reviewing those steps' own output/log summaries yourself rather than assuming a clean pipeline run means everything downstream is correct too.
- **If your cohort has multiple tissues or cores, pool them before cell typing** (`pool_tissues_for_celltyping.py`) rather than typing each one separately — otherwise each tissue gets its own independently calibrated thresholds, which can introduce apparent differences between tissues that are really just threshold noise. Whether pooling (one shared threshold per marker) or per-tissue typing is the right call for your cohort is a real methodological choice, not just a technical default — worth being deliberate about it for the paper.
