# QuPath: preprocessing QC visualization

**Where this fits:** a manual, optional QuPath step for visually reviewing which cells were kept/excluded by preprocessing's QC filters, after the preprocessing step has run with `preprocessing.qupath_export.enabled: true`.

**Status: Afrouz's own script, carried over from the manual notebook pipeline (originally `notebook/04-2_preprocessing_visualization.groovy`) -- previously used against real data there, now wired into the automated pipeline's TSV output for the first time.** One real mismatch was found and fixed on the Python side (see below) before this TSV was trustworthy as input to it; the script itself was not rewritten.

## What it does

Draws every segmented cell as a colored dot: green for cells that passed QC (`Included`), and a distinct color per exclusion reason (`Excl_SmallArea` red, `Excl_LowDAPI` orange, `Excl_SmallArea_LowDAPI` purple, `Excl_Noise` blue). Lets you visually sanity-check the preprocessing QC thresholds against the actual tissue, not just the aggregate counts in `qupath_export_summary.csv`.

Also corrects local crop coordinates to whole-slide coordinates when an `offset_x`/`offset_y` column or an `_x{X}_y{Y}_w{W}_h{H}` filename pattern is present. For the current CRC TMA pipeline, sample/image IDs (e.g. `reg003_X01_Y01_Z06`) don't carry pixel offsets at all, so this will fall back to offset `(0, 0)` -- correct as long as you're viewing each core/region as its own QuPath image (which is how this pipeline treats them), not plotting multiple cores back onto one original whole-slide scan.

## How to run it

Open the QuPath project for the specific image/core you want to review, edit `tsvPath` at the top of the script to point at that core's `{image_id}_qupath_cells.tsv` (under `results/qupath_exports/`), then run it from QuPath's Script Editor. Save the QuPath project afterward if you want to keep what was drawn.

## Things to know -- fixed before first use

- **The DAPI column name mismatch.** This script looks up a column literally named `DAPI`. The pipeline's actual nuclei-intensity column is named after whatever the panel resolves to (`nuclei_channel: "auto"` -> `panel_names[0]` -- for the CRC TMA panel this is `"HOECHST1 (C1)"`, confirmed against `experiments/crc_tma_full/channelNames.txt`), never the literal string `"DAPI"`. Unguarded, this would have made Groovy's negative-index behavior (`cols[-1]`) silently read the *last* column in the TSV (`dapi_threshold`, one scalar repeated for every cell) as if it were each cell's real DAPI value -- wrong data, no error. **Fixed on the Python side**, not by editing this script: `preprocessing.py`'s `run_qupath_export()` now also writes a literal `DAPI` alias column (a copy of the real resolved nuclei column) specifically so this script works against the TSV unmodified.
- Offset correction defaults to `(0, 0)` for this pipeline's current filenames (see above) -- treat that as expected, not a bug, unless you're deliberately viewing cores plotted back onto an original whole-slide scan.
- This is the first real run of this script against the automated pipeline's TSV output (as opposed to the manual notebook pipeline it came from) -- compare the imported cell count and classification breakdown against the export step's own printed summary the first time.

## Related

- `spatia/analysis/preprocessing.py` (`run_preprocessing` / `run_qupath_export`) -- classifies cells, writes `*_qupath_cells.tsv` (now including the `DAPI` alias column)
- `04-2_preprocessing_visualization.groovy` -- this script
- `05-5_celltyping_visualization.groovy` / `07-2_triad_visualization.groovy` -- sibling scripts, same import/render pattern, for cell typing and triads respectively
