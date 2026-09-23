# Cell-type overlay visualization

**Status (updated 2026-09-23): superseded — do not wire this module in.** The automatic, per-run version of what this module does now lives directly in `cell_typing.py` (`plot_spatial_celltype_overlay()`, folded in 2026-09-12), and already runs for every cell-typing step with no config needed — see [Cell typing](spatia_analysis_cell_typing.md). It also improves on this module for CRC specifically: an auto-generated 40-color palette instead of a fixed color map, since CRC's 23+ cell types don't fit this module's original LILRB2-study palette. This file (`spatia/analysis/visualization.py`) is kept only as a standalone reference/manual tool — it is not imported anywhere in the pipeline and is not the code path that produces your automatic overlays.

## What this does

Renders a plot of your cells, colored by cell type and sized by cell area — a quick visual sanity check of your cell-typing results without needing to open QuPath.

## How to use it

Give it a table of cells with `centroid_x`, `centroid_y`, and `cell_type` columns (your cell-typed data, or a triad step's per-cell output), and it produces a PNG. Optionally pass your own color mapping if you want specific colors for specific cell types — otherwise it picks a consistent color per cell type automatically.

## Things to know

- This module is not wired into the pipeline, and — since `cell_typing.py` already does the equivalent automatically — it shouldn't be. Use it standalone only for a one-off plot with a custom color mapping or variable dot sizing by cell area, neither of which the automatic version does.
- If you pass a partial custom color mapping that's missing one of your actual cell types, that cell type renders in flat grey rather than raising an error — worth double-checking your color mapping covers every type you expect if the plot looks like it's missing a category.
