# Cell-type overlay visualization

**Status:** a standalone tool, not yet wired into the automatic pipeline run — you call it yourself when you want a look.

## What this does

Renders a plot of your cells, colored by cell type and sized by cell area — a quick visual sanity check of your cell-typing results without needing to open QuPath.

## How to use it

Give it a table of cells with `centroid_x`, `centroid_y`, and `cell_type` columns (your cell-typed data, or a triad step's per-cell output), and it produces a PNG. Optionally pass your own color mapping if you want specific colors for specific cell types — otherwise it picks a consistent color per cell type automatically.

## Things to know

- Because this isn't wired into the pipeline yet, it won't run on its own — if you'd like a cell-type overlay generated automatically for every image as part of a normal run, that's a small follow-up worth asking for.
- If you pass a partial custom color mapping that's missing one of your actual cell types, that cell type renders in flat grey rather than raising an error — worth double-checking your color mapping covers every type you expect if the plot looks like it's missing a category.
