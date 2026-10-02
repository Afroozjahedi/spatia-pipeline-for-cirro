# QuPath: cell-typing visualization (automatic track)

**Where this fits:** a manual, optional QuPath step for visually reviewing the full cell-typing call -- including Unassigned cells -- after the cell_typing step has run with `cell_typing.qupath_export.enabled: true`.

**Status: not yet run against a real QuPath project.** Written and reviewed, but not tested on an actual slide yet -- check the number of cells it renders against the export step's own printed count (`[celltype_qupath_export] N cells exported...`) the first time you use it.

## What it does

Draws every typed cell as a colored dot directly in your QuPath project -- one color per cell type, assigned automatically and consistently within one TSV. **Unassigned cells are always grey**, so they stand out from confident calls. Click any cell to see its metadata: cell type, tissue/core, experiment group, and -- this is the part built specifically for Unassigned cells -- its full **positive marker list**, plus a **candidate_pattern** suggestion when its marker signature matches one of the top-15 recurring patterns already identified in `cell_typing_data/unassigned_candidate_cell_types.csv`. You can inspect an Unassigned cell's actual marker profile, in its real tissue context, without leaving QuPath to cross-reference a separate CSV.

Set `unassignedOnly = true` at the top of the script to render only Unassigned cells, for a focused adjudication pass.

## How to run it

Open the QuPath project for the specific image/core you want to review, edit `tsvPath` at the top of the script to point at that core's `{tissue_id}_celltypes.tsv` (under `results/celltype_qupath_exports/`), then run it from QuPath's Script Editor. Save the QuPath project afterward if you want to keep what was drawn.

## Things to know

- Make sure the file you point the script at actually matches the image you have open in QuPath -- there's no automatic check for this.
- Cell-type colors are assigned by sorting the distinct cell types found in that one TSV and cycling a fixed 20-color palette -- consistent within one run of the script, not guaranteed to match a different TSV's assignment for the same cell type (same caveat as `spatia/analysis/visualization.py`'s own color map).
- `candidate_pattern` is only ever set for Unassigned cells, and only for the ~80-95% of them covered by the top-15 patterns (see the export step's own printed coverage %) -- a blank value means either the cell is confidently typed, or its pattern was too rare/diffuse to make the top 15.
- This hasn't been run on a real project yet, so treat the first use as a test and compare the rendered count against what the export step reported.

## Related

- `spatia/analysis/cell_typing.py` (`run_cell_typing` / `export_qupath_celltypes`) -- types cells, writes `*_celltypes.tsv`
- `cell_typing_data/unassigned_candidate_cell_types.csv` -- the source-of-truth pattern table `candidate_pattern` is read back from
- `05-5_celltyping_visualization.groovy` -- this script
- `07-2_triad_visualization.groovy` -- sibling script for the triads export, same import/render pattern
