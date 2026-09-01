# `spatia/analysis/visualization.py`

**Pipeline position: NOT YET WIRED IN** — this module exists and is functional but is not called from `run_pipeline.py` as a step (stated explicitly in the module's own docstring). It's a standalone utility, ready to be adopted.

## Purpose

Python port of a QuPath Groovy script (`05-4-2_cell-typing_visualization.groovy`) that renders a cell-type overlay: each cell plotted as a colored dot at its centroid, sized by cell area, colored by cell type. This replaces a manual QuPath visualization step with a Python one that can run headlessly as part of the pipeline.

The module's docstring is unusually detailed about provenance and scope decisions — worth reading directly, but summarized:
- Ports only the **centroid + color-map** logic from the groovy script's option (b), not the full polygon cell-boundary outlines from a different groovy script (05-4-1). This was an explicit scope decision (dated 2026-07-21 / 2026-07-23 in the docstring).
- The original groovy script's cell-type→color map (`LEGACY_LILRB2_COLOR_MAP`) was built for the LILRB2 mouse triad study's cell-type vocabulary and does **not** match the CRC TMA study's vocabulary (`cell_type_definitions/crc_tma.yaml`) — porting it as-is would silently mislabel/gray out every CRC cell type. Rather than inventing a new hardcoded palette, this module defaults to an automatic categorical palette (reusing `triads.py`'s existing `tab10`-style color-cycling approach) and accepts an optional `color_map` override.

## Entry points (no config-driven pipeline wrapper yet)

```python
from spatia.analysis.visualization import plot_celltype_overlay, build_color_map, print_celltype_breakdown

ax = plot_celltype_overlay(df, output_path="overlay.png")
```

## Key functions

| Function | Role |
|---|---|
| `build_color_map(cell_types, palette=None)` | Assigns each (sorted, for reproducibility) cell type a color by cycling a palette — the same `tab10`-derived 10-color palette `triads.py` uses. |
| `compute_radius(area)` | Port of the groovy script's radius formula: `max(sqrt(area/π), 3.0)` if area present, else a fixed `5.0`. Distinct from (but consistent in spirit with) `preprocessing.py`'s separate QC-TSV radius formula. |
| `plot_celltype_overlay(df, output_path=None, color_map=None, ax=None, title=None, figsize=(12,12), dpi=200)` | Main rendering function. Requires `centroid_x`, `centroid_y`, `cell_type` columns; optionally uses `area` for per-cell dot sizing. Coordinates are assumed already global (no offset re-application needed, matching this pipeline's existing cell-typed output convention). Y-axis inverted to match image-coordinate convention used elsewhere (`triads.py`'s QC panels). |
| `print_celltype_breakdown(df)` | Console summary of cell-type counts/percentages, matching the groovy script's log output format — useful for eyeballing parity against the old QuPath-based results. |

## Inputs / Outputs

- **In:** any DataFrame with `centroid_x`, `centroid_y`, `cell_type` (e.g. the cell-typed h5ad's `.obs`, or a triads step's flagged-cells CSV)
- **Out:** a single overlay PNG (via `output_path`) or a matplotlib `Axes` for embedding in a larger figure

## Dependencies

`matplotlib`, `numpy`, `pandas`, stdlib `math`.

## Notes / risks

- **Genuinely not wired into the pipeline yet — this is a known, self-declared gap, not a hidden one (confidence: high, stated directly in the docstring).** The docstring explicitly says "whether it should become its own pipeline step vs. an option on an existing step is an open sub-question ... not decided here." This is exactly the kind of decision point worth resolving deliberately: does this become a `run_pipeline.py` step (e.g. after `cell_typing` or after `triads`), or an optional flag on an existing step's plotting? Left as-is, it's a working function nobody's calling automatically.
- **`LEGACY_LILRB2_COLOR_MAP` is dead reference code, kept intentionally (confidence: high — explicitly labeled "reference only" in the code).** Not a risk so much as a maintenance note: if the LILRB2 project's color conventions ever change, this stale copy won't be updated in lockstep, and anyone copy-pasting it for "manual reuse" (as the docstring suggests) should verify it's still current.
- **Unmapped cell types fall back to a fixed grey, not an error (confidence: medium).** If a caller passes a partial `color_map` that's missing an actual cell type present in `df`, that type silently renders as `#646464` (the groovy script's "Unknown" grey) rather than raising — consistent with the original script's behavior, but could mask a typo'd color-map key.
