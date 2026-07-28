"""
spatia.analysis.visualization
==============================

Python port of the QuPath cell-type overlay script
(``pipeline with semi-automatic cell typing/05-4-2_cell-typing_visualization.groovy``),
per Decision Q7 (2026-07-21): Afrouz chose option (b) — port the QuPath
overlay logic into Python rather than keep QuPath as a manual step — and,
per the Day 5 (continued, part 2) scope correction, specifically the
**centroid + color-map** logic from 05-4-2 only (not 05-4-1's full polygon
cell-boundary outlines).

Status: Q7 groundwork (2026-07-23 / Day 7). Logic below was extracted by
reading 05-4-2 in full (not written from memory of what it "probably does").
Verified against the real cell dataframe schema (`centroid_x`, `centroid_y`,
`cell_type`, `area`, `cell_id`) already produced by
``spatia/analysis/cell_typing.py`` and already consumed the same way in
``spatia/analysis/triads.py`` (see its `axes[0].scatter` cell-type panel).
The radius formula below also already has Python precedent elsewhere in this
codebase: ``preprocessing.py``'s QC-TSV export block uses the identical
``sqrt(area / pi)`` clipped to a floor for its own (non-cell-type) QuPath
TSV — this port reuses that same math for consistency, not a new formula.

NOT YET WIRED IN: this module is not yet called from `run_pipeline.py` as a
step, and whether it should become its own pipeline step vs. an option on
an existing step is an open sub-question from Day 5's Q7 scope note — not
decided here (operating rule 1).

IMPORTANT — cell-type/color mapping is NOT reusable verbatim:
The 05-4-2 groovy script's ``cellTypeColors`` map is hardcoded for the
LILRB2 mouse triad study's cell types (e.g. "CD8 T cells", "CD4 Tregs",
"Monocytes/Macs"). The CRC TMA cell types actually produced by
``cell_type_definitions/crc_tma.yaml`` are a completely different,
CRC-specific vocabulary (e.g. "tumor cells", "CD68+ macrophages GzmB+",
"CD163+ macrophages"), with no overlap. Porting the LILRB2 map as-is would
silently mislabel/gray-out every CRC cell type. Rather than inventing a
new hardcoded 24-color CRC map myself (an analytical/design choice that is
Afrouz's call), this port defaults to an automatic categorical palette
(the same `tab10`-style cycling approach `triads.py` already uses for its
own cell-type panel), and accepts an optional explicit `color_map` override
for when a specific palette is wanted. The original LILRB2 map is kept
below only as `LEGACY_LILRB2_COLOR_MAP`, for reference / manual reuse on
that other (pan-cancer/mouse) project — not applied by default here.
"""

from __future__ import annotations

import math
import os
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Reference only — the ORIGINAL 05-4-2 groovy color map, for the LILRB2
# mouse triad study's cell-type vocabulary. NOT used as the default here
# because CRC TMA cell types (crc_tma.yaml) don't match these names.
# Kept for traceability back to the source script and manual reuse if ever
# needed on the other (pan-cancer) project.
# ---------------------------------------------------------------------------
LEGACY_LILRB2_COLOR_MAP: Dict[str, str] = {
    "CD8 T cells":     "#1f77b4",  # blue        (31, 119, 180)
    "CD4 T cells":     "#2ca044",  # green       (44, 160, 44)
    "CD4 Tregs":       "#bcbd22",  # yellow-green(188, 189, 34)
    "Macrophages":     "#ff7f0e",  # orange      (255, 127, 14)
    "Monocytes":       "#9467bd",  # purple      (148, 103, 189)
    "Monocytes/Macs":  "#8c564b",  # brown       (140, 86, 75)
    "Mono/Mac":        "#8c564b",  # brown       (140, 86, 75)
    "Neutrophils":     "#e377c2",  # pink        (227, 119, 194)
    "Dendritic cells": "#d62728",  # red         (214, 39, 40)
    "B Cells":         "#17becf",  # teal        (23, 190, 207)
    "NK":              "#ffbb78",  # light orange(255, 187, 120)
    "NK T cells":      "#aec7e8",  # light blue  (174, 199, 232)
    "Tregs":           "#ff9896",  # light red   (255, 152, 150)
    "Artifact/CD45-":  "#969696",  # grey        (150, 150, 150)
    "Mix":             "#c49c94",  # beige       (196, 156, 148)
    "Unknown":         "#646464",  # dark grey   (100, 100, 100)
}

# Same tab10-derived cycling palette triads.py already uses for its
# cell-type panel — reused here for visual consistency across the pipeline's
# plots rather than introducing a second, different default palette.
_DEFAULT_PALETTE = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
    "#a65628", "#f781bf", "#999999", "#00ced1", "#8b0000",
]


def build_color_map(cell_types, palette=None) -> Dict[str, str]:
    """Assign each cell type a color by cycling a palette, sorted for
    deterministic (rerun-stable) assignment. Mirrors triads.py's existing
    `COLOR_PALETTE[i % len(COLOR_PALETTE)]` approach."""
    palette = palette or _DEFAULT_PALETTE
    return {
        ct: palette[i % len(palette)]
        for i, ct in enumerate(sorted(set(cell_types)))
    }


def compute_radius(area: Optional[float]) -> float:
    """Port of 05-4-2's ELLIPSE RADIUS CALCULATION:
        area available -> max(sqrt(area / pi), 3.0)
        area missing   -> 5.0
    (Distinct from preprocessing.py's QC-TSV radius, which always assumes
    area is present; this keeps 05-4-2's explicit missing-area fallback.)
    """
    if area is None or (isinstance(area, float) and math.isnan(area)):
        return 5.0
    return max(math.sqrt(area / math.pi), 3.0)


def plot_celltype_overlay(
    df: pd.DataFrame,
    output_path: Optional[str] = None,
    color_map: Optional[Dict[str, str]] = None,
    ax: Optional["plt.Axes"] = None,
    title: Optional[str] = None,
    figsize=(12, 12),
    dpi: int = 200,
) -> "plt.Axes":
    """Python port of 05-4-2_cell-typing_visualization.groovy's centroid +
    color-map rendering (option (b), centroid dots — not 05-4-1's full
    polygon cell boundaries, per Afrouz's Q7 scope correction).

    Required columns: centroid_x, centroid_y, cell_type
    Optional column:  area (drives per-cell marker size via compute_radius;
                       falls back to a fixed radius of 5.0 if absent, exactly
                       as the groovy script does for missing area).

    Coordinate handling matches the groovy script's documented behavior:
    centroid_x / centroid_y are assumed to already be GLOBAL coordinates
    (this pipeline's cell_typing.py output, like the original TSV, does not
    need offsets re-applied here).
    """
    required = {"centroid_x", "centroid_y", "cell_type"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"plot_celltype_overlay: missing required columns: {sorted(missing)}")

    cmap = color_map or build_color_map(df["cell_type"].unique())

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    has_area = "area" in df.columns
    for ct, sub in df.groupby("cell_type"):
        color = cmap.get(ct, "#646464")  # unmapped types fall back to the
                                          # groovy script's own "Unknown" grey
        if has_area:
            sizes = sub["area"].apply(compute_radius).to_numpy() ** 2 * math.pi
        else:
            sizes = compute_radius(None) ** 2 * math.pi  # scalar -> broadcasts
        ax.scatter(
            sub["centroid_x"], sub["centroid_y"],
            s=sizes, color=color, linewidths=0, alpha=0.85,
            label=f"{ct} ({len(sub):,})",
        )

    ax.invert_yaxis()  # image-coordinate convention, matches triads.py panels
    ax.set_xlabel("X (px)")
    ax.set_ylabel("Y (px)")
    if title:
        ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, markerscale=0.6, ncol=2, framealpha=0.9,
              loc="upper right", borderpad=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if output_path:
        ax.figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
        print(f"Saved cell-type overlay: {output_path}")

    return ax


def print_celltype_breakdown(df: pd.DataFrame) -> None:
    """Console summary matching 05-4-2's trailing 'Cell-type breakdown' log
    output, for parity when eyeballing results against the old QuPath run."""
    counts = df["cell_type"].value_counts()
    total = len(df)
    print("Cell-type breakdown:")
    for ct, cnt in counts.items():
        print(f"  {ct:<30} : {cnt:>6}  ({cnt * 100.0 / total:4.1f}%)")
