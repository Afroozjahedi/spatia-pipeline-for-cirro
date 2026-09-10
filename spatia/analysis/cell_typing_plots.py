"""
spatia/analysis/cell_typing_plots.py
=====================================
Diagnostic/reporting plot suite for a completed run_cell_typing() call.

Added 2026-09-10 at Afrouz's explicit request, after the first real
end-to-end cell_typing run on the pooled CRC TMA cohort (300,008 cells,
23 cell-type definitions) produced only the two plots automatic mode has
always made (cell_type_counts.png, cell_type_by_experiment_group.png) --
not enough to actually validate the GMM thresholds or the cell-type calls.
Modeled on the (already-working, already-validated) plotting code in the
LILRB2 mouse study's 05_1_a notebook, adapted for this pipeline's
column_map / short-marker-name convention and its already-z-scored data.

Entry point: generate_cell_typing_diagnostics(...), called once from
run_cell_typing() (spatia/analysis/cell_typing.py) right after
_save_and_plot() for both automatic and semi_automatic mode. Every
section below is independently wrapped in try/except -- one figure
failing (a plotting quirk, an unexpected column) never blocks the rest,
matching this pipeline's established convention (see preprocessing.py's
summary-report hook).

Deliberately NOT included in this first pass (ask if you want these too):
  - Per-marker UMAP (raw expression + binary positivity side by side) --
    the mouse notebook has this, but it's 55 more figure pairs on top of
    the 55 GMM distribution plots already generated here.
  - A faceted "one small UMAP panel per core" grid -- Afrouz chose one
    combined UMAP colored by core instead (2026-09-10 decision), which is
    what compute_and_plot_umaps() does. Revisit if that's not informative
    enough once you've looked at it.
"""

from __future__ import annotations

import os
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm

try:
    from spatia.analysis.preprocessing_summary import GROUP_COLORS, _color_for
except Exception:
    GROUP_COLORS = {}
    def _color_for(group, groups_in_order):
        cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        idx = groups_in_order.index(group) if group in groups_in_order else 0
        return cycle[idx % len(cycle)]


# ── Small shared helpers ────────────────────────────────────────────────────

def _savefig(fig, path_no_ext: str, formats=("png", "pdf")):
    """Save one figure in every requested format. 'png' always gets dpi=150."""
    for fmt in formats:
        kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            kwargs["dpi"] = 150
        fig.savefig(f"{path_no_ext}.{fmt}", **kwargs)


def _inverse_transform_scalar_local(value_t, transform: str, arcsinh_cofactor: float,
                                     transform_params: dict):
    """
    Local, minimal re-implementation of cell_typing.py's
    _inverse_transform_scalar -- duplicated (not imported) to keep this
    plotting module decoupled from cell_typing.py's private internals and
    avoid a circular import (cell_typing.py imports FROM this module).
    Covers "none"/"arcsinh"/"log1p" exactly. "yeojohnson" is intentionally
    NOT reproduced here (it needs the fitted lambda's full inverse formula)
    -- for that one transform, plot_marker_distributions() falls back to
    histogram + threshold line only, no GMM component curves, rather than
    risk a silently-wrong inverse. Not a gap for Afrouz's current config
    (gmm.transform: "none" -- confirmed from her actual run log).
    """
    if transform == "none":
        return float(value_t)
    elif transform == "arcsinh":
        return float(np.sinh(value_t) * arcsinh_cofactor)
    elif transform == "log1p":
        return float(np.expm1(value_t))
    else:
        raise ValueError(f"transform '{transform}' not supported for local inverse")


def _get_column(adata, marker: str, column_map: dict):
    real_col = column_map.get(marker, marker)
    if real_col not in adata.var_names:
        return None, real_col
    vals = adata[:, real_col].X
    if hasattr(vals, "toarray"):
        vals = vals.toarray().flatten()
    else:
        vals = np.asarray(vals).flatten()
    return vals, real_col


# ── 1. Per-marker GMM distribution plots ────────────────────────────────────

def plot_marker_distributions(adata, thresholds: dict, fit_info: dict, column_map: dict,
                               plot_dir: str, transform: str = "none",
                               arcsinh_cofactor: float = 5.0, formats=("png", "pdf")):
    """
    One figure per marker: histogram of raw expression, the fitted 2-component
    GMM (REUSING fit_info[marker]["gmm"] -- the exact fit that produced the
    real threshold, not a fresh re-fit), and the threshold line. This is more
    correct than re-fitting a naive GMM for the plot (what the reference
    mouse-study notebook does): it also respects this panel's documented
    zero-inflation handling (14/53 CRC markers have 30-50% exact-zero
    intensity -- see cell_typing.py's _fit_gmm docstring), since fit_info's
    "clean"/"gmm" already reflect that exclusion.

    Component curves are drawn as relative shapes on the raw x-axis (no
    Jacobian correction for the transform) -- this is a diagnostic overlay
    to sanity-check separation and threshold placement, not a rigorous
    raw-space density estimate.
    """
    column_map = column_map or {}
    fit_info = fit_info or {}
    sub_dir = os.path.join(plot_dir, "marker_distributions")
    os.makedirs(sub_dir, exist_ok=True)

    n_done, skipped = 0, []
    for marker, threshold in thresholds.items():
        info = fit_info.get(marker)
        vals, real_col = _get_column(adata, marker, column_map)
        if vals is None or info is None or info.get("gmm") is None:
            skipped.append(marker)
            continue

        clean = np.asarray(info.get("clean", vals)).flatten()
        clean = clean[np.isfinite(clean)]
        if len(clean) < 10:
            skipped.append(marker)
            continue

        gmm = info["gmm"]
        low_idx, high_idx = info["low_idx"], info["high_idx"]
        means_t = gmm.means_.flatten()
        stds_t = np.sqrt(gmm.covariances_.flatten())
        weights = gmm.weights_
        tparams = info.get("transform_params", {})

        try:
            lo, hi = np.percentile(clean, [0.1, 99.9])
            if lo >= hi:
                lo, hi = float(np.min(clean)), float(np.max(clean))
            x_raw = np.linspace(lo, hi, 500)

            fig, ax = plt.subplots(figsize=(9, 5.5))
            ax.hist(clean, bins=60, alpha=0.55, density=True, color="#8c8c8c",
                    edgecolor="white", linewidth=0.3, label=f"All cells (n={len(clean):,})")

            if transform == "yeojohnson":
                # No local inverse for this transform -- histogram + threshold
                # only (see _inverse_transform_scalar_local docstring).
                curves_drawn = False
            else:
                x_t = (np.arcsinh(x_raw / arcsinh_cofactor) if transform == "arcsinh"
                       else np.log1p(x_raw) if transform == "log1p" else x_raw)
                comp_names = {low_idx: "negative", high_idx: "positive"}
                mixture = np.zeros_like(x_t)
                for i in (low_idx, high_idx):
                    comp_pdf = weights[i] * norm.pdf(x_t, means_t[i], stds_t[i])
                    mixture += comp_pdf
                    ax.plot(x_raw, comp_pdf, linewidth=2,
                            label=f"{comp_names[i]} component (μ={means_t[i]:.2f}, w={weights[i]:.2f})")
                ax.plot(x_raw, mixture, "k--", linewidth=1.3, label="GMM mixture")
                curves_drawn = True

            ax.axvline(threshold, color="#C44E52", linestyle="--", linewidth=1.8,
                        label=f"Threshold: {threshold:.3f}")

            pos_col = f"{marker}_pos"
            if pos_col in adata.obs.columns:
                pct_pos = 100 * adata.obs[pos_col].mean()
                n_pos = int(adata.obs[pos_col].sum())
                title = f"{marker}  ({real_col})\nPositive: {n_pos:,} cells ({pct_pos:.1f}%)"
            else:
                title = f"{marker}  ({real_col})"
            ax.set_title(title, fontsize=11)
            ax.set_xlabel("Expression level (raw)")
            ax.set_ylabel("Density")
            ax.legend(fontsize="small", frameon=False)
            ax.spines[["top", "right"]].set_visible(False)
            plt.tight_layout()

            safe_marker = marker.replace("/", "_").replace(" ", "_")
            _savefig(fig, os.path.join(sub_dir, f"{safe_marker}_distribution"), formats)
            plt.close(fig)
            n_done += 1
        except Exception as e:
            plt.close("all")
            skipped.append(f"{marker} ({e})")

    print(f"  [diagnostics] Marker distribution plots: {n_done} written to {sub_dir}/"
          + (f"  -- skipped: {skipped}" if skipped else ""))


# ── 2. Marker-positivity-% heatmap by cell type ─────────────────────────────

def plot_positivity_heatmap(adata, markers: list, plot_dir: str, data_dir: str,
                             formats=("png", "pdf")):
    """
    Rows = cell type, columns = marker, value = % of that cell type's cells
    positive for that marker (from the existing {marker}_pos columns --
    same convention plot_marker_distributions() reads). Sequential single-hue
    colormap (0-100%), per the project's dataviz convention for magnitude.
    """
    if "cell_type" not in adata.obs.columns:
        print("  [diagnostics] Skipping positivity heatmap -- no 'cell_type' column.")
        return

    cell_types = adata.obs["cell_type"].value_counts().index.tolist()
    pos_cols = [m for m in markers if f"{m}_pos" in adata.obs.columns]
    if not pos_cols:
        print("  [diagnostics] Skipping positivity heatmap -- no *_pos columns found.")
        return

    table = pd.DataFrame(index=cell_types, columns=pos_cols, dtype=float)
    for ct in cell_types:
        sub = adata.obs.loc[adata.obs["cell_type"] == ct, [f"{m}_pos" for m in pos_cols]]
        table.loc[ct] = 100 * sub.mean().values

    table.to_csv(os.path.join(data_dir, "cell_type_marker_positivity_pct.csv"))

    fig_w = max(10, 0.35 * len(pos_cols) + 3)
    fig_h = max(6, 0.35 * len(cell_types) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(table, cmap="rocket_r", vmin=0, vmax=100, cbar_kws={"label": "% positive"}, ax=ax)
    ax.set_title("Marker Positivity (%) by Cell Type", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, "cell_type_marker_positivity_heatmap"), formats)
    plt.close(fig)
    print(f"  [diagnostics] Positivity heatmap: {len(cell_types)} cell types x {len(pos_cols)} markers")


# ── 3. Mean-expression heatmap by cell type ─────────────────────────────────

def plot_mean_expression_heatmap(adata, markers: list, column_map: dict, plot_dir: str,
                                  data_dir: str, formats=("png", "pdf")):
    """Rows = cell type, columns = marker, value = mean raw expression."""
    if "cell_type" not in adata.obs.columns:
        return
    column_map = column_map or {}
    cell_types = adata.obs["cell_type"].value_counts().index.tolist()
    usable = [(m, column_map.get(m, m)) for m in markers if column_map.get(m, m) in adata.var_names]
    if not usable:
        print("  [diagnostics] Skipping mean-expression heatmap -- no usable markers.")
        return

    table = pd.DataFrame(index=cell_types, columns=[m for m, _ in usable], dtype=float)
    for ct in cell_types:
        mask = (adata.obs["cell_type"] == ct).values
        sub = adata[mask]
        for m, real_col in usable:
            vals = sub[:, real_col].X
            vals = vals.toarray().flatten() if hasattr(vals, "toarray") else np.asarray(vals).flatten()
            table.loc[ct, m] = float(np.mean(vals)) if len(vals) else np.nan

    table.to_csv(os.path.join(data_dir, "cell_type_marker_means.csv"))

    fig_w = max(10, 0.35 * len(usable) + 3)
    fig_h = max(6, 0.35 * len(cell_types) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(table, cmap="mako_r", cbar_kws={"label": "Mean expression (raw)"}, ax=ax,
                center=0)
    ax.set_title("Mean Marker Expression by Cell Type", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, "cell_type_marker_expression_heatmap"), formats)
    plt.close(fig)
    print(f"  [diagnostics] Mean-expression heatmap: {len(cell_types)} cell types x {len(usable)} markers")


# ── 4. Restyled cell-type-counts bar chart ──────────────────────────────────

def plot_cell_type_counts_pretty(adata, plot_dir: str, formats=("png", "pdf")):
    """
    Replaces automatic mode's original cell_type_counts.png/.pdf (written
    moments earlier by _save_and_plot()) with a horizontal, sorted,
    value-labeled version -- the vertical/rotated-label version was flagged
    as unreadable with this many cell-type categories (long names collide
    at 45deg rotation). Horizontal bars sidestep that entirely.
    """
    if "cell_type" not in adata.obs.columns:
        return
    counts = adata.obs["cell_type"].value_counts().sort_values(ascending=True)
    total = int(counts.sum())
    pct = counts / total * 100

    fig_h = max(5, 0.32 * len(counts) + 1.5)
    fig, ax = plt.subplots(figsize=(9.5, fig_h))
    colors = ["#9a9a9a" if ct == "Unassigned" else "#4C72B0" for ct in counts.index]
    bars = ax.barh(counts.index.astype(str), counts.values, color=colors,
                    edgecolor="white", height=0.72, zorder=3)
    for bar, c, p in zip(bars, counts.values, pct.values):
        ax.text(bar.get_width() + total * 0.006, bar.get_y() + bar.get_height() / 2,
                f"{c:,}  ({p:.1f}%)", va="center", fontsize=8.5, color="#333333")

    ax.set_xlabel("Number of cells")
    ax.set_title(f"Cell Type Distribution  —  {total:,} cells total", fontsize=13, pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.grid(True, color="#e3e3e3", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.margins(x=0.14)
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, "cell_type_counts"), formats)
    plt.close(fig)
    print(f"  [diagnostics] Restyled cell_type_counts.png ({len(counts)} cell types)")


# ── 5. Stacked cell-type-by-group bar chart ─────────────────────────────────

def plot_stacked_by_group(adata, group_col: str, plot_dir: str, data_dir: str,
                           formats=("png", "pdf")):
    """
    100%-stacked bar chart: one bar per experiment_group (e.g. CLR, DII),
    each fully stacked into its cell-type composition -- segments sum to
    exactly 100% per bar. This is the standard reading of "stacked barplot
    of cell types" and is DIFFERENT from a first draft of this function
    that put cell_type on the x-axis and stacked CLR+DII on top of each
    other per type -- caught in verification (2026-09-10, synthetic test):
    that version's bars summed to (%-within-CLR + %-within-DII), which
    isn't bounded or meaningful (e.g. a bar could total ~96%, another
    ~28%), since each group's percentages are independently normalized to
    100% across cell types, not across groups. Group-on-x-axis, cell-type
    stacked is the only arrangement where the stack total is actually 100%.

    Adds this alongside (not replacing -- Afrouz confirmed the existing
    grouped/dodged version is fine as-is) _save_and_plot()'s existing
    cell_type_by_experiment_group.png, which shows the same per-group
    percentages as side-by-side bars instead of stacked ones.
    """
    if "cell_type" not in adata.obs.columns or group_col not in adata.obs.columns:
        return
    comp = pd.crosstab(adata.obs[group_col], adata.obs["cell_type"], normalize="index") * 100
    # Canonical, reproducible column (segment) order: most-common cell type
    # first, so colors/order are stable across re-runs of this cohort.
    order = adata.obs["cell_type"].value_counts().index.tolist()
    comp = comp[[c for c in order if c in comp.columns]]
    comp.to_csv(os.path.join(data_dir, f"cell_type_by_{group_col}_stacked100_pct.csv"))

    n_types = comp.shape[1]
    palette = sns.color_palette("husl", n_types)
    if "Unassigned" in comp.columns:
        palette[list(comp.columns).index("Unassigned")] = (0.6, 0.6, 0.6)

    fig, ax = plt.subplots(figsize=(max(6, 1.6 * comp.shape[0] + 2), 7))
    comp.plot(kind="bar", stacked=True, ax=ax, color=palette, edgecolor="white",
              linewidth=0.5, width=0.55)
    ax.set_title(f"Cell Type Composition by {group_col} (100% stacked)", fontsize=13)
    ax.set_ylabel("% of cells")
    ax.set_xlabel("")
    ax.set_ylim(0, 100)
    plt.setp(ax.get_xticklabels(), rotation=0)
    ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=7, ncol=1 if n_types <= 15 else 2)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, f"cell_type_by_{group_col}_stacked100"), formats)
    plt.close(fig)
    print(f"  [diagnostics] 100%-stacked cell-type-by-{group_col} chart written ({n_types} cell types)")


# ── 6. UMAP: cell type / group / core ───────────────────────────────────────

def compute_and_plot_umaps(adata, plot_dir: str, data_dir: str, group_col: str, tissue_col: str,
                            random_state: int = 42, formats=("png", "pdf")):
    """
    Computes PCA -> neighbors -> UMAP if not already present, then saves
    three colorings: cell_type (the "whole-TMA" map), group_col, tissue_col
    (one combined UMAP colored by core -- Afrouz's 2026-09-10 choice over a
    per-core facet grid). Also exports umap_coordinates.csv (UMAP1/UMAP2 +
    cell_type/group/tissue per cell) since the h5ad _save_and_plot() already
    wrote does NOT include these embeddings (they're computed after that
    write) -- re-writing the full 300k x 80 h5ad just to add 2 columns
    wasn't worth the extra I/O; the CSV is the cheap way to keep this
    reusable outside this script.

    IMPORTANT, and different from run_clustering() (semi_automatic mode):
    this does NOT call sc.pp.normalize_total / sc.pp.log1p / sc.pp.scale
    before PCA. This pipeline's preprocessing already per-cell z-scores
    every marker (confirmed: real thresholds on this exact dataset are
    legitimately negative, e.g. CD7=-0.0039, CD3=-0.0393, CD163=-0.0552 --
    log1p on any z-scored value below -1 produces NaN/-inf and would
    silently corrupt PCA). Going straight to PCA on the already-standardized
    matrix is correct here. NOTE: run_clustering() (used by semi_automatic
    mode) still has this exact log1p-on-already-z-scored-data issue -- not
    fixed here since automatic mode never calls it, but worth knowing
    before switching this config to semi_automatic mode.
    """
    import scanpy as sc

    if "X_umap" not in adata.obsm:
        n_comps = max(2, min(50, adata.n_vars - 1, adata.n_obs - 1))
        print(f"  [diagnostics] Computing PCA (n_comps={n_comps}) -> neighbors -> UMAP "
              f"on {adata.n_obs:,} cells (no re-normalization -- data is already z-scored)...")
        sc.pp.pca(adata, n_comps=n_comps, random_state=random_state)
        sc.pp.neighbors(adata, random_state=random_state)
        sc.tl.umap(adata, random_state=random_state)
    else:
        print("  [diagnostics] X_umap already present (semi_automatic clustering ran earlier) -- reusing it.")

    coords = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index)
    for col in ("cell_type", group_col, tissue_col):
        if col in adata.obs.columns:
            coords[col] = adata.obs[col].values
    coords.to_csv(os.path.join(data_dir, "umap_coordinates.csv"))

    colorings = [("cell_type", "UMAP — Cell Types (whole TMA)", "umap_cell_types_all")]
    if group_col in adata.obs.columns:
        colorings.append((group_col, f"UMAP — {group_col}", f"umap_by_{group_col}"))
    if tissue_col in adata.obs.columns:
        n_cores = adata.obs[tissue_col].nunique()
        colorings.append((tissue_col, f"UMAP — {tissue_col} ({n_cores} cores)", f"umap_by_{tissue_col}"))

    for color_key, title, fname in colorings:
        try:
            fig, ax = plt.subplots(figsize=(9, 8))
            legend_loc = "right margin" if adata.obs[color_key].nunique() <= 30 else None
            sc.pl.umap(adata, color=color_key, ax=ax, show=False, title=title,
                       legend_loc=legend_loc, legend_fontsize=7)
            _savefig(fig, os.path.join(plot_dir, fname), formats)
            plt.close(fig)
        except Exception as e:
            plt.close("all")
            print(f"  [diagnostics] UMAP coloring by '{color_key}' failed (non-fatal): {e}")

    print(f"  [diagnostics] UMAP: {len(colorings)} coloring(s) written, coordinates saved to "
          f"{data_dir}/umap_coordinates.csv")


# ── Entry point ──────────────────────────────────────────────────────────────

def generate_cell_typing_diagnostics(adata, thresholds: dict, fit_info: dict, column_map: dict,
                                      data_dir: str, plot_dir: str,
                                      transform: str = "none", arcsinh_cofactor: float = 5.0,
                                      group_col: str = "experiment_group",
                                      tissue_col: str = "tissue_id",
                                      random_state: int = 42, formats=("png", "pdf")):
    """
    Full diagnostic plot suite for a completed cell_typing run. Called once
    from run_cell_typing() after _save_and_plot(), for both automatic and
    semi_automatic mode. See module docstring for what's included and what
    was deliberately left out of this first pass.
    """
    print(f"\n{'=' * 80}\nGenerating cell-typing diagnostic plots...\n{'=' * 80}")

    try:
        plot_marker_distributions(adata, thresholds, fit_info, column_map, plot_dir,
                                   transform, arcsinh_cofactor, formats)
    except Exception as e:
        print(f"  WARNING: marker distribution plots failed (non-fatal): {e}")

    try:
        plot_positivity_heatmap(adata, list(thresholds.keys()), plot_dir, data_dir, formats)
    except Exception as e:
        print(f"  WARNING: positivity heatmap failed (non-fatal): {e}")

    try:
        plot_mean_expression_heatmap(adata, list(thresholds.keys()), column_map, plot_dir, data_dir, formats)
    except Exception as e:
        print(f"  WARNING: mean-expression heatmap failed (non-fatal): {e}")

    try:
        plot_cell_type_counts_pretty(adata, plot_dir, formats)
    except Exception as e:
        print(f"  WARNING: restyled cell_type_counts plot failed (non-fatal): {e}")

    try:
        plot_stacked_by_group(adata, group_col, plot_dir, data_dir, formats)
    except Exception as e:
        print(f"  WARNING: stacked cell-type-by-group plot failed (non-fatal): {e}")

    try:
        compute_and_plot_umaps(adata, plot_dir, data_dir, group_col, tissue_col, random_state, formats)
    except Exception as e:
        print(f"  WARNING: UMAP diagnostics failed (non-fatal): {e}")

    print(f"\n[diagnostics] Complete. Outputs in: {plot_dir} (figures) / {data_dir} (CSVs)")
