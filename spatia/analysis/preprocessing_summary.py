"""
spatia/analysis/preprocessing_summary.py
=========================================
Aggregate QC/summary report across a completed run_preprocessing() run.

This is the importable library: report-generation logic only, no CLI. It is
invoked automatically at the end of every run_preprocessing() call (see the
"Generating QC summary report" block near the end of that function in
spatia/analysis/preprocessing.py) -- so `--steps preprocessing` alone
regenerates summary_report/, with no separate command required. That call is
wrapped in try/except there: a bug in report generation is printed as a
warning and never fails the preprocessing step itself or blocks downstream
steps.

For a *manual* re-render without re-running (often hours-long) image
processing -- e.g. immediately after a report-code fix, or to regenerate a
report you deleted -- run the thin CLI wrapper at the repo root instead:

    python generate_preprocessing_summary.py --config experiments/crc_tma_full_pipeline.yaml

That CLI just loads the YAML config and calls generate_summary_report(cfg)
below -- the exact same function run_preprocessing() calls internally.

Reads the same config used for the pipeline run, finds the outputs
run_preprocessing() already wrote, and produces five framings:

  1. Tissue-level        : per-tissue cell counts and QC removal, from the
                            processing_stats_<timestamp>.csv row-per-image table.
  2. Group-level          : CLR vs DII (or whatever cfg["experiment"]["groups"]
                            is) comparison of the same QC metrics.
  3. Removal reasons      : coarse (size+DAPI filter vs. noise removal) cell
                            counts by group, from processing_stats.csv --
                            always available.
  4. Removal reasons (fine): small-area vs low-DAPI vs noise, individually --
                            only if qupath_exports/qupath_export_summary.csv
                            exists (written by the separate, optional
                            run_qupath_export(cfg) step -- not part of
                            run_pipeline.py --steps, so most runs won't have
                            it yet; this report is skipped with an
                            explanatory message, not an error, if absent).
  5. Marker-level         : per-marker mean expression by experiment_group,
                            read directly from the combined
                            *_combined_all_experiment_groups.h5ad files
                            (processing_stats.csv has no marker columns).
                            These values are already per-cell z-scores (see
                            run_preprocessing()'s zscore normalization step),
                            so markers are already on one comparable scale --
                            for exactly 2 groups this is plotted as a signed
                            group-mean-difference bar chart; for >2 groups, a
                            heatmap of the raw (not re-normalized) group means.

Does NOT re-run any pipeline step. Purely reads existing outputs.

Outputs (written under <output_dir>/combined_processed_data/summary_report/):
    tissue_level_qc.png / .csv
    group_comparison_qc.png / .csv
    removal_reasons_qc.png / .csv                (always)
    removal_reasons_fine_qc.png / .csv            (only if qupath_exports/ exists)
    marker_level_by_group.png / .csv
    SUMMARY.txt
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

try:
    import anndata as ad
except ImportError:
    ad = None


GROUP_COLORS = {
    # Falls back to matplotlib's default cycle for any group name not listed
    # here, so this isn't a hard dependency on exactly two groups named
    # CLR/DII -- it degrades gracefully for other configs.
}
_DEFAULT_CYCLE = plt.rcParams["axes.prop_cycle"].by_key()["color"]


def _color_for(group: str, groups_in_order: list) -> str:
    if group in GROUP_COLORS:
        return GROUP_COLORS[group]
    idx = groups_in_order.index(group) if group in groups_in_order else 0
    return _DEFAULT_CYCLE[idx % len(_DEFAULT_CYCLE)]


def _latest_stats_csv(log_dir: str) -> str:
    candidates = sorted(glob.glob(os.path.join(log_dir, "processing_stats_*.csv")))
    if not candidates:
        raise FileNotFoundError(
            f"No processing_stats_*.csv found in {log_dir} -- "
            f"has run_preprocessing() completed at least once?"
        )
    return candidates[-1]  # lexically-sortable timestamp -> last is newest


def _load_stats(log_dir: str) -> pd.DataFrame:
    path = _latest_stats_csv(log_dir)
    print(f"Reading tissue/image-level stats: {path}")
    df = pd.read_csv(path)
    return df


# ── 1. Tissue-level ─────────────────────────────────────────────────────

def make_tissue_level_report(stats: pd.DataFrame, groups_in_order: list, out_dir: str):
    done = stats[stats.get("status", "") == "PROCESSED"].copy()
    if done.empty:
        print("  ⚠️  No rows with status == PROCESSED -- skipping tissue-level report")
        return

    done = done.sort_values("total_percent_removed", ascending=False)
    csv_path = os.path.join(out_dir, "tissue_level_qc.csv")
    done[[
        "tissue_id", "image_id", "experiment_group", "slide_folder",
        "original_cells", "final_cells",
        "percent_removed_by_filter", "percent_removed_by_noise", "total_percent_removed",
    ]].to_csv(csv_path, index=False)
    print(f"  ✓ Saved {csv_path}")

    fig, axes = plt.subplots(2, 1, figsize=(max(12, len(done) * 0.18), 10), sharex=True)
    x = np.arange(len(done))
    colors = [_color_for(g, groups_in_order) for g in done["experiment_group"]]

    axes[0].bar(x, done["original_cells"], color="#cfcfcf", label="Original cells")
    axes[0].bar(x, done["final_cells"], color=colors, label="Final cells (post-QC)")
    axes[0].set_ylabel("Cell count")
    axes[0].set_title("Per-image cell counts, original vs. after QC (sorted by % removed)")
    axes[0].legend(loc="upper right")

    axes[1].bar(x, done["total_percent_removed"], color=colors)
    axes[1].set_ylabel("Total % cells removed")
    axes[1].set_xlabel("Image (sorted, most cells removed → least)")
    axes[1].set_title("QC removal rate per image")
    axes[1].axhline(done["total_percent_removed"].mean(), color="black", linestyle="--",
                     linewidth=1, label=f"Mean: {done['total_percent_removed'].mean():.1f}%")
    axes[1].legend(loc="upper right")

    if len(done) <= 60:
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(done["image_id"], rotation=90, fontsize=6)
    else:
        axes[1].set_xticks([])

    handles = [plt.Rectangle((0, 0), 1, 1, color=_color_for(g, groups_in_order)) for g in groups_in_order]
    fig.legend(handles, groups_in_order, loc="upper center", ncol=len(groups_in_order),
               bbox_to_anchor=(0.5, 1.02), title="experiment_group")

    plt.tight_layout()
    png_path = os.path.join(out_dir, "tissue_level_qc.png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved {png_path}")


# ── 2. Group-level (CLR vs DII) ──────────────────────────────────────────

def make_group_comparison_report(stats: pd.DataFrame, groups_in_order: list, out_dir: str):
    done = stats[stats.get("status", "") == "PROCESSED"].copy()
    if done.empty:
        print("  ⚠️  No rows with status == PROCESSED -- skipping group comparison report")
        return

    metrics = ["original_cells", "final_cells", "percent_removed_by_filter",
               "percent_removed_by_noise", "total_percent_removed"]
    summary = done.groupby("experiment_group")[metrics].agg(["count", "mean", "median", "std"])
    csv_path = os.path.join(out_dir, "group_comparison_qc.csv")
    summary.to_csv(csv_path)
    print(f"  ✓ Saved {csv_path}")

    plot_metrics = ["final_cells", "total_percent_removed"]
    fig, axes = plt.subplots(1, len(plot_metrics), figsize=(6 * len(plot_metrics), 5))
    if len(plot_metrics) == 1:
        axes = [axes]

    for ax, metric in zip(axes, plot_metrics):
        data_by_group = [done.loc[done["experiment_group"] == g, metric].dropna().values
                          for g in groups_in_order if g in done["experiment_group"].unique()]
        labels = [g for g in groups_in_order if g in done["experiment_group"].unique()]
        colors = [_color_for(g, groups_in_order) for g in labels]

        bp = ax.boxplot(data_by_group, labels=labels, patch_artist=True, showmeans=True)
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
        for i, d in enumerate(data_by_group, start=1):
            jitter = np.random.default_rng(0).normal(0, 0.04, size=len(d))
            ax.scatter(np.full(len(d), i) + jitter, d, s=8, color="black", alpha=0.4, zorder=3)

        ax.set_title(metric.replace("_", " "))
        ax.set_ylabel(metric)

    n_tissues = done.groupby("experiment_group")["tissue_id"].nunique()
    n_images = done.groupby("experiment_group").size()
    subtitle = "  |  ".join(
        f"{g}: {n_images.get(g, 0)} images, {n_tissues.get(g, 0)} tissues" for g in groups_in_order
    )
    fig.suptitle(f"CLR/DII comparison — {subtitle}", fontsize=11)
    plt.tight_layout()
    png_path = os.path.join(out_dir, "group_comparison_qc.png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved {png_path}")

    return summary


# ── 3. Cell-removal reasons ────────────────────────────────────────────────

def make_removal_reasons_report(stats: pd.DataFrame, groups_in_order: list, out_dir: str):
    """
    Coarse removal-reason breakdown, from processing_stats.csv alone -- always
    available after any run_preprocessing() run, no extra step needed. Two
    reasons are tracked at this level: the combined size+DAPI QC filter, and
    z-score noise removal. For the finer per-reason split (small-area vs
    low-DAPI vs noise, individually -- see preprocessing.py's
    run_qupath_export()/_classify_cells()), see
    make_fine_removal_reasons_report() below, which needs that separate,
    optional step to have been run first.
    """
    done = stats[stats.get("status", "") == "PROCESSED"].copy()
    if done.empty:
        print("  ⚠️  No rows with status == PROCESSED -- skipping removal-reasons report")
        return

    reason_cols = ["final_cells", "cells_removed_by_filter", "cells_removed_by_noise"]
    agg = done.groupby("experiment_group")[reason_cols].sum()
    agg = agg.reindex([g for g in groups_in_order if g in agg.index])
    agg["total_raw_cells"] = agg[reason_cols].sum(axis=1)
    for col in reason_cols:
        agg[f"{col}_pct"] = agg[col] / agg["total_raw_cells"] * 100

    csv_path = os.path.join(out_dir, "removal_reasons_qc.csv")
    agg.to_csv(csv_path)
    print(f"  ✓ Saved {csv_path}")

    stack_cols = ["final_cells", "cells_removed_by_filter", "cells_removed_by_noise"]
    labels     = ["Kept", "Removed: size+DAPI filter", "Removed: noise removal"]
    colors     = ["#4C9F70", "#E8A33D", "#C1666B"]

    fig, (ax_n, ax_pct) = plt.subplots(1, 2, figsize=(11, 5))
    x = np.arange(len(agg.index))
    bottom_n, bottom_pct = np.zeros(len(agg)), np.zeros(len(agg))
    for col, label, color in zip(stack_cols, labels, colors):
        vals_n, vals_pct = agg[col].values, agg[f"{col}_pct"].values
        ax_n.bar(x, vals_n, bottom=bottom_n, label=label, color=color)
        ax_pct.bar(x, vals_pct, bottom=bottom_pct, label=label, color=color)
        bottom_n += vals_n
        bottom_pct += vals_pct

    for ax, title, ylabel in [(ax_n, "Cell counts", "Cells"), (ax_pct, "As % of raw cells", "% of raw cells")]:
        ax.set_xticks(x)
        ax.set_xticklabels(agg.index)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
    ax_pct.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=1, fontsize=8)
    fig.suptitle("Cell removal reasons by experiment_group (coarse: filter vs. noise)", fontsize=12)
    plt.tight_layout()
    png_path = os.path.join(out_dir, "removal_reasons_qc.png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved {png_path}")

    return agg


def make_fine_removal_reasons_report(output_dir: str, groups_in_order: list, out_dir: str):
    """
    Fine-grained removal-reason breakdown (Excl_SmallArea / Excl_LowDAPI /
    Excl_SmallArea_LowDAPI / Excl_Noise / Included), read from
    qupath_exports/qupath_export_summary.csv -- written by
    run_qupath_export(cfg), a SEPARATE, optional function that is NOT part of
    run_pipeline.py's --steps registry and is not called automatically by
    run_preprocessing(). Most runs won't have this file yet. Skips (with an
    explanatory message, not an error) if it's missing.
    """
    summary_path = os.path.join(output_dir, "qupath_exports", "qupath_export_summary.csv")
    if not os.path.exists(summary_path):
        print(f"  ⚠️  {summary_path} not found -- skipping fine-grained removal reasons. "
              f"This needs run_qupath_export(cfg) to be run separately first (it's not "
              f"wired into run_pipeline.py --steps): "
              f"python -c \"import yaml; from spatia.analysis.preprocessing import run_qupath_export; "
              f"run_qupath_export(yaml.safe_load(open('<your config>')))\"")
        return

    df = pd.read_csv(summary_path)
    reason_cols = [c for c in ["Included", "Excl_SmallArea", "Excl_LowDAPI",
                                "Excl_SmallArea_LowDAPI", "Excl_Noise"] if c in df.columns]
    if not reason_cols or "experiment_group" not in df.columns:
        print(f"  ⚠️  {summary_path} doesn't have the expected columns -- skipping")
        return

    agg = df.groupby("experiment_group")[reason_cols].sum(min_count=1).fillna(0)
    agg = agg.reindex([g for g in groups_in_order if g in agg.index])
    total = agg.sum(axis=1)
    pct = agg.div(total, axis=0) * 100

    csv_path = os.path.join(out_dir, "removal_reasons_fine_qc.csv")
    agg.to_csv(csv_path)
    print(f"  ✓ Saved {csv_path}")

    # Same palette run_qupath_export()'s own QC plots use, for visual consistency.
    colors_map = {
        "Included":               "#00D200",
        "Excl_SmallArea":         "#FF3C3C",
        "Excl_LowDAPI":           "#FFA500",
        "Excl_SmallArea_LowDAPI": "#B400B4",
        "Excl_Noise":             "#1E90FF",
    }

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(pct.index))
    bottom = np.zeros(len(pct))
    for col in reason_cols:
        vals = pct[col].values
        ax.bar(x, vals, bottom=bottom, label=col, color=colors_map.get(col, "grey"))
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(pct.index)
    ax.set_ylabel("% of raw cells")
    ax.set_title("Cell classification by experiment_group (fine-grained, from run_qupath_export)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=8)
    plt.tight_layout()
    png_path = os.path.join(out_dir, "removal_reasons_fine_qc.png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved {png_path}")

    return agg


# ── 3. Marker-level ───────────────────────────────────────────────────────

def make_marker_level_report(tissues_dir: str, groups_in_order: list, out_dir: str):
    if ad is None:
        print("  ⚠️  anndata not importable in this environment -- skipping marker-level "
              "report. Run this script inside the spacec venv on JupyterHub "
              "(the same one run_pipeline.py runs in), not a plain python3.")
        return

    h5ad_files = sorted(glob.glob(os.path.join(tissues_dir, "*_combined_all_experiment_groups.h5ad")))
    if not h5ad_files:
        print(f"  ⚠️  No *_combined_all_experiment_groups.h5ad found in {tissues_dir} -- "
              f"skipping marker-level report")
        return

    print(f"  Reading {len(h5ad_files)} combined h5ad file(s) for marker-level means "
          f"(this loads every tissue's expression matrix -- may take a few minutes)...")

    per_tissue_group_means = []  # list of (marker -> mean) rows, tagged by group
    row_meta = []
    for i, path in enumerate(h5ad_files, start=1):
        try:
            adata = ad.read_h5ad(path)
        except Exception as e:
            print(f"    ⚠️  Could not open {os.path.basename(path)}: {e}")
            continue
        if "experiment_group" not in adata.obs.columns or adata.n_obs == 0:
            continue
        df = adata.to_df()  # cells x markers, obs index aligned
        df["experiment_group"] = adata.obs["experiment_group"].values
        group_means = df.groupby("experiment_group").mean(numeric_only=True)
        for g in group_means.index:
            per_tissue_group_means.append(group_means.loc[g])
            row_meta.append({"tissue_file": os.path.basename(path), "experiment_group": g,
                              "n_cells": int((adata.obs["experiment_group"] == g).sum())})
        if i % 20 == 0 or i == len(h5ad_files):
            print(f"    ...{i}/{len(h5ad_files)} files read")

    if not per_tissue_group_means:
        print("  ⚠️  No usable per-group marker means extracted -- skipping marker-level report")
        return

    means_df = pd.DataFrame(per_tissue_group_means)
    meta_df = pd.DataFrame(row_meta)
    means_df["experiment_group"] = meta_df["experiment_group"].values

    # Weighted-by-cell-count mean of means across tissues, per group per marker.
    marker_cols = [c for c in means_df.columns if c != "experiment_group"]
    weights = meta_df["n_cells"].values
    rows = []
    for g in groups_in_order:
        mask = (means_df["experiment_group"] == g).values
        if not mask.any():
            continue
        w = weights[mask]
        sub = means_df.loc[mask, marker_cols]
        weighted_mean = np.average(sub.values, axis=0, weights=w)
        rows.append(pd.Series(weighted_mean, index=marker_cols, name=g))
    if not rows:
        print("  ⚠️  No groups had usable data -- skipping marker-level report")
        return
    group_marker_means = pd.DataFrame(rows)  # groups x markers

    csv_path = os.path.join(out_dir, "marker_level_by_group.csv")
    group_marker_means.T.to_csv(csv_path)
    print(f"  ✓ Saved {csv_path}")

    # IMPORTANT: these are already per-cell z-scores, not raw intensities --
    # run_preprocessing() z-score normalizes every marker (sp.pp.format(...,
    # method="zscore")) before this h5ad is ever written. Markers are already
    # on one common, comparable scale, so no further re-normalization is
    # correct here.
    #
    # Fixed 2026-09-10: the first version of this function re-z-scored these
    # already-z-scored group means AGAIN, per marker row, across groups. For
    # the common case of exactly 2 groups that's mathematically degenerate --
    # z-scoring any 2 points always produces the same magnitude (+-1/sqrt(2)),
    # regardless of how different the two groups actually are. Only the SIGN
    # carried real information; every tile in that heatmap had the same
    # visual intensity, which is exactly what made it uninterpretable. Fixed
    # by dropping the re-z-scoring: for 2 groups, plot the actual group-mean
    # difference (which IS a meaningful, magnitude-preserving quantity once
    # the inputs are already z-scores); for >2 groups, a heatmap of the raw
    # (not re-normalized) group means.
    present_groups = [g for g in groups_in_order if g in group_marker_means.index]

    if len(present_groups) == 2:
        g1, g2 = present_groups
        diff = (group_marker_means.loc[g1] - group_marker_means.loc[g2]).sort_values()
        bar_colors = ["#C1666B" if v < 0 else "#4C9F70" for v in diff.values]
        fig, ax = plt.subplots(figsize=(7, max(8, len(diff) * 0.22)))
        ax.barh(diff.index, diff.values, color=bar_colors)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel(f"Δ mean per-cell z-score ({g1} − {g2})")
        ax.set_title(f"Marker expression difference by experiment_group\n"
                     f"({g1} higher →       ← {g2} higher)")
        ax.tick_params(axis="y", labelsize=7)
        plt.tight_layout()
    else:
        heat = group_marker_means.T.reindex(columns=present_groups)
        vmax = float(np.nanmax(np.abs(heat.values))) if heat.size else 1.0
        vmax = vmax if vmax > 0 else 1.0
        fig, ax = plt.subplots(figsize=(max(6, len(present_groups) * 1.5), max(8, len(heat) * 0.28)))
        im = ax.imshow(heat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(len(heat.columns)))
        ax.set_xticklabels(heat.columns)
        ax.set_yticks(range(len(heat.index)))
        ax.set_yticklabels(heat.index, fontsize=7)
        ax.set_title("Mean per-cell z-score by experiment_group\n(already on a common scale -- not re-normalized)")
        plt.colorbar(im, ax=ax, label="mean z-score", fraction=0.05, pad=0.02)
        plt.tight_layout()

    png_path = os.path.join(out_dir, "marker_level_by_group.png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved {png_path}")

    return group_marker_means


# ── Main ──────────────────────────────────────────────────────────────────


def generate_summary_report(cfg: dict) -> str:
    """
    Generate the full 5-report QC summary for a completed run_preprocessing()
    output, driven by the same cfg dict run_preprocessing() itself received.

    Called automatically from the end of run_preprocessing() -- every
    `--steps preprocessing` run regenerates summary_report/ with no separate
    command needed (see spatia/analysis/preprocessing.py). Can also be
    called standalone via generate_preprocessing_summary.py at the repo
    root, to cheaply re-render the report -- e.g. right after a report-code
    fix -- without re-running the (often hours-long) image processing.

    Returns the summary_report/ output directory path.
    """
    groups_in_order = list(cfg["experiment"]["groups"])
    base_out = cfg["paths"]["output_dir"]
    combined_dir = os.path.join(base_out, "combined_processed_data")
    log_dir = os.path.join(combined_dir, "processing_logs")
    tissues_dir = os.path.join(combined_dir, "individual_processed_data")
    out_dir = os.path.join(combined_dir, "summary_report")
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 80)
    print("PREPROCESSING SUMMARY REPORT")
    print("=" * 80)
    print(f"Groups:        {groups_in_order}")
    print(f"Stats source:  {log_dir}")
    print(f"h5ad source:   {tissues_dir}")
    print(f"Output:        {out_dir}")
    print("=" * 80)

    stats = _load_stats(log_dir)

    print("\n[1/5] Tissue-level report...")
    make_tissue_level_report(stats, groups_in_order, out_dir)

    print("\n[2/5] Group comparison report...")
    group_summary = make_group_comparison_report(stats, groups_in_order, out_dir)

    print("\n[3/5] Removal-reasons report (coarse: filter vs. noise)...")
    removal_summary = make_removal_reasons_report(stats, groups_in_order, out_dir)

    print("\n[4/5] Removal-reasons report (fine-grained, if available)...")
    fine_removal_summary = make_fine_removal_reasons_report(base_out, groups_in_order, out_dir)

    print("\n[5/5] Marker-level report (reads all combined h5ad files)...")
    marker_summary = make_marker_level_report(tissues_dir, groups_in_order, out_dir)

    summary_txt_path = os.path.join(out_dir, "SUMMARY.txt")
    with open(summary_txt_path, "w") as f:
        f.write("PREPROCESSING SUMMARY REPORT\n")
        f.write("=" * 80 + "\n")
        done = stats[stats.get("status", "") == "PROCESSED"]
        f.write(f"Images processed: {len(done)} / {len(stats)}\n")
        f.write(f"Unique tissues:   {done['tissue_id'].nunique() if len(done) else 0}\n")
        if group_summary is not None:
            f.write("\nBy experiment_group (QC metrics):\n")
            f.write(group_summary.to_string())
            f.write("\n")
        if removal_summary is not None:
            f.write("\nRemoval reasons by experiment_group (coarse):\n")
            f.write(removal_summary.to_string())
            f.write("\n")
        if fine_removal_summary is not None:
            f.write("\nRemoval reasons by experiment_group (fine-grained):\n")
            f.write(fine_removal_summary.to_string())
            f.write("\n")
        else:
            f.write("\nFine-grained removal reasons: not available -- run_qupath_export(cfg) "
                    "hasn't been run for this output_dir. See removal_reasons_qc.png/.csv for "
                    "the coarse (filter vs. noise) breakdown instead.\n")
        if marker_summary is not None:
            f.write(f"\nMarkers summarized: {marker_summary.shape[1]}\n")
    print(f"\n✓ Saved {summary_txt_path}")

    print(f"\nDone. All outputs in: {out_dir}")
    return out_dir
