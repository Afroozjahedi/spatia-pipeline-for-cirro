#!/usr/bin/env python3
"""
generate_preprocessing_summary.py
==================================
Aggregate QC/summary report across a completed run_preprocessing() run.

Reads the same config used for the pipeline run (--config), finds the
outputs run_preprocessing() already wrote, and produces three framings:

  1. Tissue-level : per-tissue cell counts and QC removal, from the
                     processing_stats_<timestamp>.csv row-per-image table.
  2. Group-level   : CLR vs DII (or whatever cfg["experiment"]["groups"] is)
                     comparison of the same QC metrics.
  3. Marker-level  : per-marker mean expression by experiment_group, read
                     directly from the combined *_combined_all_experiment_groups.h5ad
                     files (processing_stats.csv has no marker columns --
                     only the h5ad/csv outputs carry per-marker intensities).

Does NOT re-run any pipeline step. Purely reads existing outputs.

Usage
-----
    python generate_preprocessing_summary.py --config experiments/crc_tma_full_pipeline.yaml

Outputs (written under <output_dir>/combined_processed_data/summary_report/):
    tissue_level_qc.png
    tissue_level_qc.csv
    group_comparison_qc.png
    group_comparison_qc.csv
    marker_level_by_group.png
    marker_level_by_group.csv
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

    # Heatmap: markers (rows) x groups (cols), z-scored across groups per marker
    # so markers with very different absolute scales are still visually comparable.
    heat = group_marker_means.T.copy()
    heat_z = heat.sub(heat.mean(axis=1), axis=0).div(heat.std(axis=1).replace(0, np.nan), axis=0)
    heat_z = heat_z.dropna(how="all")

    fig, ax = plt.subplots(figsize=(max(6, len(groups_in_order) * 1.5), max(8, len(heat_z) * 0.28)))
    im = ax.imshow(heat_z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(heat_z.columns)))
    ax.set_xticklabels(heat_z.columns)
    ax.set_yticks(range(len(heat_z.index)))
    ax.set_yticklabels(heat_z.index, fontsize=7)
    ax.set_title("Mean marker expression by experiment_group\n(row z-scored across groups)")
    plt.colorbar(im, ax=ax, label="z-score", fraction=0.05, pad=0.02)
    plt.tight_layout()
    png_path = os.path.join(out_dir, "marker_level_by_group.png")
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ Saved {png_path}")

    return group_marker_means


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, help="Same YAML config used for the pipeline run")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Config not found: {args.config}")
        sys.exit(2)

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

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
    print(f"Config:        {args.config}")
    print(f"Groups:        {groups_in_order}")
    print(f"Stats source:  {log_dir}")
    print(f"h5ad source:   {tissues_dir}")
    print(f"Output:        {out_dir}")
    print("=" * 80)

    stats = _load_stats(log_dir)

    print("\n[1/3] Tissue-level report...")
    make_tissue_level_report(stats, groups_in_order, out_dir)

    print("\n[2/3] Group comparison report (CLR vs DII)...")
    group_summary = make_group_comparison_report(stats, groups_in_order, out_dir)

    print("\n[3/3] Marker-level report (reads all combined h5ad files)...")
    marker_summary = make_marker_level_report(tissues_dir, groups_in_order, out_dir)

    # ── Plain-text top-line summary ──────────────────────────────────────
    summary_txt_path = os.path.join(out_dir, "SUMMARY.txt")
    with open(summary_txt_path, "w") as f:
        f.write("PREPROCESSING SUMMARY REPORT\n")
        f.write("=" * 80 + "\n")
        done = stats[stats.get("status", "") == "PROCESSED"]
        f.write(f"Images processed: {len(done)} / {len(stats)}\n")
        f.write(f"Unique tissues:   {done['tissue_id'].nunique() if len(done) else 0}\n")
        if group_summary is not None:
            f.write("\nBy experiment_group:\n")
            f.write(group_summary.to_string())
            f.write("\n")
        if marker_summary is not None:
            f.write(f"\nMarkers summarized: {marker_summary.shape[1]}\n")
    print(f"\n✓ Saved {summary_txt_path}")

    print(f"\nDone. All outputs in: {out_dir}")


if __name__ == "__main__":
    main()
