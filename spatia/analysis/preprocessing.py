"""
spatia/analysis/preprocessing.py
=================================
Config-driven preprocessing for spatial proteomics data.

Refactored from: 04-0_preprocessing.ipynb

Entry points
------------
run_preprocessing(cfg)         -- main pipeline: filter, normalise, noise-remove, save
run_qupath_export(cfg)         -- optional QC export for QuPath visualisation

Both accept the parsed YAML config dict produced by:
    import yaml
    with open("config_example.yaml") as f:
        cfg = yaml.safe_load(f)
"""

from __future__ import annotations

import os
import re
import sys
import pickle
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import spacec as sp
from kneed import KneeLocator


# ─────────────────────────────────────────────────────────────────────────────
# LOGGING
# ─────────────────────────────────────────────────────────────────────────────

class _DualLogger:
    """Tee stdout → console + log file."""
    def __init__(self, filename: str):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding="utf-8")

    def write(self, message: str):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()

    def close(self):
        self.log.close()


# ─────────────────────────────────────────────────────────────────────────────
# TISSUE IDENTIFIER UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def extract_tissue_identifier(image_id: str, conditions: List[str]) -> str:
    """
    Strip condition labels and coordinate blocks from an image_id.

    Example
    -------
    'WT_PirB_D14-2_3_Scan1_x16356_y1375_masked'
    → 'PirB_D14-2_3_Scan1_masked'
    """
    tissue_id = image_id.strip().rstrip("_").strip()

    for cond in conditions:
        patterns = [
            (rf"^{cond}[_\s]+", ""),        # prefix:  KO_…
            (rf"[_\s]+{cond}[_\s]*$", ""),  # suffix:  …_KO
            (rf"[_\s]+{cond}[_\s]+", "_"),  # middle:  …_KO_…
        ]
        for pattern, replacement in patterns:
            tissue_id = re.sub(pattern, replacement, tissue_id, flags=re.IGNORECASE)

    # Remove coordinate blocks: _x1234_y5678  or  _x1234_y5678_w100_h200
    tissue_id = re.sub(r"_x\d+_y\d+(?:_w\d+_h\d+)?", "", tissue_id)

    # Collapse multiple underscores / spaces
    tissue_id = re.sub(r"[_\s]+", "_", tissue_id).strip("_")
    return tissue_id


def detect_condition(image_id: str, conditions: List[str]) -> str:
    """
    Return the condition label found in *image_id*, or 'Unknown'.

    Detection order (most → least reliable):
      1. Prefix match  (KO_…  or  KO …)
      2. Suffix match  (…_KO  or  … KO)
      3. Unambiguous substring match
    """
    uid = image_id.upper()
    for cond in conditions:
        cu = cond.upper()
        if uid.startswith(f"{cu}_") or uid.startswith(f"{cu} "):
            return cond

    for cond in conditions:
        cu = cond.upper()
        if uid.endswith(f"_{cu}") or uid.endswith(f" {cu}"):
            return cond

    matches = [c for c in conditions if c.upper() in uid]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        print(f"  ⚠️  Ambiguous condition for '{image_id}' (found: {matches}) → 'Unknown'")
    else:
        print(f"  ⚠️  No condition detected for '{image_id}' → 'Unknown'")
    return "Unknown"


# ─────────────────────────────────────────────────────────────────────────────
# QC HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _get_last_marker_col(df: pd.DataFrame, last_marker: str) -> int:
    """
    Return the integer column index of *last_marker*.
    Falls back to the rightmost non-metadata column if the marker is absent.

    Why this matters
    ----------------
    spacec's remove_noise() and make_anndata() use a column index (col_sum /
    col_num) to split the DataFrame into marker columns vs metadata columns.
    Hardcoding 'SIGLEC F' breaks on any panel that doesn't include that marker.
    """
    METADATA_COLS = {
        "label", "area", "x", "y", "slide_folder", "image_ID", "condition",
        "eccentricity", "perimeter", "convex_area",
        "axis_major_length", "axis_minor_length", "DAPI",
    }
    if last_marker and last_marker in df.columns:
        return df.columns.get_loc(last_marker)

    marker_cols = [c for c in df.columns if c not in METADATA_COLS]
    if not marker_cols:
        raise ValueError("No marker columns found in DataFrame.")
    return df.columns.get_loc(marker_cols[-1])


def auto_detect_cutoffs(
    df: pd.DataFrame,
    col_num: int,
    cut_off: float = 0.01,
    count_bin: int = 50,
) -> Tuple[float, float]:
    """
    Knee-based automatic detection of z_count and z_sum noise cutoffs.

    Parameters
    ----------
    df        : normalised marker DataFrame
    col_num   : index of the last marker column
    cut_off   : percentile fallback fraction (default 0.01 → 99th percentile)
    count_bin : number of histogram bins

    Returns
    -------
    (z_count_cutoff, z_sum_cutoff)
    """
    marker_cols = df.iloc[:, : col_num + 1].columns
    z = np.array([df[c].values for c in marker_cols]).T

    z_count = np.sum(z > 0, axis=1)
    z_sum = np.sum(z, axis=1)

    def _knee_or_percentile(values, n_bins, fallback_pct):
        hist, bins = np.histogram(values, bins=n_bins)
        centers = 0.5 * (bins[1:] + bins[:-1])
        try:
            kl = KneeLocator(centers, hist, S=1.0, curve="convex", direction="decreasing")
            knee = kl.knee
        except Exception:
            knee = None
        return knee if knee else np.percentile(values, 100 * (1 - fallback_pct))

    z_count_cutoff = _knee_or_percentile(z_count, count_bin, cut_off)
    z_sum_cutoff   = _knee_or_percentile(z_sum,   count_bin, cut_off)

    # Floor sanity checks
    if z_count_cutoff < 10:
        z_count_cutoff = np.percentile(z_count, 99)
    if z_sum_cutoff < 20:
        z_sum_cutoff = np.percentile(z_sum, 99)

    return float(z_count_cutoff), float(z_sum_cutoff)


def _processed_files_exist(tissue_id: str, condition: Optional[str], out_dir: str) -> bool:
    """Return True if the per-condition (or combined) CSV + h5ad both exist."""
    if condition:
        base = os.path.join(out_dir, f"{tissue_id}_{condition}")
    else:
        base = os.path.join(out_dir, f"{tissue_id}_combined_all_conditions")
    return os.path.exists(base + ".csv") and os.path.exists(base + ".h5ad")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PREPROCESSING PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def run_preprocessing(cfg: dict) -> dict:
    """
    Full preprocessing pipeline driven by the SPATIA config dict.

    Steps
    -----
    1. Scan segmentation_results_dir for mesmer_result.csv files
    2. For each image: load → filter (size + DAPI) → z-score normalise →
       auto-detect noise cutoffs → remove_noise → save per-condition files
    3. Combine per-condition CSVs into a tissue-level combined h5ad
    4. Generate marker-expression overlay visualisations

    Parameters
    ----------
    cfg : parsed YAML config (yaml.safe_load)

    Returns
    -------
    dict with keys:
        all_processed_tissues : {tissue_id: [tissue_info_dict, ...]}
        processing_stats      : list of per-image dicts
        output_dir            : path to combined_processed_data/
    """

    # ── Config extraction ─────────────────────────────────────────────────
    conditions     = cfg["experiment"]["conditions"]
    seg_dir        = cfg["paths"]["segmentation_results_dir"]
    base_out       = cfg["paths"]["output_dir"]

    pp_cfg         = cfg.get("preprocessing", {})
    last_marker    = pp_cfg.get("last_marker", None)   # e.g. "SIGLEC F"
    size_pct       = pp_cfg.get("qc_filter", {}).get("size_percentile",  1)
    dapi_pct       = pp_cfg.get("qc_filter", {}).get("dapi_percentile",  1)
    noise_cut_off  = pp_cfg.get("noise", {}).get("cut_off",   0.01)
    noise_bins     = pp_cfg.get("noise", {}).get("count_bin",  50)

    # ── Output directories ────────────────────────────────────────────────
    out_dir       = os.path.join(base_out, "combined_processed_data")
    tissues_dir   = os.path.join(out_dir,  "individual_processed_data")
    log_dir       = os.path.join(out_dir,  "processing_logs")
    viz_dir       = os.path.join(out_dir,  "marker_visualizations")

    for d in [out_dir, tissues_dir, log_dir, viz_dir]:
        os.makedirs(d, exist_ok=True)

    # ── Logging ───────────────────────────────────────────────────────────
    timestamp  = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file   = os.path.join(log_dir, f"processing_log_{timestamp}.txt")
    dual_log   = _DualLogger(log_file)
    sys.stdout = dual_log

    print("=" * 80)
    print("SPATIA PREPROCESSING PIPELINE")
    print("=" * 80)
    print(f"Started:    {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Conditions: {conditions}")
    print(f"Seg dir:    {seg_dir}")
    print(f"Output dir: {out_dir}")
    print(f"Log file:   {log_file}")
    print("=" * 80)

    # ── Discover slide folders ────────────────────────────────────────────
    if not os.path.isdir(seg_dir):
        raise FileNotFoundError(f"segmentation_results_dir not found: {seg_dir}")

    slide_folders = [
        f for f in os.listdir(seg_dir)
        if os.path.isdir(os.path.join(seg_dir, f))
    ]
    if not slide_folders:
        raise FileNotFoundError(f"No slide folders found in {seg_dir}")

    print(f"\nFound {len(slide_folders)} slide folder(s)")

    # ── Main processing loop ──────────────────────────────────────────────
    all_processed_tissues: Dict[str, list] = {}
    overlay_mapping: Dict[str, dict]        = {}
    processing_stats: List[dict]            = []

    for slide_folder in slide_folders:
        slide_dir = os.path.join(seg_dir, slide_folder)
        print(f"\n{'=' * 80}")
        print(f"Slide folder: {slide_folder}")
        print("=" * 80)

        csv_files = [f for f in os.listdir(slide_dir) if f.endswith("mesmer_result.csv")]
        if not csv_files:
            print(f"  No mesmer_result.csv found — skipping")
            continue

        # Load overlay pickle files for visualisation
        for pf in [f for f in os.listdir(slide_dir) if f.endswith("_seg_output.pickle")]:
            try:
                with open(os.path.join(slide_dir, pf), "rb") as fh:
                    data = pickle.load(fh)
                image_id = pf.replace("_seg_output.pickle", "")
                overlay_mapping[image_id] = data
                overlay_mapping[image_id + "_"] = data
                print(f"  Loaded overlay: {image_id}")
            except Exception as e:
                print(f"  Error loading overlay {pf}: {e}")

        for csv_file in csv_files:
            file_path = os.path.join(slide_dir, csv_file)
            print(f"\n--- {csv_file} ---")

            tissue_stats: dict = {
                "filename":    csv_file,
                "slide_folder": slide_folder,
                "timestamp":   datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            }

            try:
                df = pd.read_csv(file_path)
                df["slide_folder"] = slide_folder

                base_name = csv_file.replace("_mesmer_result.csv", "").replace("mesmer_result.csv", "")
                df["image_ID"] = base_name
                image_id  = base_name
                condition = detect_condition(image_id, conditions)
                df["condition"] = condition

                tissue_id = extract_tissue_identifier(image_id, conditions)
                tissue_stats.update({
                    "image_id":       image_id,
                    "condition":      condition,
                    "tissue_id":      tissue_id,
                    "original_cells": df.shape[0],
                })

                print(f"  image_ID:  {image_id}")
                print(f"  tissue_id: {tissue_id}")
                print(f"  condition: {condition}")

                # ── Skip if already processed ───────────────────────────
                if _processed_files_exist(tissue_id, condition, tissues_dir):
                    print(f"  ⏭️  Already processed — loading from disk")
                    tissue_stats["status"] = "SKIPPED - Already processed"
                    processing_stats.append(tissue_stats)

                    df_clean = pd.read_csv(
                        os.path.join(tissues_dir, f"{tissue_id}_{condition}.csv")
                    )
                    col_num = _get_last_marker_col(df_clean, last_marker)

                    all_processed_tissues.setdefault(tissue_id, []).append({
                        "data":      df_clean,
                        "condition": condition,
                        "image_id":  image_id,
                        "col_num":   col_num,
                        "skipped":   True,
                    })
                    continue

                # ── QC: size + DAPI filter ──────────────────────────────
                area_thresh = np.percentile(df["area"], size_pct)
                dapi_thresh = np.percentile(df["DAPI"],  dapi_pct)
                print(f"  {size_pct}% thresholds — area: {area_thresh:.2f}, DAPI: {dapi_thresh:.2f}")

                tissue_stats["area_threshold"] = area_thresh
                tissue_stats["dapi_threshold"] = dapi_thresh

                df_filt = sp.pp.filter_data(
                    df,
                    nuc_thres=dapi_thresh,
                    size_thres=area_thresh,
                    nuc_marker="DAPI",
                    cell_size="area",
                    log_scale=False,
                )
                n_after_filt = df_filt.shape[0]
                n_removed_filt = df.shape[0] - n_after_filt
                print(f"  After filter: {n_after_filt} cells  "
                      f"(removed {n_removed_filt}, "
                      f"{n_removed_filt / df.shape[0] * 100:.2f}%)")

                tissue_stats.update({
                    "cells_after_size_dapi_filter": n_after_filt,
                    "cells_removed_by_filter":      n_removed_filt,
                    "percent_removed_by_filter":    n_removed_filt / df.shape[0] * 100,
                })

                # ── Normalisation ───────────────────────────────────────
                print("  Normalising (z-score)…")
                df_norm = sp.pp.format(
                    data=df_filt,
                    list_out=[
                        "eccentricity", "perimeter", "convex_area",
                        "axis_major_length", "axis_minor_length", "label",
                    ],
                    list_keep=["DAPI", "x", "y", "area", "image_ID", "condition", "slide_folder"],
                    method="zscore",
                )

                col_num = _get_last_marker_col(df_norm, last_marker)

                # ── Noise removal ───────────────────────────────────────
                print(f"  Detecting noise cutoffs (last marker col: {col_num})…")
                z_count_cut, z_sum_cut = auto_detect_cutoffs(
                    df_norm, col_num,
                    cut_off=noise_cut_off,
                    count_bin=noise_bins,
                )
                print(f"  Cutoffs — z_count: {z_count_cut:.2f}, z_sum: {z_sum_cut:.2f}")

                tissue_stats["z_count_cutoff"] = z_count_cut
                tissue_stats["z_sum_cutoff"]   = z_sum_cut

                df_clean, _ = sp.pp.remove_noise(
                    df=df_norm,
                    col_num=col_num,
                    z_count_thres=z_count_cut,
                    z_sum_thres=z_sum_cut,
                )
                n_after_noise  = df_clean.shape[0]
                n_removed_noise = df_norm.shape[0] - n_after_noise
                print(f"  After noise removal: {n_after_noise} cells  "
                      f"(removed {n_removed_noise}, "
                      f"{n_removed_noise / df_norm.shape[0] * 100:.2f}%)")

                total_removed = df.shape[0] - n_after_noise
                tissue_stats.update({
                    "cells_after_noise_removal":  n_after_noise,
                    "cells_removed_by_noise":     n_removed_noise,
                    "percent_removed_by_noise":   n_removed_noise / df_norm.shape[0] * 100,
                    "total_cells_removed":        total_removed,
                    "total_percent_removed":      total_removed / df.shape[0] * 100,
                    "final_cells":                n_after_noise,
                    "status":                     "PROCESSED",
                })
                processing_stats.append(tissue_stats)

                all_processed_tissues.setdefault(tissue_id, []).append({
                    "data":      df_clean,
                    "condition": condition,
                    "image_id":  image_id,
                    "col_num":   col_num,
                    "skipped":   False,
                })
                print(f"  ✓ Processed {image_id}")

            except Exception as exc:
                print(f"  ✗ Error: {exc}")
                tissue_stats["status"] = f"ERROR: {exc}"
                processing_stats.append(tissue_stats)
                import traceback; traceback.print_exc()

    # ── Save per-condition + combined files ───────────────────────────────
    print(f"\n{'=' * 80}\nSAVING PROCESSED TISSUES\n{'=' * 80}")

    for tissue_id, tissue_list in all_processed_tissues.items():
        print(f"\n--- {tissue_id} ---")
        any_new = any(not t.get("skipped", False) for t in tissue_list)

        # Per-condition
        for tinfo in tissue_list:
            if tinfo.get("skipped"):
                print(f"  ⏭️  {tinfo['condition']} already on disk")
                continue
            cond     = tinfo["condition"]
            df_cond  = tinfo["data"]
            col_num  = tinfo["col_num"]

            csv_out = os.path.join(tissues_dir, f"{tissue_id}_{cond}.csv")
            df_cond.to_csv(csv_out, index=False)

            h5ad_out = os.path.join(tissues_dir, f"{tissue_id}_{cond}.h5ad")
            adata = sp.hf.make_anndata(df_nn=df_cond, col_sum=col_num, nonFuncAb_list=[])
            adata.write_h5ad(h5ad_out)
            print(f"  ✓ Saved {cond}: {df_cond.shape[0]} cells")

        # Combined — scan disk for all condition CSVs (handles partial runs)
        combined_csv = os.path.join(tissues_dir, f"{tissue_id}_combined_all_conditions.csv")
        combined_exists = os.path.exists(combined_csv)
        unique_conditions = list({t["condition"] for t in tissue_list})

        if combined_exists and not any_new and len(unique_conditions) == 1:
            print(f"  ⏭️  Combined already exists, no new data")
            continue

        disk_dfs = []
        for cond in conditions + ["Unknown"]:
            cond_csv = os.path.join(tissues_dir, f"{tissue_id}_{cond}.csv")
            if os.path.exists(cond_csv):
                df_c = pd.read_csv(cond_csv)
                df_c["condition"] = cond   # ensure correct label
                disk_dfs.append(df_c)
                print(f"  + {cond}: {len(df_c)} cells")

        if not disk_dfs:
            print(f"  ⚠️  No condition files found on disk — skipping combined")
            continue

        combined_df = pd.concat(disk_dfs, ignore_index=True)
        combined_df.to_csv(combined_csv, index=False)

        col_num = tissue_list[0]["col_num"]
        adata_combined = sp.hf.make_anndata(
            df_nn=combined_df, col_sum=col_num, nonFuncAb_list=[]
        )
        adata_combined_path = os.path.join(
            tissues_dir, f"{tissue_id}_combined_all_conditions.h5ad"
        )
        adata_combined.write_h5ad(adata_combined_path)
        print(f"  ✓ Combined: {len(combined_df)} cells → {adata_combined_path}")

    # ── Marker visualisations ─────────────────────────────────────────────
    print(f"\n{'=' * 80}\nGENERATING VISUALISATIONS\n{'=' * 80}")

    for tissue_id, tissue_list in all_processed_tissues.items():
        for tinfo in tissue_list:
            image_id  = tinfo["image_id"]
            condition = tinfo["condition"]
            df_tissue = tinfo["data"]
            safe_id   = tissue_id.replace("/", "_").replace(" ", "_")

            viz_path = os.path.join(viz_dir, f"{safe_id}_{condition}_all_markers.png")
            if os.path.exists(viz_path):
                print(f"  ⏭️  Viz already exists: {safe_id}_{condition}")
                continue

            if image_id not in overlay_mapping:
                print(f"  ⚠️  No overlay data for {image_id} — skipping viz")
                continue

            overlay = overlay_mapping[image_id]
            if "img" not in overlay:
                print(f"  ⚠️  No 'img' key in overlay for {image_id} — skipping viz")
                continue

            METADATA = {
                "label", "area", "x", "y", "slide_folder", "image_ID", "condition",
                "eccentricity", "perimeter", "convex_area",
                "axis_major_length", "axis_minor_length", "DAPI",
            }
            all_markers = [c for c in df_tissue.columns if c not in METADATA]
            if not all_markers:
                continue

            n_cols = 5
            n_rows = int(np.ceil(len(all_markers) / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 4, n_rows * 4))
            fig.suptitle(f"{tissue_id} – {condition} – All Markers", fontsize=16)
            axes_flat = np.array(axes).flatten()

            for i, marker in enumerate(all_markers):
                ax = axes_flat[i]
                ax.imshow(overlay["img"], cmap="gray")
                sc = ax.scatter(
                    df_tissue["x"], df_tissue["y"],
                    c=df_tissue[marker], s=2, cmap="viridis", alpha=0.7,
                    vmin=df_tissue[marker].quantile(0.01),
                    vmax=df_tissue[marker].quantile(0.99),
                )
                plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
                ax.set_title(marker, fontsize=10)
                ax.axis("off")

            for j in range(i + 1, len(axes_flat)):
                axes_flat[j].axis("off")

            plt.tight_layout()
            plt.savefig(viz_path, dpi=150, bbox_inches="tight")
            plt.close()
            print(f"  ✓ Saved viz: {viz_path}")

    # ── Save processing statistics ────────────────────────────────────────
    stats_df   = pd.DataFrame(processing_stats)
    stats_path = os.path.join(log_dir, f"processing_stats_{timestamp}.csv")
    stats_df.to_csv(stats_path, index=False)

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n{'=' * 80}\nSUMMARY\n{'=' * 80}")
    print(f"Unique tissues processed: {len(all_processed_tissues)}")
    if len(stats_df) > 0:
        done  = stats_df[stats_df["status"] == "PROCESSED"]
        skip  = stats_df[stats_df.get("status", pd.Series()).str.contains("SKIPPED", na=False)]
        errs  = stats_df[stats_df.get("status", pd.Series()).str.contains("ERROR",   na=False)]
        print(f"  Processed: {len(done)}   Skipped: {len(skip)}   Errors: {len(errs)}")
        if len(done) > 0:
            print(f"  Avg total removal: {done['total_percent_removed'].mean():.2f}%")

    print(f"\nOutput locations:")
    print(f"  Tissues:  {tissues_dir}")
    print(f"  Viz:      {viz_dir}")
    print(f"  Logs:     {log_dir}")
    print(f"\nProcessing complete: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 80)

    # Restore stdout
    dual_log.close()
    sys.stdout = dual_log.terminal

    return {
        "all_processed_tissues": all_processed_tissues,
        "processing_stats":      processing_stats,
        "output_dir":            out_dir,
    }


# ─────────────────────────────────────────────────────────────────────────────
# QUPATH EXPORT (optional second step)
# ─────────────────────────────────────────────────────────────────────────────

def _classify_cells(
    df_raw: pd.DataFrame,
    included_keys: set,
    area_thresh: float,
    dapi_thresh: float,
) -> pd.Series:
    """
    Label each raw cell with a QC classification:
      Included | Excl_SmallArea | Excl_LowDAPI |
      Excl_SmallArea_LowDAPI | Excl_Noise

    Cells are matched to processed output by (x, y) coordinates
    rounded to 2 decimal places.
    """
    coord_key = list(zip(
        df_raw["x"].round(2).astype(str),
        df_raw["y"].round(2).astype(str),
    ))

    small_area = df_raw["area"] < area_thresh
    low_dapi   = df_raw["DAPI"] < dapi_thresh

    labels = []
    for i, key in enumerate(coord_key):
        if key in included_keys:
            labels.append("Included")
        elif small_area.iloc[i] and low_dapi.iloc[i]:
            labels.append("Excl_SmallArea_LowDAPI")
        elif small_area.iloc[i]:
            labels.append("Excl_SmallArea")
        elif low_dapi.iloc[i]:
            labels.append("Excl_LowDAPI")
        else:
            labels.append("Excl_Noise")
    return pd.Series(labels, index=df_raw.index)


def run_qupath_export(cfg: dict) -> str:
    """
    Generate per-image QuPath TSV files and QC visualisations.

    Requires run_preprocessing() to have already been executed
    (reads per-condition CSVs from individual_processed_data/).

    Returns
    -------
    Path to qupath_exports/ directory.
    """
    conditions = cfg["experiment"]["conditions"]
    seg_dir    = cfg["paths"]["segmentation_results_dir"]
    base_out   = cfg["paths"]["output_dir"]

    pp_cfg     = cfg.get("preprocessing", {})
    last_marker = pp_cfg.get("last_marker", None)
    size_pct    = pp_cfg.get("qc_filter", {}).get("size_percentile", 1)
    dapi_pct    = pp_cfg.get("qc_filter", {}).get("dapi_percentile", 1)

    tissues_dir = os.path.join(base_out, "combined_processed_data", "individual_processed_data")
    qupath_dir  = os.path.join(base_out, "qupath_exports")
    qcviz_dir   = os.path.join(qupath_dir, "qc_visualizations")
    os.makedirs(qupath_dir, exist_ok=True)
    os.makedirs(qcviz_dir,  exist_ok=True)

    COLOR_MAP = {
        "Included":             "#00D200",
        "Excl_SmallArea":       "#FF3C3C",
        "Excl_LowDAPI":         "#FFA500",
        "Excl_SmallArea_LowDAPI": "#B400B4",
        "Excl_Noise":           "#1E90FF",
    }

    slide_folders = [
        f for f in os.listdir(seg_dir)
        if os.path.isdir(os.path.join(seg_dir, f))
    ]

    summary_records = []

    for slide_folder in slide_folders:
        slide_dir = os.path.join(seg_dir, slide_folder)
        csv_files = [f for f in os.listdir(slide_dir) if f.endswith("mesmer_result.csv")]

        # Load overlay images for this slide
        overlay_mapping: Dict[str, dict] = {}
        for pf in [f for f in os.listdir(slide_dir) if f.endswith("_seg_output.pickle")]:
            try:
                with open(os.path.join(slide_dir, pf), "rb") as fh:
                    data = pickle.load(fh)
                image_id = pf.replace("_seg_output.pickle", "")
                overlay_mapping[image_id] = data
            except Exception:
                pass

        for csv_file in csv_files:
            df_raw = pd.read_csv(os.path.join(slide_dir, csv_file))
            image_id  = csv_file.replace("_mesmer_result.csv", "").replace("mesmer_result.csv", "")
            condition = detect_condition(image_id, conditions)
            tissue_id = extract_tissue_identifier(image_id, conditions)

            df_raw["image_ID"]  = image_id
            df_raw["condition"] = condition

            # Load processed (included) cells
            proc_csv = os.path.join(tissues_dir, f"{tissue_id}_{condition}.csv")
            if not os.path.exists(proc_csv):
                print(f"⚠️  Processed file not found: {proc_csv}  — run preprocessing first")
                continue

            df_proc = pd.read_csv(proc_csv)
            area_thresh = np.percentile(df_raw["area"], size_pct)
            dapi_thresh = np.percentile(df_raw["DAPI"],  dapi_pct)

            included_keys = set(zip(
                df_proc["x"].round(2).astype(str),
                df_proc["y"].round(2).astype(str),
            ))

            df_raw["classification"] = _classify_cells(
                df_raw, included_keys, area_thresh, dapi_thresh
            )
            df_raw["area_threshold"] = area_thresh
            df_raw["dapi_threshold"] = dapi_thresh
            df_raw["tissue_id"]      = tissue_id

            # Build QuPath TSV
            import math
            df_raw["radius"]    = np.sqrt(df_raw["area"] / math.pi).clip(lower=3.0)
            df_raw["centroid_x"] = df_raw["x"]
            df_raw["centroid_y"] = df_raw["y"]
            df_raw["roi_x"]      = df_raw["x"] - df_raw["radius"]
            df_raw["roi_y"]      = df_raw["y"] - df_raw["radius"]
            df_raw["roi_width"]  = 2 * df_raw["radius"]
            df_raw["roi_height"] = 2 * df_raw["radius"]

            tsv_cols = [
                "centroid_x", "centroid_y", "roi_x", "roi_y",
                "roi_width", "roi_height", "area", "DAPI",
                "condition", "image_ID", "tissue_id",
                "classification", "area_threshold", "dapi_threshold",
            ]
            safe_id  = image_id.replace("/", "_").replace(" ", "_")
            tsv_path = os.path.join(qupath_dir, f"{safe_id}_qupath_cells.tsv")
            df_raw[tsv_cols].to_csv(tsv_path, sep="\t", index=False)

            counts = df_raw["classification"].value_counts().to_dict()
            print(f"  ✓ QuPath TSV: {os.path.basename(tsv_path)}  |  "
                  + "  ".join(f"{k}: {v}" for k, v in counts.items()))

            summary_records.append({
                "tissue_id": tissue_id,
                "condition": condition,
                "total_raw": len(df_raw),
                "tsv_file":  os.path.basename(tsv_path),
                **counts,
            })

            # ── QC visualisations ──────────────────────────────────────
            overlay_img = overlay_mapping.get(image_id, {}).get("img", None)

            # Panel 1: spatial inclusion/exclusion
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
            fig.suptitle(f"{tissue_id} – {condition} QC", fontsize=14)
            for ax, subset, title in [
                (ax1, df_raw[df_raw["classification"] == "Included"], "Included"),
                (ax2, df_raw[df_raw["classification"] != "Included"], "Excluded"),
            ]:
                if overlay_img is not None:
                    ax.imshow(overlay_img, cmap="gray")
                colors = [COLOR_MAP.get(c, "grey") for c in subset["classification"]]
                ax.scatter(subset["x"], subset["y"], c=colors, s=1, alpha=0.6)
                ax.set_title(title, fontsize=12)
                ax.axis("off")
            plt.tight_layout()
            plt.savefig(os.path.join(qcviz_dir, f"{safe_id}_QC_spatial.png"), dpi=120, bbox_inches="tight")
            plt.close()

            # Panel 2: bar chart
            fig, ax = plt.subplots(figsize=(8, 5))
            labels_order = list(COLOR_MAP.keys())
            vals = [counts.get(l, 0) for l in labels_order]
            colors_bar = [COLOR_MAP[l] for l in labels_order]
            bars = ax.bar(labels_order, vals, color=colors_bar)
            total = sum(vals)
            for bar, val in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + total * 0.005,
                        f"{val}\n({val/total*100:.1f}%)", ha="center", va="bottom", fontsize=9)
            ax.set_ylabel("Cell count")
            ax.set_title(f"{tissue_id} – {condition} QC breakdown")
            plt.xticks(rotation=15, ha="right")
            plt.tight_layout()
            plt.savefig(os.path.join(qcviz_dir, f"{safe_id}_QC_barchart.png"), dpi=120, bbox_inches="tight")
            plt.close()

            # Panel 3: DAPI vs Area
            fig, ax = plt.subplots(figsize=(8, 6))
            for cls, color in COLOR_MAP.items():
                sub = df_raw[df_raw["classification"] == cls]
                if len(sub):
                    ax.scatter(sub["area"], sub["DAPI"], c=color, s=1,
                               alpha=0.4, label=cls)
            ax.axvline(area_thresh, color="black", linestyle="--", linewidth=1, label="area threshold")
            ax.axhline(dapi_thresh, color="gray",  linestyle="--", linewidth=1, label="DAPI threshold")
            ax.set_xlabel("Cell area")
            ax.set_ylabel("DAPI intensity")
            ax.set_title(f"{tissue_id} – {condition} DAPI vs Area")
            ax.legend(markerscale=5, fontsize=8, loc="upper right")
            plt.tight_layout()
            plt.savefig(os.path.join(qcviz_dir, f"{safe_id}_QC_dapi_vs_area.png"), dpi=120, bbox_inches="tight")
            plt.close()

    # Summary CSV
    if summary_records:
        pd.DataFrame(summary_records).to_csv(
            os.path.join(qupath_dir, "qupath_export_summary.csv"), index=False
        )

    print(f"\n✓ QuPath exports saved to: {qupath_dir}")
    return qupath_dir
