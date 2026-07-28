"""
spatia.analysis.functional
==========================
Standalone functional marker analysis step.

Reads the per-image outputs already written by the triads step:
  {output_dir}/{image_id}_cells_with_triad_flags.csv  — cell-level data + triad membership
  {output_dir}/{image_id}_triad_pairs.csv             — pairwise distances for each triad

Recomputes in-triad membership at the functional reporting radius (which can be
tighter than the search radius used during detection) and then, for each configured
cell type, compares marker expression:

  1. In-triad vs not-in-triad  (violin + Mann-Whitney U, all conditions pooled)
  2. Condition A vs condition B within in-triad cells only  (CLR vs DII, etc.)

Outputs land in {output_dir}/functional_analysis/:
  functional_{CellType}_intriad_vs_out.png
  functional_{CellType}_condition_compare.png
  functional_marker_summary.csv

Usage (standalone pipeline step)
---------------------------------
    python run_pipeline.py --config experiments/crc_tma.yaml \
        --steps functional

Config shape (crc_tma.yaml)
----------------------------
    analysis:
      functional:
        enabled: true
        report_radius_um: 20.0
        markers:
          GzmB:   "Granzyme B - cytotoxicity:Cyc_13_ch_2"
          PD-1:   "PD-1 - checkpoint:Cyc_12_ch_4"
          ...
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu


# ── Config helpers (mirrors triads.py) ───────────────────────────────────────

def _get_condition(image_id: str, image_condition_map: dict, conditions: list) -> str:
    if image_id in image_condition_map:
        return image_condition_map[image_id]
    for cond in conditions:
        if image_id.upper().startswith(cond.upper()):
            return cond
    return "Unknown"


# ── Plotting helpers ──────────────────────────────────────────────────────────

def _violin_pair(ax, vals_a, vals_b, label_a, label_b, color_a, color_b, title):
    """Draw a two-group violin on ax. Handles empty groups gracefully."""
    if len(vals_a) >= 3 and len(vals_b) >= 3:
        stat, pval = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
    else:
        stat, pval = np.nan, np.nan

    if len(vals_a) > 0 and len(vals_b) > 0:
        parts = ax.violinplot([vals_a, vals_b], positions=[0, 1],
                              showmedians=True, showextrema=False)
        for pc, color in zip(parts["bodies"], [color_a, color_b]):
            pc.set_facecolor(color)
            pc.set_alpha(0.75)
        parts["cmedians"].set_color("black")
        parts["cmedians"].set_linewidth(2)

    sig = ("***" if pval < 0.001 else "**" if pval < 0.01
           else "*" if pval < 0.05 else "ns")
    pval_str = f"{pval:.2e}" if not np.isnan(pval) else "n/a"
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{label_a}\n(n={len(vals_a):,})",
                        f"{label_b}\n(n={len(vals_b):,})"], fontsize=8)
    ax.set_title(f"{title}\n{sig}  p={pval_str}", fontsize=9, fontweight="bold")
    ax.set_ylabel("Expression (intensity)", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return stat, pval, sig


# ── Main entry point ──────────────────────────────────────────────────────────

def run_functional_analysis(cfg: dict) -> None:
    """
    Standalone functional marker analysis.  Called by run_pipeline.py.
    Reads triads step outputs from output_dir — run after the triads step.
    """
    exp_name     = cfg["experiment"]["name"]
    conditions   = cfg["experiment"]["conditions"]
    img_cond_map = cfg["experiment"].get("image_condition_map", {})
    output_dir   = cfg["paths"]["output_dir"]

    func_cfg = cfg.get("analysis", {}).get("functional", {})
    if not func_cfg.get("enabled", False):
        print("[functional] disabled in config — skipping.")
        return

    markers        = func_cfg.get("markers", {})
    func_radius    = func_cfg.get("report_radius_um", 20.0)
    t_cfg          = cfg.get("analysis", {}).get("triad", {})
    anchor_type    = t_cfg.get("anchor_type",    "anchor")
    partner1_type  = t_cfg.get("partner_type_1", "partner1")
    partner2_type  = t_cfg.get("partner_type_2", "partner2")
    anchor_name    = t_cfg.get("anchor_name",    anchor_type)
    partner1_name  = t_cfg.get("partner_1_name", partner1_type)
    partner2_name  = t_cfg.get("partner_2_name", partner2_type)

    if not markers:
        print("[functional] No markers configured — add analysis.functional.markers to YAML.")
        return

    func_out = os.path.join(output_dir, "functional_analysis")
    os.makedirs(func_out, exist_ok=True)

    print(f"[functional] Experiment  : {exp_name}")
    print(f"[functional] Conditions  : {conditions}")
    print(f"[functional] Func radius : {func_radius} µm")
    print(f"[functional] Markers     : {list(markers.keys())}")
    print(f"[functional] Input dir   : {output_dir}\n")

    # ── Discover per-image cells files ────────────────────────────────────────
    cell_files = sorted(glob.glob(
        os.path.join(output_dir, "*_cells_with_triad_flags.csv")
    ))
    if not cell_files:
        print(f"[functional] ⚠️  No *_cells_with_triad_flags.csv found in {output_dir}")
        print("[functional]    Run the triads step first.")
        return
    print(f"[functional] Found {len(cell_files)} cell files.\n")

    # ── Load and accumulate per-image data ───────────────────────────────────
    all_cell_dfs = []

    for cell_path in cell_files:
        basename = os.path.basename(cell_path)
        # e.g. CLR_reg001_A_cells_with_triad_flags.csv → CLR_reg001_A
        image_id = basename.replace("_cells_with_triad_flags.csv", "")
        condition = _get_condition(image_id, img_cond_map, conditions)

        # Load cell-level data (contains all original columns if KEEP_INTENSITIES=True)
        try:
            df = pd.read_csv(cell_path)
        except Exception as e:
            print(f"  ⚠️  Cannot read {basename}: {e} — skipping.")
            continue

        if "cell_id" not in df.columns:
            df.insert(0, "cell_id", df.index.astype(str))

        # Load corresponding triad pairs (may not exist if no triads were found)
        pairs_path = os.path.join(output_dir, f"{image_id}_triad_pairs.csv")
        triad_pairs = pd.DataFrame()
        if os.path.exists(pairs_path):
            try:
                triad_pairs = pd.read_csv(pairs_path)
            except Exception:
                pass

        # Recompute in-triad flags at func_radius (tighter than search radius)
        if not triad_pairs.empty and "dist_anchor_p1_um" in triad_pairs.columns:
            pairs_report = triad_pairs[
                triad_pairs[["dist_anchor_p1_um", "dist_anchor_p2_um"]].max(axis=1)
                <= func_radius
            ]
        else:
            pairs_report = pd.DataFrame()

        df["_in_triad_anchor"]   = False
        df["_in_triad_partner1"] = False
        df["_in_triad_partner2"] = False

        if not pairs_report.empty:
            a_ids  = set(pairs_report["anchor_cell_id"].astype(str))
            p1_ids = set(pairs_report["partner1_cell_id"].astype(str))
            p2_ids = set(pairs_report["partner2_cell_id"].astype(str))
            cid    = df["cell_id"].astype(str)
            df["_in_triad_anchor"]   = cid.isin(a_ids)
            df["_in_triad_partner1"] = cid.isin(p1_ids)
            df["_in_triad_partner2"] = cid.isin(p2_ids)

        df["image_id"]  = image_id
        df["condition"] = condition
        all_cell_dfs.append(df)

    if not all_cell_dfs:
        print("[functional] ⚠️  No cell data loaded — check that the triads step ran first.")
        return

    all_cells = pd.concat(all_cell_dfs, ignore_index=True)
    print(f"[functional] Aggregated {len(all_cells):,} cells from {len(all_cell_dfs)} images\n")

    # ── Per-cell-type analysis ────────────────────────────────────────────────
    COND_COLORS = dict(zip(
        conditions,
        ["#2196F3", "#F44336", "#4CAF50", "#FF9800", "#9C27B0"],
    ))
    COLOR_OUT = "#1f77b4"   # blue — out of triad
    COLOR_IN  = "#ff7f0e"   # orange — in triad

    cell_role_map = [
        (anchor_type,   anchor_name,   "_in_triad_anchor"),
        (partner1_type, partner1_name, "_in_triad_partner1"),
        (partner2_type, partner2_name, "_in_triad_partner2"),
    ]

    all_summary = []

    for cell_type, display_name, flag_col in cell_role_map:
        cells = all_cells[all_cells["cell_type"] == cell_type].copy()
        if cells.empty:
            print(f"[functional] ⚠️  No {display_name} cells — skipping.")
            continue

        avail   = {n: c for n, c in markers.items() if c in cells.columns}
        missing = [n for n in markers if n not in avail]
        if missing:
            print(f"[functional]   Marker columns not found for {display_name}: {missing}")
            print(f"[functional]   (Re-run prepare_crc_data.py with KEEP_INTENSITIES=True, "
                  f"then re-run the triads step.)")
        if not avail:
            print(f"[functional] ⚠️  No usable marker columns for {display_name} — skipping.")
            continue

        in_triad     = cells[cells[flag_col] == True]
        not_in_triad = cells[cells[flag_col] == False]
        n_in, n_out  = len(in_triad), len(not_in_triad)
        print(f"[functional] {display_name}: {n_in:,} in-triad  |  {n_out:,} out-of-triad")

        n_m = len(avail)

        # ── Plot 1: in-triad vs out-of-triad ─────────────────────────────────
        fig1, axes1 = plt.subplots(1, n_m, figsize=(max(4 * n_m, 8), 6))
        axes1 = [axes1] if n_m == 1 else list(axes1)

        for ax, (marker_name, col) in zip(axes1, avail.items()):
            vals_in  = in_triad[col].dropna().values.astype(float)
            vals_out = not_in_triad[col].dropna().values.astype(float)

            stat, pval, sig = _violin_pair(
                ax, vals_out, vals_in,
                label_a=f"Out of triad", label_b=f"In triad",
                color_a=COLOR_OUT, color_b=COLOR_IN,
                title=marker_name,
            )
            all_summary.append({
                "cell_type":     display_name,
                "comparison":    "in_triad_vs_out",
                "marker":        marker_name,
                "marker_col":    col,
                "group_A":       "in_triad",
                "group_B":       "not_in_triad",
                "n_A":           len(vals_in),
                "n_B":           len(vals_out),
                "mean_A":        float(np.mean(vals_in))    if len(vals_in)  else np.nan,
                "mean_B":        float(np.mean(vals_out))   if len(vals_out) else np.nan,
                "median_A":      float(np.median(vals_in))  if len(vals_in)  else np.nan,
                "median_B":      float(np.median(vals_out)) if len(vals_out) else np.nan,
                "mannwhitney_U": stat,
                "p_value":       pval,
                "significance":  sig,
            })

        fig1.suptitle(
            f"{display_name}: In-triad vs Out-of-triad Marker Expression\n"
            f"(triad threshold = {func_radius} µm)",
            fontsize=12, fontweight="bold",
        )
        plt.tight_layout()
        fig1.savefig(
            os.path.join(func_out, f"functional_{display_name}_intriad_vs_out.png"),
            dpi=150, bbox_inches="tight",
        )
        plt.close(fig1)

        # ── Plot 2: condition comparison within in-triad cells ────────────────
        if len(conditions) >= 2 and "condition" in cells.columns and n_in >= 6:
            cond_sub    = {c: in_triad[in_triad["condition"] == c] for c in conditions}
            valid_conds = [c for c in conditions if len(cond_sub[c]) >= 3]

            if len(valid_conds) >= 2:
                fig2, axes2 = plt.subplots(1, n_m, figsize=(max(4 * n_m, 8), 6))
                axes2 = [axes2] if n_m == 1 else list(axes2)

                cA, cB = valid_conds[0], valid_conds[1]
                for ax, (marker_name, col) in zip(axes2, avail.items()):
                    vA = cond_sub[cA][col].dropna().values.astype(float)
                    vB = cond_sub[cB][col].dropna().values.astype(float)

                    stat2, pval2, sig2 = _violin_pair(
                        ax, vA, vB,
                        label_a=cA, label_b=cB,
                        color_a=COND_COLORS.get(cA, "#888"),
                        color_b=COND_COLORS.get(cB, "#888"),
                        title=marker_name,
                    )
                    all_summary.append({
                        "cell_type":     display_name,
                        "comparison":    f"{cA}_vs_{cB}_in_triad",
                        "marker":        marker_name,
                        "marker_col":    col,
                        "group_A":       cA,
                        "group_B":       cB,
                        "n_A":           len(vA),
                        "n_B":           len(vB),
                        "mean_A":        float(np.mean(vA))    if len(vA) else np.nan,
                        "mean_B":        float(np.mean(vB))    if len(vB) else np.nan,
                        "median_A":      float(np.median(vA))  if len(vA) else np.nan,
                        "median_B":      float(np.median(vB))  if len(vB) else np.nan,
                        "mannwhitney_U": stat2,
                        "p_value":       pval2,
                        "significance":  sig2,
                    })

                fig2.suptitle(
                    f"{display_name}: {cA} vs {cB} — In-triad cells only\n"
                    f"(triad threshold = {func_radius} µm)",
                    fontsize=12, fontweight="bold",
                )
                plt.tight_layout()
                fig2.savefig(
                    os.path.join(func_out, f"functional_{display_name}_condition_compare.png"),
                    dpi=150, bbox_inches="tight",
                )
                plt.close(fig2)

    # ── Summary CSV ───────────────────────────────────────────────────────────
    if all_summary:
        summary_df = pd.DataFrame(all_summary)
        n_tests = len(summary_df)
        bonf_alpha = 0.05 / max(n_tests, 1)
        summary_df["bonferroni_alpha"] = round(bonf_alpha, 8)
        summary_df["sig_bonferroni"]   = summary_df["p_value"] < bonf_alpha
        out_csv = os.path.join(func_out, "functional_marker_summary.csv")
        summary_df.to_csv(out_csv, index=False)
        sig_hits = summary_df["sig_bonferroni"].sum()
        print(f"\n[functional] ✓  {n_tests} tests  |  Bonferroni α={bonf_alpha:.3e}  "
              f"|  {sig_hits} hits after correction")
        print(f"[functional]    Outputs → {func_out}")
    else:
        print("[functional] ⚠️  No results produced.")
