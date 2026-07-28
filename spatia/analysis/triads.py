"""
spatia.analysis.triads
======================
Generalizable triad detection for spatial proteomics data.

All parameters come from a config dict (loaded from config_example.yaml).
No hardcoded paths, conditions, cell types, or pixel sizes.

Usage
-----
    from spatia.analysis.triads import run_triad_analysis
    run_triad_analysis(config)   # config is a dict loaded from YAML

The config dict shape mirrors config_example.yaml.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from itertools import combinations, permutations


# ── Config helpers ────────────────────────────────────────────────────────────

def _get_condition_areas(cfg: dict) -> dict:
    """
    Returns a flat {condition: area_um2} dict.
    Handles both flat and timepoint-nested config formats.
    """
    raw = cfg["imaging"]["condition_areas_um2"]
    # If values are dicts, it's timepoint-nested — flatten by summing areas.
    if raw and isinstance(next(iter(raw.values())), dict):
        flat = {}
        for tp_areas in raw.values():
            for cond, area in tp_areas.items():
                flat[cond] = flat.get(cond, 0) + area
        return flat
    return dict(raw)


def _get_condition(image_id: str, image_condition_map: dict, conditions: list) -> str:
    """
    Determine condition from image_id.
    Priority: explicit map → prefix match → 'Unknown'.
    """
    if image_id in image_condition_map:
        return image_condition_map[image_id]
    for cond in conditions:
        if image_id.upper().startswith(cond.upper()):
            return cond
    return "Unknown"


# ── Core detection ────────────────────────────────────────────────────────────

def find_triads(
    df: pd.DataFrame,
    anchor_type: str,
    partner1_type: str,
    partner2_type: str,
    radius_px: float,
    tree: cKDTree,
    microns_per_pixel: float,
) -> pd.DataFrame:
    """
    Find ALL valid triads (anchor + partner1 + partner2) within radius_px.

    For each anchor cell, every combination of partner1 × partner2
    within the radius is recorded (not just the nearest pair).

    Returns a DataFrame of triad records, one row per triad.
    """
    anchors = df[df["cell_type"] == anchor_type]
    if anchors.empty:
        return pd.DataFrame()

    coords = df[["centroid_x", "centroid_y"]].values
    ct_arr = df["cell_type"].values
    cell_ids = (
        df["cell_id"].values
        if "cell_id" in df.columns
        else df.index.astype(str).values
    )

    combo_label = f"{anchor_type}__{partner1_type}__{partner2_type}"
    triad_records = []
    anchor_indices = anchors.index.tolist()
    anchor_coords = anchors[["centroid_x", "centroid_y"]].values

    neighbor_lists = tree.query_ball_point(anchor_coords, r=radius_px, workers=-1)

    for ai, (anchor_idx, neighbors_idx) in enumerate(zip(anchor_indices, neighbor_lists)):
        neighbors_idx = [n for n in neighbors_idx if n != anchor_idx]
        if not neighbors_idx:
            continue

        neighbor_types = ct_arr[neighbors_idx]
        p1_mask = neighbor_types == partner1_type
        p2_mask = neighbor_types == partner2_type

        if not p1_mask.any() or not p2_mask.any():
            continue

        ax, ay = anchor_coords[ai]
        p1_idxs = np.array(neighbors_idx)[p1_mask]
        p1_dists = np.linalg.norm(coords[p1_idxs] - [ax, ay], axis=1)
        p2_idxs = np.array(neighbors_idx)[p2_mask]
        p2_dists = np.linalg.norm(coords[p2_idxs] - [ax, ay], axis=1)

        for p1_idx, p1_dist in zip(p1_idxs, p1_dists):
            for p2_idx, p2_dist in zip(p2_idxs, p2_dists):
                if p1_idx == p2_idx:
                    continue
                p1p2_dist_px = np.linalg.norm(coords[p1_idx] - coords[p2_idx])
                triad_records.append({
                    "anchor_cell_id":    cell_ids[anchor_idx],
                    "anchor_x":          ax,
                    "anchor_y":          ay,
                    "anchor_type":       anchor_type,
                    "partner1_cell_id":  cell_ids[p1_idx],
                    "partner1_x":        coords[p1_idx, 0],
                    "partner1_y":        coords[p1_idx, 1],
                    "dist_anchor_p1_um": p1_dist * microns_per_pixel,
                    "partner1_type":     partner1_type,
                    "partner2_cell_id":  cell_ids[p2_idx],
                    "partner2_x":        coords[p2_idx, 0],
                    "partner2_y":        coords[p2_idx, 1],
                    "dist_anchor_p2_um": p2_dist * microns_per_pixel,
                    "dist_p1_p2_um":     p1p2_dist_px * microns_per_pixel,
                    "partner2_type":     partner2_type,
                    "triad_combo_label": combo_label,
                    "n_p1_neighbors":    int(p1_mask.sum()),
                    "n_p2_neighbors":    int(p2_mask.sum()),
                })

    return pd.DataFrame(triad_records)


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_triad_qc(df, triad_df, image_id, combo_label, save_path, radius_um,
                  report_radius_um=None):
    """
    Two-panel QC plot for a single image.
    Left : all cells, coloured by type, with triad anchors overlaid.
    Right: only triad cells (filtered to report_radius_um if set), with connecting lines.
    """
    # Use tab10 (10 vivid colours) cycling — much more visible than tab20
    COLOR_PALETTE = [
        "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
        "#a65628", "#f781bf", "#999999", "#00ced1", "#8b0000",
    ]
    ANCHOR_COLOR   = "#FFD700"   # gold — pops on any background
    PARTNER1_COLOR = "#1f77b4"   # bright blue
    PARTNER2_COLOR = "#d62728"   # vivid red

    # Filter triads to the reporting threshold for the right-hand panel.
    # Use max(anchor→p1, anchor→p2) — "both partners within X µm of DC".
    # The p1↔p2 distance is NOT included: two T cells can be far apart while
    # both sitting within contact range of the same DC.
    if report_radius_um is not None and not triad_df.empty:
        triad_df_report = triad_df[
            triad_df[["dist_anchor_p1_um", "dist_anchor_p2_um"]].max(axis=1)
            <= report_radius_um
        ].copy()
        panel_radius = report_radius_um
    else:
        triad_df_report = triad_df.copy()
        panel_radius = radius_um

    fig, axes = plt.subplots(1, 2, figsize=(20, 9))
    cell_types = sorted(df["cell_type"].unique())

    # ── Left panel: all cells coloured by type ────────────────────────────────
    for i, ct in enumerate(cell_types):
        sub = df[df["cell_type"] == ct]
        is_key = ct in (combo_label.split("__")[0], combo_label.split("__")[1], combo_label.split("__")[2])
        axes[0].scatter(
            sub["centroid_x"], sub["centroid_y"],
            s=10 if is_key else 6,
            alpha=0.75 if is_key else 0.45,
            linewidths=0,
            color=COLOR_PALETTE[i % len(COLOR_PALETTE)],
            label=f"{ct} ({len(sub):,})",
        )
    if not triad_df_report.empty:
        axes[0].scatter(
            triad_df_report["anchor_x"], triad_df_report["anchor_y"],
            s=60, color=ANCHOR_COLOR, edgecolors="black", linewidths=0.5,
            zorder=6, marker="*",
            label=f"Triad anchors @ {panel_radius} µm ({len(triad_df_report):,})",
        )

    axes[0].set_title(f"{image_id} — All cells\n(★ = triad anchors within {panel_radius} µm)", fontsize=10, fontweight="bold")
    axes[0].invert_yaxis()
    axes[0].set_xlabel("X (px)", fontsize=10)
    axes[0].set_ylabel("Y (px)", fontsize=10)
    axes[0].legend(fontsize=7, markerscale=2, ncol=2, framealpha=0.9,
                   loc="upper right", borderpad=0.5)
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    # ── Right panel: triad cells only ────────────────────────────────────────
    names = combo_label.split("__")
    anchor_lbl   = names[0] if len(names) > 0 else "Anchor"
    partner1_lbl = names[1] if len(names) > 1 else "Partner 1"
    partner2_lbl = names[2] if len(names) > 2 else "Partner 2"

    # Background: all cells dimmed (always rendered so axes scale correctly)
    triad_ids = set()
    if not triad_df_report.empty:
        triad_ids = (
            set(triad_df_report["anchor_cell_id"].astype(str))
            | set(triad_df_report["partner1_cell_id"].astype(str))
            | set(triad_df_report["partner2_cell_id"].astype(str))
        )
    bg = df[~df["cell_id"].astype(str).isin(triad_ids)]
    axes[1].scatter(bg["centroid_x"], bg["centroid_y"],
                    s=4, alpha=0.15, color="lightgray", linewidths=0, zorder=1)

    if not triad_df_report.empty:
        # Connecting lines
        sample = triad_df_report.sample(min(500, len(triad_df_report)), random_state=42)
        for _, row in sample.iterrows():
            axes[1].plot([row["anchor_x"], row["partner1_x"]], [row["anchor_y"], row["partner1_y"]],
                         color=PARTNER1_COLOR, alpha=0.5, linewidth=0.8, zorder=2)
            axes[1].plot([row["anchor_x"], row["partner2_x"]], [row["anchor_y"], row["partner2_y"]],
                         color=PARTNER2_COLOR, alpha=0.5, linewidth=0.8, zorder=2)

        # Triad cell dots — drawn last so they're on top
        axes[1].scatter(triad_df_report["partner1_x"], triad_df_report["partner1_y"],
                        s=35, color=PARTNER1_COLOR, alpha=0.85, linewidths=0.4,
                        edgecolors="white", zorder=4, label=f"{partner1_lbl}")
        axes[1].scatter(triad_df_report["partner2_x"], triad_df_report["partner2_y"],
                        s=35, color=PARTNER2_COLOR, alpha=0.85, linewidths=0.4,
                        edgecolors="white", zorder=4, label=f"{partner2_lbl}")
        axes[1].scatter(triad_df_report["anchor_x"], triad_df_report["anchor_y"],
                        s=80, color=ANCHOR_COLOR, edgecolors="black", linewidths=0.6,
                        zorder=5, marker="*", label=f"{anchor_lbl} (anchor)")

    axes[1].set_title(
        f"Triads: {combo_label}\n"
        f"n = {len(triad_df_report):,} triads  |  radius = {panel_radius} µm",
        fontsize=10, fontweight="bold"
    )
    axes[1].invert_yaxis()
    axes[1].set_xlabel("X (px)", fontsize=10)
    axes[1].set_ylabel("Y (px)", fontsize=10)
    axes[1].legend(fontsize=9, markerscale=1.5, framealpha=0.9)
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    fig.suptitle(f"{image_id}", fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(save_path, dpi=180, bbox_inches="tight")
    plt.close()


def plot_trajectory(all_triads_combined, output_dir, radius_um, condition_areas, conditions,
                    report_radius_um=None, trajectory_min_um=0):
    """
    Distance accumulation trajectory plot.

    Parameters
    ----------
    radius_um           : search radius (trajectory x-axis max)
    report_radius_um    : primary reporting threshold — drawn as a vertical line
    trajectory_min_um   : where the trajectory x-axis starts (default 0)
    """
    if all_triads_combined is None or all_triads_combined.empty:
        return

    colors     = {c: col for c, col in zip(conditions, ["steelblue", "tomato", "seagreen", "darkorange"])}
    linestyles = {c: ls  for c, ls  in zip(conditions, ["-", "--", "-.", ":"])}
    thresholds = np.linspace(trajectory_min_um, radius_um, 300)

    atc = all_triads_combined.copy()
    # "Triad at distance d" = both partners within d µm of the anchor DC.
    # We use max(anchor→p1, anchor→p2), not the p1↔p2 distance.
    atc["max_anchor_dist_um"] = atc[["dist_anchor_p1_um", "dist_anchor_p2_um"]].max(axis=1)

    fig, axes = plt.subplots(2, 1, figsize=(10, 12))

    for cond in conditions:
        sub = atc[atc["condition"] == cond]
        if sub.empty:
            continue
        area_mm2 = condition_areas.get(cond, np.nan) / 1e6
        raw_counts = [int((sub["max_anchor_dist_um"] <= t).sum()) for t in thresholds]
        densities  = [c / area_mm2 if area_mm2 > 0 else 0 for c in raw_counts]
        axes[0].plot(thresholds, raw_counts,  color=colors.get(cond, "gray"), linestyle=linestyles.get(cond, "-"), linewidth=2.5, label=f"{cond}  (n={len(sub)})")
        axes[1].plot(thresholds, densities,   color=colors.get(cond, "gray"), linestyle=linestyles.get(cond, "-"), linewidth=2.5, label=f"{cond}  (area={area_mm2:.4f} mm²)")

    # Vertical lines: search radius + reporting threshold
    for ax in axes:
        ax.axvline(radius_um, color="gray", linestyle=":", linewidth=1.2, label=f"Search radius ({radius_um} µm)")
        if report_radius_um and report_radius_um != radius_um:
            ax.axvline(report_radius_um, color="black", linestyle="--", linewidth=1.5, label=f"Report threshold ({report_radius_um} µm)")
        ax.set_xlabel("Distance threshold (µm)", fontsize=12)
        ax.set_xlim(trajectory_min_um, radius_um)
        ax.legend(fontsize=10)
        ax.grid(alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("Cumulative number of triads", fontsize=12)
    axes[0].set_title(f"Triad Accumulation by Distance\n({'  vs  '.join(conditions)})", fontsize=11)
    axes[1].set_ylabel("Cumulative triads per mm²", fontsize=12)
    axes[1].set_title("Triad Density Accumulation by Distance", fontsize=11)

    plt.tight_layout(pad=3)
    fig.savefig(os.path.join(output_dir, "trajectory_triads_by_distance.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_condition_comparison(summary_df, all_triads_combined, output_dir, radius_um, condition_areas, conditions,
                              report_radius_um=None, trajectory_min_um=0):
    colors = {c: col for c, col in zip(conditions, ["steelblue", "tomato", "seagreen", "darkorange"])}
    summary_df = summary_df[summary_df["condition"].isin(conditions)].copy()
    if summary_df.empty:
        return

    grp = summary_df.groupby(["condition", "combo_label"])["n_triads"].sum().reset_index()
    grp["area_mm2"]       = grp["condition"].map(lambda c: condition_areas.get(c, np.nan) / 1e6)
    grp["triads_per_mm2"] = grp["n_triads"] / grp["area_mm2"]

    combos = sorted(grp["combo_label"].unique())
    x, width = np.arange(len(combos)), 0.8 / max(len(conditions), 1)

    fig1, axes1 = plt.subplots(2, 1, figsize=(max(8, len(combos) * 3), 12))
    for i, cond in enumerate(conditions):
        sub_grp = grp[grp["condition"] == cond]
        count_vals   = [sub_grp.loc[sub_grp["combo_label"] == c, "n_triads"].sum()       if c in sub_grp["combo_label"].values else 0 for c in combos]
        density_vals = [sub_grp.loc[sub_grp["combo_label"] == c, "triads_per_mm2"].sum() if c in sub_grp["combo_label"].values else 0 for c in combos]
        for ax, vals in zip(axes1, [count_vals, density_vals]):
            bars = ax.bar(x + i * width, vals, width, label=cond, color=colors.get(cond, "gray"), alpha=0.85, edgecolor="white")
            for bar in bars:
                h = bar.get_height()
                if h > 0:
                    ax.text(bar.get_x() + bar.get_width() / 2, h * 1.02,
                            f"{int(h):,}" if ax == axes1[0] else f"{h:.2f}",
                            ha="center", va="bottom", fontsize=10, fontweight="bold")

    for ax, ylabel, title in zip(
        axes1,
        ["Number of Triads", "Triads per mm²"],
        [f"Triad Counts (radius = {radius_um} µm)", "Triad Density"],
    ):
        ax.set_xticks(x + width * (len(conditions) - 1) / 2)
        ax.set_xticklabels(combos, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout(pad=3)
    fig1.savefig(os.path.join(output_dir, "condition_comparison_counts_density.png"), dpi=150, bbox_inches="tight")
    plt.close(fig1)

    if all_triads_combined is None or all_triads_combined.empty:
        return

    atc = all_triads_combined[all_triads_combined["condition"].isin(conditions)].copy()
    if atc.empty:
        return

    dist_pairs = {
        "Anchor – Partner1": "dist_anchor_p1_um",
        "Anchor – Partner2": "dist_anchor_p2_um",
        "Partner1 – Partner2": "dist_p1_p2_um",
    }
    pair_labels, dist_cols = list(dist_pairs.keys()), list(dist_pairs.values())
    xd = np.arange(len(pair_labels))

    fig2, ax2 = plt.subplots(figsize=(10, 7))
    for i, cond in enumerate(conditions):
        sub   = atc[atc["condition"] == cond]
        means = [sub[col].mean() if not sub[col].empty else 0 for col in dist_cols]
        sems  = [sub[col].sem()  if not sub[col].empty else 0 for col in dist_cols]
        bars = ax2.bar(xd + i * width, means, width, label=f"{cond} (n={len(sub):,})",
                       color=colors.get(cond, "gray"), alpha=0.85, yerr=sems, capsize=5,
                       edgecolor="white", error_kw=dict(elinewidth=1.5, ecolor="black"))
        for bar, m, s in zip(bars, means, sems):
            if m > 0:
                ax2.text(bar.get_x() + bar.get_width() / 2, m + s + 0.5,
                         f"{m:.1f} µm", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax2.set_xticks(xd + width * (len(conditions) - 1) / 2)
    ax2.set_xticklabels(pair_labels, fontsize=11)
    ax2.set_ylabel("Mean Distance (µm)", fontsize=12)
    ax2.set_title(f"Mean Pairwise Distances in Triads (radius = {radius_um} µm, error bars = SEM)", fontsize=11)
    ax2.legend(fontsize=10)
    ax2.grid(axis="y", alpha=0.3)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    plt.tight_layout()
    fig2.savefig(os.path.join(output_dir, "condition_comparison_distances_um.png"), dpi=150, bbox_inches="tight")
    plt.close(fig2)

    grp.to_csv(os.path.join(output_dir, "condition_comparison_counts.csv"), index=False)
    dist_summary = atc.groupby("condition")[dist_cols].agg(["mean", "std", "sem", "count"]).round(3)
    dist_summary.columns = ["__".join(c) for c in dist_summary.columns]
    dist_summary.to_csv(os.path.join(output_dir, "condition_comparison_distances.csv"))

    plot_trajectory(atc, output_dir, radius_um, condition_areas, conditions,
                    report_radius_um=report_radius_um,
                    trajectory_min_um=trajectory_min_um)


# ── Functional marker analysis ────────────────────────────────────────────────

def _run_functional_marker_analysis(
    all_cells_df: pd.DataFrame,
    output_dir: str,
    func_cfg: dict,
    report_radius_um: float,
    anchor_type: str,
    partner1_type: str,
    partner2_type: str,
    anchor_name: str,
    partner1_name: str,
    partner2_name: str,
    conditions: list,
) -> None:
    """
    Compare functional marker expression between in-triad vs not-in-triad cells.

    For each cell type (anchor, partner1, partner2):
      1. In-triad vs not-in-triad violin plots + Mann-Whitney U
      2. Condition A vs B within in-triad cells

    Parameters
    ----------
    all_cells_df : pd.DataFrame
        Aggregated cell table across all images.  Must contain:
          cell_type, condition, image_id,
          _in_triad_anchor, _in_triad_partner1, _in_triad_partner2,
          and all marker intensity columns listed in func_cfg["markers"].
    """
    from scipy.stats import mannwhitneyu

    markers = func_cfg.get("markers", {})
    if not markers:
        print("[functional] No markers configured — skipping.")
        return

    func_out = os.path.join(output_dir, "functional_analysis")
    os.makedirs(func_out, exist_ok=True)

    COLOR_IN  = "#ff7f0e"   # orange  — in-triad
    COLOR_OUT = "#1f77b4"   # blue    — out-of-triad
    COND_COLORS = dict(zip(conditions, ["#2196F3", "#F44336", "#4CAF50", "#FF9800"]))

    cell_role_map = [
        (anchor_type,   anchor_name,   "_in_triad_anchor"),
        (partner1_type, partner1_name, "_in_triad_partner1"),
        (partner2_type, partner2_name, "_in_triad_partner2"),
    ]

    all_summary = []

    for cell_type, display_name, flag_col in cell_role_map:
        if flag_col not in all_cells_df.columns:
            print(f"[functional] ⚠️  Flag column '{flag_col}' missing — skipping {display_name}")
            continue

        cells = all_cells_df[all_cells_df["cell_type"] == cell_type].copy()
        if cells.empty:
            print(f"[functional] ⚠️  No {display_name} cells in dataset.")
            continue

        # Only keep markers whose columns actually exist in the data
        avail = {name: col for name, col in markers.items() if col in cells.columns}
        missing = [name for name in markers if name not in avail]
        if missing:
            print(f"[functional]   Columns not found for {display_name}: {missing}")
        if not avail:
            print(f"[functional] ⚠️  No marker columns found for {display_name} — skipping.")
            continue

        in_triad     = cells[cells[flag_col] == True]
        not_in_triad = cells[cells[flag_col] == False]
        n_in  = len(in_triad)
        n_out = len(not_in_triad)
        print(f"\n[functional] {display_name}: {n_in:,} in-triad  |  {n_out:,} not-in-triad")

        n_m = len(avail)

        # ── Plot 1: In-triad vs Not-in-triad ─────────────────────────────────
        fig1, axes1 = plt.subplots(1, n_m, figsize=(max(4 * n_m, 8), 6))
        axes1 = [axes1] if n_m == 1 else list(axes1)

        for ax, (marker_name, col) in zip(axes1, avail.items()):
            vals_in  = in_triad[col].dropna().values.astype(float)
            vals_out = not_in_triad[col].dropna().values.astype(float)

            if len(vals_in) >= 3 and len(vals_out) >= 3:
                stat, pval = mannwhitneyu(vals_in, vals_out, alternative="two-sided")
            else:
                stat, pval = np.nan, np.nan

            if len(vals_out) > 0 and len(vals_in) > 0:
                parts = ax.violinplot([vals_out, vals_in], positions=[0, 1],
                                      showmedians=True, showextrema=False)
                for pc, color in zip(parts["bodies"], [COLOR_OUT, COLOR_IN]):
                    pc.set_facecolor(color)
                    pc.set_alpha(0.75)
                parts["cmedians"].set_color("black")
                parts["cmedians"].set_linewidth(2)

            sig = ("***" if pval < 0.001 else "**" if pval < 0.01
                   else "*" if pval < 0.05 else "ns")
            ax.set_xticks([0, 1])
            ax.set_xticklabels(
                [f"Out of triad\n(n={len(vals_out):,})", f"In triad\n(n={len(vals_in):,})"],
                fontsize=8,
            )
            pval_str = f"{pval:.2e}" if not np.isnan(pval) else "n/a"
            ax.set_title(f"{marker_name}\n{sig}  p={pval_str}", fontsize=9, fontweight="bold")
            ax.set_ylabel("Expression (intensity)", fontsize=9)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            all_summary.append({
                "cell_type":       display_name,
                "comparison":      "in_triad_vs_out",
                "marker":          marker_name,
                "marker_col":      col,
                "group_A":         "in_triad",
                "group_B":         "not_in_triad",
                "n_A":             len(vals_in),
                "n_B":             len(vals_out),
                "mean_A":          float(np.mean(vals_in))   if len(vals_in)  else np.nan,
                "mean_B":          float(np.mean(vals_out))  if len(vals_out) else np.nan,
                "median_A":        float(np.median(vals_in)) if len(vals_in)  else np.nan,
                "median_B":        float(np.median(vals_out))if len(vals_out) else np.nan,
                "mannwhitney_U":   stat,
                "p_value":         pval,
                "significance":    sig,
            })

        fig1.suptitle(
            f"{display_name}: In-triad vs Out-of-triad Marker Expression\n"
            f"(triad threshold = {report_radius_um} µm)",
            fontsize=12, fontweight="bold",
        )
        plt.tight_layout()
        fig1.savefig(
            os.path.join(func_out, f"functional_{display_name}_intriad_vs_out.png"),
            dpi=150, bbox_inches="tight",
        )
        plt.close(fig1)

        # ── Plot 2: Condition comparison (CLR vs DII) within in-triad ────────
        if len(conditions) >= 2 and "condition" in cells.columns and n_in >= 2:
            cond_sub = {c: in_triad[in_triad["condition"] == c] for c in conditions}
            valid_conds = [c for c in conditions if len(cond_sub[c]) >= 3]

            if len(valid_conds) >= 2:
                fig2, axes2 = plt.subplots(1, n_m, figsize=(max(4 * n_m, 8), 6))
                axes2 = [axes2] if n_m == 1 else list(axes2)

                for ax, (marker_name, col) in zip(axes2, avail.items()):
                    vals_per_cond = [
                        cond_sub[c][col].dropna().values.astype(float)
                        for c in valid_conds
                    ]

                    if len(vals_per_cond[0]) >= 3 and len(vals_per_cond[1]) >= 3:
                        stat2, pval2 = mannwhitneyu(vals_per_cond[0], vals_per_cond[1],
                                                    alternative="two-sided")
                    else:
                        stat2, pval2 = np.nan, np.nan

                    if all(len(v) > 0 for v in vals_per_cond):
                        parts2 = ax.violinplot(vals_per_cond,
                                               positions=range(len(valid_conds)),
                                               showmedians=True, showextrema=False)
                        for pc, cond in zip(parts2["bodies"], valid_conds):
                            pc.set_facecolor(COND_COLORS.get(cond, "#888888"))
                            pc.set_alpha(0.75)
                        parts2["cmedians"].set_color("black")
                        parts2["cmedians"].set_linewidth(2)

                    sig2 = ("***" if pval2 < 0.001 else "**" if pval2 < 0.01
                            else "*" if pval2 < 0.05 else "ns")
                    ax.set_xticks(range(len(valid_conds)))
                    ax.set_xticklabels(
                        [f"{c}\n(n={len(cond_sub[c]):,})" for c in valid_conds],
                        fontsize=8,
                    )
                    pval2_str = f"{pval2:.2e}" if not np.isnan(pval2) else "n/a"
                    ax.set_title(f"{marker_name}\n{sig2}  p={pval2_str}",
                                 fontsize=9, fontweight="bold")
                    ax.set_ylabel("Expression (intensity)", fontsize=9)
                    ax.spines["top"].set_visible(False)
                    ax.spines["right"].set_visible(False)

                    all_summary.append({
                        "cell_type":     display_name,
                        "comparison":    f"{'_vs_'.join(valid_conds)}_in_triad",
                        "marker":        marker_name,
                        "marker_col":    col,
                        "group_A":       valid_conds[0],
                        "group_B":       valid_conds[1],
                        "n_A":           len(vals_per_cond[0]),
                        "n_B":           len(vals_per_cond[1]),
                        "mean_A":        float(np.mean(vals_per_cond[0])) if len(vals_per_cond[0]) else np.nan,
                        "mean_B":        float(np.mean(vals_per_cond[1])) if len(vals_per_cond[1]) else np.nan,
                        "median_A":      float(np.median(vals_per_cond[0])) if len(vals_per_cond[0]) else np.nan,
                        "median_B":      float(np.median(vals_per_cond[1])) if len(vals_per_cond[1]) else np.nan,
                        "mannwhitney_U": stat2,
                        "p_value":       pval2,
                        "significance":  sig2,
                    })

                fig2.suptitle(
                    f"{display_name}: {' vs '.join(valid_conds)} — In-triad cells only\n"
                    f"(triad threshold = {report_radius_um} µm)",
                    fontsize=12, fontweight="bold",
                )
                plt.tight_layout()
                fig2.savefig(
                    os.path.join(func_out,
                                 f"functional_{display_name}_condition_compare.png"),
                    dpi=150, bbox_inches="tight",
                )
                plt.close(fig2)

    if all_summary:
        summary_df = pd.DataFrame(all_summary)
        n_tests = len(summary_df)
        bonferroni_alpha = 0.05 / max(n_tests, 1)
        summary_df["bonferroni_alpha"] = round(bonferroni_alpha, 8)
        summary_df["sig_bonferroni"]   = summary_df["p_value"] < bonferroni_alpha
        summary_df.to_csv(os.path.join(func_out, "functional_marker_summary.csv"), index=False)
        sig_hits = summary_df["sig_bonferroni"].sum()
        print(f"\n[functional] ✓  {n_tests} tests, Bonferroni α={bonferroni_alpha:.3e}, "
              f"{sig_hits} hits after correction")
        print(f"[functional]    Outputs → {func_out}")
    else:
        print("[functional] No results to save.")


# ── Main entry point ──────────────────────────────────────────────────────────

def run_triad_analysis(cfg: dict) -> None:
    """
    Run triad analysis for an entire experiment defined by cfg.

    Parameters
    ----------
    cfg : dict
        Loaded from config_example.yaml (via yaml.safe_load).
    """
    # ── Unpack config ─────────────────────────────────────────
    exp_name      = cfg["experiment"]["name"]
    conditions    = cfg["experiment"]["conditions"]
    img_cond_map  = cfg["experiment"].get("image_condition_map", {})
    mpp           = cfg["imaging"]["microns_per_pixel"]
    cond_areas    = _get_condition_areas(cfg)

    input_dir     = cfg["paths"]["input_dir"]
    output_dir    = cfg["paths"]["output_dir"]
    qc_plot_dir   = os.path.join(output_dir, "qc_plots_triads")

    t_cfg          = cfg["analysis"]["triad"]
    radius_um      = t_cfg["radius_um"]
    radius_px      = radius_um / mpp
    anchor_type    = t_cfg.get("anchor_type")
    partner1_type  = t_cfg.get("partner_type_1")
    partner2_type  = t_cfg.get("partner_type_2")
    matched_only   = t_cfg.get("matched_only", True)
    min_triad_size = t_cfg.get("min_triad_size", 1)

    # Short display names — used in labels, plot titles, output filenames.
    # Fallback to the full cell type string if not set.
    anchor_name   = t_cfg.get("anchor_name")   or (anchor_type   or "anchor")
    partner1_name = t_cfg.get("partner_1_name") or (partner1_type or "partner1")
    partner2_name = t_cfg.get("partner_2_name") or (partner2_type or "partner2")

    # Distance trajectory range
    report_radius_um   = t_cfg.get("report_radius_um", radius_um)
    trajectory_min_um  = t_cfg.get("trajectory_min_um", 0)

    # Functional analysis is now a separate pipeline step (spatia/analysis/functional.py).
    # triads.py no longer runs it inline — just set func_enabled=False to avoid the
    # legacy embedded path.
    func_enabled = False

    os.makedirs(output_dir,   exist_ok=True)
    os.makedirs(qc_plot_dir,  exist_ok=True)

    print(f"[SPATIA] Experiment : {exp_name}")
    print(f"[SPATIA] Conditions : {conditions}")
    print(f"[SPATIA] Radius     : {radius_um} µm  →  {radius_px:.2f} px")
    print(f"[SPATIA] Input      : {input_dir}")
    print(f"[SPATIA] Output     : {output_dir}\n")

    # ── Discover input files ──────────────────────────────────
    csv_files = sorted([
        f for f in os.listdir(input_dir)
        if f.endswith("_matched_with_boundaries.csv") and not f.startswith("._")
    ])
    if not csv_files:
        print(f"[SPATIA] ⚠️  No matched CSV files found in {input_dir}")
        return
    print(f"[SPATIA] Found {len(csv_files)} file(s).\n")

    all_summary_rows = []
    all_triad_dfs    = []
    all_func_cell_dfs = []   # accumulated for functional analysis (only if enabled)

    for csv_file in csv_files:
        image_id  = csv_file.replace("_matched_with_boundaries.csv", "")
        condition = _get_condition(image_id, img_cond_map, conditions)
        print(f"{'='*65}")
        print(f"Processing : {image_id}  [condition: {condition}]")

        # ── Load CSV (try multiple encodings) ────────────────
        df = None
        for enc in ["utf-8", "latin-1", "cp1252", "utf-8-sig"]:
            try:
                df = pd.read_csv(os.path.join(input_dir, csv_file), encoding=enc)
                break
            except (UnicodeDecodeError, Exception):
                continue
        if df is None:
            print(f"  ⚠️  Could not read {csv_file} — skipping.")
            continue

        print(f"  Total cells : {len(df)}")
        if matched_only and "matched" in df.columns:
            df = df[df["matched"] == True].reset_index(drop=True)
            print(f"  After matched filter : {len(df)}")
        if df.empty:
            print("  ⚠️  No cells — skipping.")
            continue
        if "cell_id" not in df.columns:
            df.insert(0, "cell_id", df.index.astype(str))

        tree = cKDTree(df[["centroid_x", "centroid_y"]].values)
        cell_types = sorted(df["cell_type"].dropna().unique().tolist())
        print(f"  Cell types ({len(cell_types)}): {cell_types}")

        # ── Build combo list ──────────────────────────────────
        if anchor_type and partner1_type and partner2_type:
            combos = [(anchor_type, partner1_type, partner2_type)]
        else:
            combos = sorted(set(
                perm
                for trio in combinations(cell_types, 3)
                for perm in permutations(trio)
            ))

        image_triad_dfs = []
        image_summary   = []

        for (at, p1t, p2t) in combos:
            if at not in cell_types or p1t not in cell_types or p2t not in cell_types:
                continue
            triad_df = find_triads(df, at, p1t, p2t, radius_px, tree, mpp)
            n = len(triad_df)
            if n < min_triad_size:
                continue

            # Use short names when anchor/partners match the configured types
            def _short(full, configured, short):
                return short if full == configured else full
            al = _short(at,  anchor_type,   anchor_name)
            p1l = _short(p1t, partner1_type, partner1_name)
            p2l = _short(p2t, partner2_type, partner2_name)
            combo_label = f"{al}__{p1l}__{p2l}"
            print(f"    {combo_label:<60}: {n:>6} triads")

            # Overwrite the full-name label set inside find_triads with short names
            triad_df["triad_combo_label"] = combo_label
            triad_df["image_id"]  = image_id
            triad_df["condition"] = condition
            image_triad_dfs.append(triad_df)
            image_summary.append({
                "image_id":             image_id,
                "condition":            condition,
                "anchor_type":          at,
                "partner1_type":        p1t,
                "partner2_type":        p2t,
                "combo_label":          combo_label,
                "n_triads":             n,
                "radius_um":            radius_um,
                "radius_px":            round(radius_px, 3),
                "mean_dist_anchor_p1_um":  round(triad_df["dist_anchor_p1_um"].mean(), 3),
                "mean_dist_anchor_p2_um":  round(triad_df["dist_anchor_p2_um"].mean(), 3),
                "mean_dist_p1_p2_um":      round(triad_df["dist_p1_p2_um"].mean(), 3),
            })

        if image_triad_dfs:
            all_triads_df = pd.concat(image_triad_dfs, ignore_index=True)
            all_triad_dfs.append(all_triads_df)
            all_triads_df.to_csv(os.path.join(output_dir, f"{image_id}_triad_pairs.csv"), index=False)

            # Flag cells
            triad_ids    = set(all_triads_df["anchor_cell_id"].astype(str)) | set(all_triads_df["partner1_cell_id"].astype(str)) | set(all_triads_df["partner2_cell_id"].astype(str))
            df["is_triad_anchor"]   = df["cell_id"].astype(str).isin(set(all_triads_df["anchor_cell_id"].astype(str)))
            df["is_triad_partner1"] = df["cell_id"].astype(str).isin(set(all_triads_df["partner1_cell_id"].astype(str)))
            df["is_triad_partner2"] = df["cell_id"].astype(str).isin(set(all_triads_df["partner2_cell_id"].astype(str)))
            df["in_any_triad"]      = df["cell_id"].astype(str).isin(triad_ids)
            df.to_csv(os.path.join(output_dir, f"{image_id}_cells_with_triad_flags.csv"), index=False)

            # Accumulate cell-level data for functional analysis (in-triad flagged
            # at func_radius so the comparison uses the tighter reporting threshold)
            if func_enabled:
                report_triads = all_triads_df[
                    all_triads_df[["dist_anchor_p1_um", "dist_anchor_p2_um"]].max(axis=1)
                    <= func_radius
                ] if not all_triads_df.empty else all_triads_df

                df_func = df.copy()
                df_func["image_id"]  = image_id
                df_func["condition"] = condition
                df_func["_in_triad_anchor"]   = df_func["cell_id"].astype(str).isin(
                    set(report_triads["anchor_cell_id"].astype(str))
                )
                df_func["_in_triad_partner1"] = df_func["cell_id"].astype(str).isin(
                    set(report_triads["partner1_cell_id"].astype(str))
                )
                df_func["_in_triad_partner2"] = df_func["cell_id"].astype(str).isin(
                    set(report_triads["partner2_cell_id"].astype(str))
                )
                all_func_cell_dfs.append(df_func)

            # QC plot for top combo — right panel filtered to report_radius_um (e.g. 15 µm)
            top_combo    = max(image_summary, key=lambda x: x["n_triads"])
            top_triad_df = all_triads_df[all_triads_df["triad_combo_label"] == top_combo["combo_label"]]
            plot_triad_qc(
                df, top_triad_df, image_id, top_combo["combo_label"],
                os.path.join(qc_plot_dir, f"{image_id}_{top_combo['combo_label']}_triad_qc.png"),
                radius_um,
                report_radius_um=report_radius_um,
            )
        else:
            print(f"  ⚠️  No triads found at radius={radius_um} µm")

        all_summary_rows.extend(image_summary)
        print()

    # ── Global outputs ────────────────────────────────────────
    if not all_summary_rows:
        print("[SPATIA] ⚠️  No triads found in any image.")
        return

    summary_df = pd.DataFrame(all_summary_rows)
    summary_df.to_csv(os.path.join(output_dir, "triad_summary.csv"), index=False)

    all_triads_combined = pd.concat(all_triad_dfs, ignore_index=True) if all_triad_dfs else pd.DataFrame()

    top20 = summary_df.groupby("combo_label")["n_triads"].sum().sort_values(ascending=False).head(20).reset_index()
    print("=" * 65)
    print(f"[SPATIA] TRIAD ANALYSIS COMPLETE — {exp_name}")
    print("=" * 65)
    print(top20.to_string(index=False))

    plot_condition_comparison(summary_df, all_triads_combined, output_dir, radius_um, cond_areas, conditions,
                              report_radius_um=report_radius_um,
                              trajectory_min_um=trajectory_min_um)

    # ── Functional marker analysis ─────────────────────────────
    if func_enabled and all_func_cell_dfs:
        print(f"\n{'='*65}")
        print("[SPATIA] FUNCTIONAL MARKER ANALYSIS")
        print(f"{'='*65}")
        all_cells_df = pd.concat(all_func_cell_dfs, ignore_index=True)
        print(f"[functional] Aggregated {len(all_cells_df):,} cells from "
              f"{len(all_func_cell_dfs)} images")
        _run_functional_marker_analysis(
            all_cells_df     = all_cells_df,
            output_dir       = output_dir,
            func_cfg         = func_cfg,
            report_radius_um = func_radius,
            anchor_type      = anchor_type,
            partner1_type    = partner1_type,
            partner2_type    = partner2_type,
            anchor_name      = anchor_name,
            partner1_name    = partner1_name,
            partner2_name    = partner2_name,
            conditions       = conditions,
        )
    elif func_enabled:
        print("[functional] No triad images — skipping functional analysis.")

    print(f"\n[SPATIA] All outputs saved to: {output_dir}")
