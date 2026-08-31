"""
spatia.analysis.triads
======================
Generalizable triad detection for spatial proteomics data.

All parameters come from a config dict (loaded from config_example.yaml).

Usage
-----
    from spatia.analysis.triads import run_triad_analysis
    run_triad_analysis(config)   # config is a dict loaded from YAML

The config dict shape mirrors config_example.yaml.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from itertools import combinations, permutations


# ── Config helpers ────────────────────────────────────────────────────────────

def _get_experiment_group_areas(cfg: dict) -> dict:
    """
    Returns a flat {experiment_group: area_um2} dict.
    Handles both flat and timepoint-nested config formats.
    """
    raw = cfg["imaging"]["experiment_group_areas_um2"]
    # If values are dicts, it's timepoint-nested — flatten by summing areas.
    if raw and isinstance(next(iter(raw.values())), dict):
        flat = {}
        for tp_areas in raw.values():
            for cond, area in tp_areas.items():
                flat[cond] = flat.get(cond, 0) + area
        return flat
    return dict(raw)


def _get_experiment_group(image_id: str, image_experiment_group_map: dict, experiment_groups: list) -> str:
    """
    Determine experiment_group from image_id.
    Priority: explicit map → prefix match → 'Unknown'.
    """
    if image_id in image_experiment_group_map:
        return image_experiment_group_map[image_id]
    for cond in experiment_groups:
        if image_id.upper().startswith(cond.upper()):
            return cond
    return "Unknown"


# ── Per-image tissue area from QuPath ROI extraction (optional) ────────────────

def _load_roi_areas(roi_labels_dir: str, mpp: float, annotation_class: str = None) -> dict:
    """
    Reads every *_roi_labels.txt produced by the QuPath ROI-extraction script
    (00_ROI_extract_mask_project.groovy) in roi_labels_dir, sums the Area_px2
    column per image (optionally filtered to a single annotation Class), and
    converts to um^2 using the pipeline's own microns_per_pixel — NOT QuPath's
    internal pixel calibration, which may be unset or unverified. This keeps
    area and triad-radius-in-pixels conversion using the same single source
    of truth for calibration.

    Returns {qupath_image_name: area_um2}. Returns {} if roi_labels_dir is
    unset/missing (callers should treat that as "feature not in use").
    """
    areas = {}
    if not roi_labels_dir or not os.path.isdir(roi_labels_dir):
        return areas

    label_files = sorted([
        f for f in os.listdir(roi_labels_dir)
        if f.endswith("_roi_labels.txt") and not f.startswith("._")
    ])
    for fname in label_files:
        image_name = fname.replace("_roi_labels.txt", "")
        path = os.path.join(roi_labels_dir, fname)
        try:
            tbl = pd.read_csv(path, sep="\t")
        except Exception as e:
            print(f"[SPATIA] ⚠️  Could not read ROI labels file {fname}: {e}")
            continue

        if "Area_px2" not in tbl.columns:
            print(f"[SPATIA] ⚠️  {fname} has no Area_px2 column — re-run the "
                  f"updated QuPath extraction script. Skipping this image.")
            continue

        sub = tbl
        if annotation_class is not None:
            if "Class" not in tbl.columns:
                print(f"[SPATIA] ⚠️  {fname} has no Class column — cannot filter "
                      f"to '{annotation_class}'. Skipping this image.")
                continue
            sub = tbl[tbl["Class"] == annotation_class]
            if sub.empty:
                print(f"[SPATIA] ⚠️  {fname}: no annotations with Class == "
                      f"'{annotation_class}' — area for this image will be 0.")

        areas[image_name] = float(sub["Area_px2"].sum()) * (mpp ** 2)

    return areas


def _normalize_name(s: str) -> str:
    """Lowercase, strip everything but letters/digits — for fuzzy name matching."""
    return re.sub(r"[^a-z0-9]", "", s.lower())


def _resolve_image_areas(image_ids: list, roi_areas_by_qupath_name: dict,
                          image_id_map: dict = None):
    """
    Matches each image_id (from *_matched_with_boundaries.csv) to a
    QuPath-derived area, in priority order:
      1. Explicit imaging.image_id_map override ({qupath_name: image_id}).
      2. Exact string match (image_id == qupath image name).
      3. Normalized match (lowercase, alphanumeric-only).

    Returns (matched: {image_id: area_um2}, unmatched: [image_id, ...]).
    """
    image_id_map = image_id_map or {}
    override_lookup = {v: k for k, v in image_id_map.items()}  # image_id -> qupath_name

    norm_lookup = {}
    for qname in roi_areas_by_qupath_name:
        norm_lookup.setdefault(_normalize_name(qname), qname)

    matched, unmatched = {}, []
    for image_id in image_ids:
        qname = None
        if image_id in override_lookup and override_lookup[image_id] in roi_areas_by_qupath_name:
            qname = override_lookup[image_id]
        elif image_id in roi_areas_by_qupath_name:
            qname = image_id
        else:
            qname = norm_lookup.get(_normalize_name(image_id))

        if qname is not None:
            matched[image_id] = roi_areas_by_qupath_name[qname]
        else:
            unmatched.append(image_id)

    return matched, unmatched


def _resolve_per_image_group_areas(
    image_ids: list,
    image_id_to_group: dict,
    roi_labels_dir: str,
    mpp: float,
    area_annotation_class,
    image_id_map: dict,
    fallback_group_areas: dict,
) -> dict:
    """
    Builds per-experiment_group total area from QuPath-derived per-image ROI
    areas, with a 3-tier per-image fallback:
      1. Matched QuPath ROI-labels area (measured).
      2. No match, but fallback_group_areas has a constant for that image's
         group -> impute area = group_constant / n_images_in_group_this_run
         (spreads the existing group total evenly across the images actually
         present in this run, rather than the images assumed when the
         constant was originally typed into the config).
      3. Neither -> NaN; that image's area is excluded from the group total.
         Its triads are NOT excluded from the count elsewhere, so that
         group's density will be a slight overestimate — flagged loudly
         below rather than silently.

    Returns {experiment_group: area_um2}. Prints a per-image audit line so
    the provenance of every density number is visible in the run log.
    """
    roi_areas_by_name = _load_roi_areas(roi_labels_dir, mpp, area_annotation_class)
    matched, unmatched = _resolve_image_areas(image_ids, roi_areas_by_name, image_id_map)

    if unmatched:
        print(f"[SPATIA] ⚠️  {len(unmatched)} image(s) had no matching QuPath ROI "
              f"area (tried exact name, normalized name, and imaging.image_id_map):")
        for u in unmatched:
            print(f"           - {u}")
        print(f"           Add an entry to imaging.image_id_map ({{qupath_name: image_id}}) "
              f"to fix a specific one; unresolved images fall back to a group average.")

    # Count images per group (all images this run, matched + unmatched) — the
    # denominator used to spread a group constant evenly when imputing.
    group_counts = {}
    for iid in image_ids:
        g = image_id_to_group[iid]
        group_counts[g] = group_counts.get(g, 0) + 1

    per_image_area = {}
    for iid in image_ids:
        g = image_id_to_group[iid]
        if iid in matched:
            per_image_area[iid] = matched[iid]
            print(f"[SPATIA] Area  {iid:<30s}: {matched[iid]:>14,.1f} µm²  (measured, QuPath ROI)")
            continue

        constant = fallback_group_areas.get(g)
        has_constant = constant is not None and not (isinstance(constant, float) and np.isnan(constant))
        if has_constant and group_counts.get(g, 0) > 0:
            imputed = constant / group_counts[g]
            per_image_area[iid] = imputed
            print(f"[SPATIA] Area  {iid:<30s}: {imputed:>14,.1f} µm²  "
                  f"(no QuPath match — imputed as group '{g}' average, "
                  f"{constant:,.1f} / {group_counts[g]} images)")
        else:
            per_image_area[iid] = float("nan")
            print(f"[SPATIA] Area  {iid:<30s}: {'NaN':>14s}  ⚠️  no QuPath match and no "
                  f"imaging.experiment_group_areas_um2 for group '{g}' — excluded from "
                  f"'{g}' density (triad counts for this image are NOT excluded, so "
                  f"'{g}' density will be a slight overestimate).")

    resolved_group_areas = {}
    for iid, area in per_image_area.items():
        if not np.isnan(area):
            g = image_id_to_group[iid]
            resolved_group_areas[g] = resolved_group_areas.get(g, 0.0) + area

    return resolved_group_areas


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


def plot_trajectory(all_triads_combined, output_dir, radius_um, experiment_group_areas, experiment_groups,
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

    colors     = {c: col for c, col in zip(experiment_groups, ["steelblue", "tomato", "seagreen", "darkorange"])}
    linestyles = {c: ls  for c, ls  in zip(experiment_groups, ["-", "--", "-.", ":"])}
    thresholds = np.linspace(trajectory_min_um, radius_um, 300)

    atc = all_triads_combined.copy()
    # "Triad at distance d" = both partners within d µm of the anchor DC.
    # We use max(anchor→p1, anchor→p2), not the p1↔p2 distance.
    atc["max_anchor_dist_um"] = atc[["dist_anchor_p1_um", "dist_anchor_p2_um"]].max(axis=1)

    fig, axes = plt.subplots(2, 1, figsize=(10, 12))

    for cond in experiment_groups:
        sub = atc[atc["experiment_group"] == cond]
        if sub.empty:
            continue
        area_mm2 = experiment_group_areas.get(cond, np.nan) / 1e6
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
    axes[0].set_title(f"Triad Accumulation by Distance\n({'  vs  '.join(experiment_groups)})", fontsize=11)
    axes[1].set_ylabel("Cumulative triads per mm²", fontsize=12)
    axes[1].set_title("Triad Density Accumulation by Distance", fontsize=11)

    plt.tight_layout(pad=3)
    fig.savefig(os.path.join(output_dir, "trajectory_triads_by_distance.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_experiment_group_comparison(summary_df, all_triads_combined, output_dir, radius_um, experiment_group_areas, experiment_groups,
                              report_radius_um=None, trajectory_min_um=0):
    colors = {c: col for c, col in zip(experiment_groups, ["steelblue", "tomato", "seagreen", "darkorange"])}
    summary_df = summary_df[summary_df["experiment_group"].isin(experiment_groups)].copy()
    if summary_df.empty:
        return

    grp = summary_df.groupby(["experiment_group", "combo_label"])["n_triads"].sum().reset_index()
    grp["area_mm2"]       = grp["experiment_group"].map(lambda c: experiment_group_areas.get(c, np.nan) / 1e6)
    grp["triads_per_mm2"] = grp["n_triads"] / grp["area_mm2"]

    combos = sorted(grp["combo_label"].unique())
    x, width = np.arange(len(combos)), 0.8 / max(len(experiment_groups), 1)

    fig1, axes1 = plt.subplots(2, 1, figsize=(max(8, len(combos) * 3), 12))
    for i, cond in enumerate(experiment_groups):
        sub_grp = grp[grp["experiment_group"] == cond]
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
        ax.set_xticks(x + width * (len(experiment_groups) - 1) / 2)
        ax.set_xticklabels(combos, rotation=30, ha="right", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.legend(fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout(pad=3)
    fig1.savefig(os.path.join(output_dir, "experiment_group_comparison_counts_density.png"), dpi=150, bbox_inches="tight")
    plt.close(fig1)

    if all_triads_combined is None or all_triads_combined.empty:
        return

    atc = all_triads_combined[all_triads_combined["experiment_group"].isin(experiment_groups)].copy()
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
    for i, cond in enumerate(experiment_groups):
        sub   = atc[atc["experiment_group"] == cond]
        means = [sub[col].mean() if not sub[col].empty else 0 for col in dist_cols]
        sems  = [sub[col].sem()  if not sub[col].empty else 0 for col in dist_cols]
        bars = ax2.bar(xd + i * width, means, width, label=f"{cond} (n={len(sub):,})",
                       color=colors.get(cond, "gray"), alpha=0.85, yerr=sems, capsize=5,
                       edgecolor="white", error_kw=dict(elinewidth=1.5, ecolor="black"))
        for bar, m, s in zip(bars, means, sems):
            if m > 0:
                ax2.text(bar.get_x() + bar.get_width() / 2, m + s + 0.5,
                         f"{m:.1f} µm", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax2.set_xticks(xd + width * (len(experiment_groups) - 1) / 2)
    ax2.set_xticklabels(pair_labels, fontsize=11)
    ax2.set_ylabel("Mean Distance (µm)", fontsize=12)
    ax2.set_title(f"Mean Pairwise Distances in Triads (radius = {radius_um} µm, error bars = SEM)", fontsize=11)
    ax2.legend(fontsize=10)
    ax2.grid(axis="y", alpha=0.3)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    plt.tight_layout()
    fig2.savefig(os.path.join(output_dir, "experiment_group_comparison_distances_um.png"), dpi=150, bbox_inches="tight")
    plt.close(fig2)

    grp.to_csv(os.path.join(output_dir, "experiment_group_comparison_counts.csv"), index=False)
    dist_summary = atc.groupby("experiment_group")[dist_cols].agg(["mean", "std", "sem", "count"]).round(3)
    dist_summary.columns = ["__".join(c) for c in dist_summary.columns]
    dist_summary.to_csv(os.path.join(output_dir, "experiment_group_comparison_distances.csv"))

    plot_trajectory(atc, output_dir, radius_um, experiment_group_areas, experiment_groups,
                    report_radius_um=report_radius_um,
                    trajectory_min_um=trajectory_min_um)


# ── Per-image result caching ────────────────────────────────────────────────

def _load_cached_triad_results(image_id: str, output_dir: str, radius_um: float, radius_px: float):
    """
    Reconstructs the per-combo summary rows run_triad_analysis() would have
    produced for one image, from its already-written {image_id}_triad_pairs.csv
    — without re-running find_triads(). Backs the skip-if-already-processed
    cache in run_triad_analysis() so adding one new image to a cohort doesn't
    force every existing image's triad detection (which can be expensive —
    see the exhaustive-search warning for unconfigured anchor/partner types)
    to be redone.

    Returns (image_summary: list[dict], triads_df: pd.DataFrame) — same shape
    as a freshly-computed image would produce.
    """
    triads_df = pd.read_csv(os.path.join(output_dir, f"{image_id}_triad_pairs.csv"))

    image_summary = []
    group_cols = ["experiment_group", "anchor_type", "partner1_type", "partner2_type", "triad_combo_label"]
    for keys, grp in triads_df.groupby(group_cols):
        experiment_group, at, p1t, p2t, combo_label = keys
        image_summary.append({
            "image_id":               image_id,
            "experiment_group":       experiment_group,
            "anchor_type":            at,
            "partner1_type":          p1t,
            "partner2_type":          p2t,
            "combo_label":            combo_label,
            "n_triads":               len(grp),
            "radius_um":              radius_um,
            "radius_px":              round(radius_px, 3),
            "mean_dist_anchor_p1_um": round(grp["dist_anchor_p1_um"].mean(), 3),
            "mean_dist_anchor_p2_um": round(grp["dist_anchor_p2_um"].mean(), 3),
            "mean_dist_p1_p2_um":     round(grp["dist_p1_p2_um"].mean(), 3),
        })
    return image_summary, triads_df


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
    experiment_groups    = cfg["experiment"]["groups"]
    img_group_map  = cfg["experiment"].get("image_experiment_group_map", {})
    mpp           = cfg["imaging"]["microns_per_pixel"]
    group_areas    = _get_experiment_group_areas(cfg)

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

    # Distance trajectory range
    report_radius_um   = t_cfg.get("report_radius_um", radius_um)
    trajectory_min_um  = t_cfg.get("trajectory_min_um", 0)

    os.makedirs(output_dir,   exist_ok=True)
    os.makedirs(qc_plot_dir,  exist_ok=True)

    print(f"[SPATIA] Experiment : {exp_name}")
    print(f"[SPATIA] Experiment groups : {experiment_groups}")
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

    # Resolve experiment_group once per image up front (also reused inside the
    # main loop below, instead of recomputing per image).
    image_ids = [f.replace("_matched_with_boundaries.csv", "") for f in csv_files]
    image_id_to_group = {
        iid: _get_experiment_group(iid, img_group_map, experiment_groups)
        for iid in image_ids
    }

    # ── Optional: per-image tissue area from QuPath ROI extraction ─────────
    # If imaging.roi_labels_dir is set, area is measured per image from the
    # QuPath ROI-extraction script's output (with a per-image fallback to a
    # group-average estimate, see _resolve_per_image_group_areas). If unset,
    # behavior is unchanged from before — group_areas stays the static
    # imaging.experiment_group_areas_um2 constant read above.
    roi_labels_dir         = cfg["imaging"].get("roi_labels_dir")
    area_annotation_class  = cfg["imaging"].get("area_annotation_class")
    image_id_map           = cfg["imaging"].get("image_id_map", {})
    if roi_labels_dir:
        print(f"[SPATIA] imaging.roi_labels_dir set — resolving per-image tissue area from QuPath ROI extraction ({roi_labels_dir})")
        group_areas = _resolve_per_image_group_areas(
            image_ids, image_id_to_group, roi_labels_dir, mpp,
            area_annotation_class, image_id_map, group_areas,
        )
        print()

    all_summary_rows = []
    all_triad_dfs    = []

    for csv_file in csv_files:
        image_id  = csv_file.replace("_matched_with_boundaries.csv", "")
        experiment_group = image_id_to_group[image_id]
        print(f"{'='*65}")
        print(f"Processing : {image_id}  [experiment_group: {experiment_group}]")

        # ── Skip if already processed ──────────────────────────
        # Mirrors the skip-if-exists pattern tif_conversion.py / roi_masking.py /
        # segmentation.py / preprocessing.py already use: if both per-image
        # output files exist, load the cached triad_pairs.csv instead of
        # re-running find_triads(), so adding one new image to a cohort
        # doesn't force every existing image to be redetected. No staleness
        # check against the source CSV — delete both cached files for an
        # image to force it to be reprocessed. Images with zero triads never
        # had these files written in the first place (see the "no triads
        # found" branch below), so they're always reprocessed, not cached.
        pairs_path = os.path.join(output_dir, f"{image_id}_triad_pairs.csv")
        flags_path = os.path.join(output_dir, f"{image_id}_cells_with_triad_flags.csv")
        if os.path.exists(pairs_path) and os.path.exists(flags_path):
            print(f"  ⏭️  Already processed — loading cached results")
            cached_summary, cached_triads_df = _load_cached_triad_results(
                image_id, output_dir, radius_um, radius_px
            )
            if cached_summary:
                all_triad_dfs.append(cached_triads_df)
            all_summary_rows.extend(cached_summary)
            print()
            continue

        # ── Load CSV (try multiple encodings) ────────────────
        df = None
        for enc in ["utf-8", "latin-1", "cp1252", "utf-8-sig"]:
            try:
                df = pd.read_csv(os.path.join(input_dir, csv_file), encoding=enc)
                break
            except UnicodeDecodeError:
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
            print(f"  ⚠️  No anchor_type/partner_type_1/partner_type_2 configured — "
                  f"running exhaustive search over {len(combos)} combo(s) "
                  f"({len(cell_types)} cell types), this may be slow.")

        image_triad_dfs = []
        image_summary   = []

        for (at, p1t, p2t) in combos:
            if at not in cell_types or p1t not in cell_types or p2t not in cell_types:
                continue
            triad_df = find_triads(df, at, p1t, p2t, radius_px, tree, mpp)
            n = len(triad_df)
            if n < min_triad_size:
                continue

            combo_label = f"{at}__{p1t}__{p2t}"
            print(f"    {combo_label:<60}: {n:>6} triads")

            triad_df["triad_combo_label"] = combo_label
            triad_df["image_id"]  = image_id
            triad_df["experiment_group"] = experiment_group
            image_triad_dfs.append(triad_df)
            image_summary.append({
                "image_id":             image_id,
                "experiment_group":            experiment_group,
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

    plot_experiment_group_comparison(summary_df, all_triads_combined, output_dir, radius_um, group_areas, experiment_groups,
                              report_radius_um=report_radius_um,
                              trajectory_min_um=trajectory_min_um)

    print(f"\n[SPATIA] All outputs saved to: {output_dir}")
