"""
spatia.analysis.functional
==========================
Standalone functional marker analysis step.

Reads the per-image outputs already written by the triads step:
  {output_dir}/triad_data/{image_id}_cells_with_triad_flags.csv  — cell-level data + triad membership
  {output_dir}/triad_data/{image_id}_triad_pairs.csv             — pairwise distances for each triad
(falls back to {output_dir}/{image_id}_*.csv directly for older runs predating
the triad_data/ subfolder, added 2026-10-09)

Recomputes in-triad membership at the functional reporting radius (which can be
tighter than the search radius used during detection) and then, for each configured
cell type, compares marker expression:

  1. In-triad vs not-in-triad  (violin + Mann-Whitney U, all experiment_groups pooled)
  2. Experiment group A vs experiment group B within in-triad cells only  (CLR vs DII, etc.)

Outputs land in {output_dir}/functional_analysis/:
  functional_{CellType}_intriad_vs_out.png
  functional_{CellType}_experiment_group_compare.png
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

def _triad_file(output_dir: str, filename: str) -> str:
    """Resolves a triads-step output file: prefers output_dir/triad_data/<filename>
    (where triads.py has written it since 2026-10-09), falling back to
    output_dir/<filename> directly for runs produced before that change."""
    new_path = os.path.join(output_dir, "triad_data", filename)
    if os.path.exists(new_path):
        return new_path
    return os.path.join(output_dir, filename)


def _get_experiment_group(image_id: str, image_experiment_group_map: dict, experiment_groups: list) -> str:
    if image_id in image_experiment_group_map:
        return image_experiment_group_map[image_id]
    for cond in experiment_groups:
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
    experiment_groups   = cfg["experiment"]["groups"]
    img_group_map = cfg["experiment"].get("image_experiment_group_map", {})
    output_dir   = cfg["paths"]["output_dir"]

    func_cfg = cfg.get("analysis", {}).get("functional", {})
    if not func_cfg.get("enabled", False):
        print("[functional] disabled in config — skipping.")
        return

    # `markers` is the shared/default marker set. Optional per-role overrides
    # (markers_anchor/markers_partner1/markers_partner2) let you check a
    # different marker panel per cell-type role — e.g. exhaustion markers on
    # a T-cell partner but maturation markers on a DC anchor. A role with no
    # override falls back to the shared `markers` dict unchanged.
    markers           = func_cfg.get("markers", {})
    markers_anchor    = func_cfg.get("markers_anchor",   markers)
    markers_partner1  = func_cfg.get("markers_partner1", markers)
    markers_partner2  = func_cfg.get("markers_partner2", markers)
    func_radius    = func_cfg.get("report_radius_um", 20.0)
    t_cfg          = cfg.get("analysis", {}).get("triad", {})
    anchor_type    = t_cfg.get("anchor_type",    "anchor")
    partner1_type  = t_cfg.get("partner_type_1", "partner1")
    partner2_type  = t_cfg.get("partner_type_2", "partner2")
    anchor_name    = t_cfg.get("anchor_name",    anchor_type)
    partner1_name  = t_cfg.get("partner_1_name", partner1_type)
    partner2_name  = t_cfg.get("partner_2_name", partner2_type)

    if not (t_cfg.get("anchor_type") and t_cfg.get("partner_type_1") and t_cfg.get("partner_type_2")):
        print("[functional] ⚠️  analysis.triad.anchor_type/partner_type_1/partner_type_2 are not "
              "all set — cell_type will never match the placeholder values "
              f"('{anchor_type}'/'{partner1_type}'/'{partner2_type}'), so this run will find 0 "
              "cells for every role and produce no output. Set all three in the config.")

    if not (markers_anchor or markers_partner1 or markers_partner2):
        print("[functional] No markers configured — add analysis.functional.markers (shared) "
              "or markers_anchor/markers_partner1/markers_partner2 (per-role) to YAML.")
        return

    func_out = os.path.join(output_dir, "functional_analysis")
    os.makedirs(func_out, exist_ok=True)

    print(f"[functional] Experiment  : {exp_name}")
    print(f"[functional] Experiment groups  : {experiment_groups}")
    print(f"[functional] Func radius : {func_radius} µm")
    print(f"[functional] Markers (shared) : {list(markers.keys())}")
    if markers_anchor is not markers:
        print(f"[functional]   anchor override   : {list(markers_anchor.keys())}")
    if markers_partner1 is not markers:
        print(f"[functional]   partner1 override : {list(markers_partner1.keys())}")
    if markers_partner2 is not markers:
        print(f"[functional]   partner2 override : {list(markers_partner2.keys())}")
    print(f"[functional] Input dir   : {output_dir}\n")

    # ── Discover per-image cells files ────────────────────────────────────────
    cell_files = sorted(glob.glob(
        os.path.join(output_dir, "triad_data", "*_cells_with_triad_flags.csv")
    )) or sorted(glob.glob(
        os.path.join(output_dir, "*_cells_with_triad_flags.csv")
    ))
    if not cell_files:
        print(f"[functional] ⚠️  No *_cells_with_triad_flags.csv found in {output_dir} or {output_dir}/triad_data")
        print("[functional]    Run the triads step first.")
        return
    print(f"[functional] Found {len(cell_files)} cell files.\n")

    # ── Load and accumulate per-image data ───────────────────────────────────
    all_cell_dfs = []

    for cell_path in cell_files:
        basename = os.path.basename(cell_path)
        # e.g. CLR_reg001_A_cells_with_triad_flags.csv → CLR_reg001_A
        image_id = basename.replace("_cells_with_triad_flags.csv", "")
        experiment_group = _get_experiment_group(image_id, img_group_map, experiment_groups)

        # Load cell-level data (contains all original columns if KEEP_INTENSITIES=True)
        try:
            df = pd.read_csv(cell_path)
        except Exception as e:
            print(f"  ⚠️  Cannot read {basename}: {e} — skipping.")
            continue

        if "cell_id" not in df.columns:
            df.insert(0, "cell_id", df.index.astype(str))
        # Same subtype grouping the triads step used (no-op if already applied
        # or no analysis.triad.type_groups configured).
        from spatia.analysis.triads import apply_triad_type_groups
        df = apply_triad_type_groups(df, t_cfg)

        # Load corresponding triad pairs (may not exist if no triads were found)
        pairs_path = _triad_file(output_dir, f"{image_id}_triad_pairs.csv")
        triad_pairs = pd.DataFrame()
        if os.path.exists(pairs_path):
            try:
                triad_pairs = pd.read_csv(pairs_path)
            except Exception as e:
                print(f"  ⚠️  Could not read {os.path.basename(pairs_path)}: {e} — "
                      f"treating {image_id} as having no in-triad cells (not 'no triads found').")

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
        df["experiment_group"] = experiment_group
        all_cell_dfs.append(df)

    if not all_cell_dfs:
        print("[functional] ⚠️  No cell data loaded — check that the triads step ran first.")
        return

    all_cells = pd.concat(all_cell_dfs, ignore_index=True)
    print(f"[functional] Aggregated {len(all_cells):,} cells from {len(all_cell_dfs)} images\n")

    # ── Per-cell-type analysis ────────────────────────────────────────────────
    # 2026-11-05 rewrite of the statistics:
    #  * The unit of replication is the IMAGE (core / ROI), not the cell.
    #    Each image contributes one median per marker. Old plots tested
    #    cells (e.g. "CLR n=42, DII n=15" were 57 in-triad CD4 *cells*),
    #    which treats neighbouring cells of one core as independent samples.
    #  * In-triad vs out-of-triad is a PAIRED test within each image that has
    #    in-triad cells (Wilcoxon signed-rank), so core-to-core staining
    #    differences cancel out.
    #  * Out-of-triad cells from cores with NO triads are now loaded too
    #    (the triads step writes their flags file) and shown as their own box.
    #  * CLR vs DII compares per-image medians of in-triad cells (Mann-Whitney).
    #  * BH-FDR across all per-image tests; per-cell p-values kept in the CSV.
    from scipy.stats import wilcoxon
    from spatia.analysis.plot_style import group_color, add_sig_bracket, bh_fdr, p_to_stars
    COLOR_OUT_NOTRIAD = "#cfcfcf"
    COLOR_OUT = "#8c8c8c"
    COLOR_IN = "#6A3D9A"
    min_cells_img = int(func_cfg.get("min_cells_per_image", 3))

    imgs_with_triads = set(all_cells.loc[
        all_cells[["_in_triad_anchor", "_in_triad_partner1", "_in_triad_partner2"]].any(axis=1),
        "image_id"])

    cell_role_map = [
        (anchor_type,   anchor_name,   "_in_triad_anchor",   markers_anchor),
        (partner1_type, partner1_name, "_in_triad_partner1", markers_partner1),
        (partner2_type, partner2_name, "_in_triad_partner2", markers_partner2),
    ]

    def _img_medians(sub, col):
        g = sub.groupby("image_id")[col]
        med = g.median()
        n = g.size()
        return med[n >= min_cells_img].dropna()

    def _box(ax, data, positions, colors, labels):
        bp = ax.boxplot(data, positions=positions, widths=0.55, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black", linewidth=2))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.8)
        rng = np.random.default_rng(0)
        for pos, d in zip(positions, data):
            ax.scatter(pos + rng.uniform(-0.12, 0.12, len(d)), d, s=14, color="black",
                       alpha=0.6, zorder=3, linewidths=0)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    all_summary = []

    for cell_type, display_name, flag_col, role_markers in cell_role_map:
        cells = all_cells[all_cells["cell_type"] == cell_type].copy()
        if cells.empty:
            print(f"[functional] ⚠️  No {display_name} cells — skipping.")
            continue

        avail   = {n: c for n, c in role_markers.items() if c in cells.columns}
        missing = [n for n in role_markers if n not in avail]
        if missing:
            print(f"[functional]   Marker columns not found for {display_name}: {missing}")
        if not avail:
            print(f"[functional] ⚠️  No usable marker columns for {display_name} — skipping.")
            continue

        in_triad     = cells[cells[flag_col] == True]
        out_tri_img  = cells[(cells[flag_col] == False) & cells["image_id"].isin(imgs_with_triads)]
        out_no_img   = cells[(cells[flag_col] == False) & ~cells["image_id"].isin(imgs_with_triads)]
        n_in = len(in_triad)
        print(f"[functional] {display_name}: {n_in:,} in-triad cells in {in_triad['image_id'].nunique()} images  |  "
              f"{len(out_tri_img):,} out-of-triad (same images)  |  "
              f"{len(out_no_img):,} in {out_no_img['image_id'].nunique()} images without triads")

        n_m = len(avail)

        # ── Plot 1: in vs out of triad, per image ───────────────────────────
        fig1, axes1 = plt.subplots(1, n_m, figsize=(max(3.6 * n_m, 8), 5.8))
        axes1 = [axes1] if n_m == 1 else list(axes1)
        for ax, (marker_name, col) in zip(axes1, avail.items()):
            m_in = _img_medians(in_triad, col)
            m_out = _img_medians(out_tri_img, col)
            m_none = _img_medians(out_no_img, col)
            paired = pd.concat([m_out.rename("out"), m_in.rename("in")], axis=1).dropna()
            for _, r in paired.iterrows():
                ax.plot([1, 2], [r["out"], r["in"]], color="#bbbbbb", linewidth=0.7, zorder=1)
            _box(ax, [m_none.values, m_out.values, m_in.values], [0, 1, 2],
                 [COLOR_OUT_NOTRIAD, COLOR_OUT, COLOR_IN],
                 [f"Out\n(cores w/o\ntriads)\nn={len(m_none)}", f"Out\n(same\ncores)\nn={len(m_out)}",
                  f"In triad\nn={len(m_in)}"])
            p_pair = (wilcoxon(paired["in"], paired["out"]).pvalue
                      if len(paired) >= 5 and (paired["in"] != paired["out"]).any() else np.nan)
            p_all = (mannwhitneyu(m_in, pd.concat([m_out, m_none]), alternative="two-sided").pvalue
                     if len(m_in) >= 3 and (len(m_out) + len(m_none)) >= 3 else np.nan)
            top = np.nanmax([np.nanmax(d) if len(d) else np.nan for d in (m_none, m_out, m_in)] + [0])
            if np.isfinite(p_pair):
                add_sig_bracket(ax, 1, 2, top * 1.04 + 1e-9, p_pair, label_fmt="{stars} paired p={p:.2g}",
                                fontsize=7)
            ax.set_title(marker_name, fontsize=10, fontweight="bold")
            ax.set_ylabel("Median expression per image", fontsize=9)
            # per-cell test kept for the CSV only
            vi, vo = in_triad[col].dropna().astype(float), out_tri_img[col].dropna().astype(float)
            p_cell = (mannwhitneyu(vi, vo, alternative="two-sided").pvalue
                      if len(vi) >= 3 and len(vo) >= 3 else np.nan)
            all_summary.append({
                "cell_type": display_name, "comparison": "in_triad_vs_out", "marker": marker_name,
                "marker_col": col, "group_A": "in_triad", "group_B": "not_in_triad",
                "unit": "image (paired)", "n_images_A": len(m_in), "n_images_B": len(m_out),
                "n_images_paired": len(paired), "n_images_no_triads": len(m_none),
                "n_cells_A": len(vi), "n_cells_B": len(vo),
                "median_A": float(m_in.median()) if len(m_in) else np.nan,
                "median_B": float(m_out.median()) if len(m_out) else np.nan,
                "median_no_triad_cores": float(m_none.median()) if len(m_none) else np.nan,
                "p_value": p_pair, "test": "Wilcoxon signed-rank on per-image medians",
                "p_in_vs_all_out_cores": p_all, "p_per_cell_mannwhitney": p_cell,
            })
        fig1.suptitle(f"{display_name}: in-triad vs out-of-triad — one dot per image (core/ROI)\n"
                      f"(triad radius {func_radius} µm; images need ≥{min_cells_img} cells; "
                      f"lines join the same image)", fontsize=11, fontweight="bold")
        plt.tight_layout()
        fig1.savefig(os.path.join(func_out, f"functional_{display_name}_intriad_vs_out.png"),
                     dpi=150, bbox_inches="tight")
        plt.close(fig1)

        # ── Plot 2: experiment groups, in-triad cells, per image ─────────────
        if len(experiment_groups) >= 2 and n_in >= 6:
            gA, gB = experiment_groups[0], experiment_groups[1]
            fig2, axes2 = plt.subplots(1, n_m, figsize=(max(3.2 * n_m, 8), 5.6))
            axes2 = [axes2] if n_m == 1 else list(axes2)
            for ax, (marker_name, col) in zip(axes2, avail.items()):
                mA = _img_medians(in_triad[in_triad["experiment_group"] == gA], col)
                mB = _img_medians(in_triad[in_triad["experiment_group"] == gB], col)
                nA_c = int((in_triad["experiment_group"] == gA).sum())
                nB_c = int((in_triad["experiment_group"] == gB).sum())
                _box(ax, [mA.values, mB.values], [0, 1], [group_color(gA), group_color(gB)],
                     [f"{gA}\n{len(mA)} images\n({nA_c} cells)", f"{gB}\n{len(mB)} images\n({nB_c} cells)"])
                p2 = (mannwhitneyu(mA, mB, alternative="two-sided").pvalue
                      if len(mA) >= 3 and len(mB) >= 3 else np.nan)
                top = np.nanmax([np.nanmax(d) if len(d) else np.nan for d in (mA, mB)] + [0])
                if np.isfinite(p2):
                    add_sig_bracket(ax, 0, 1, top * 1.04 + 1e-9, p2, fontsize=8)
                else:
                    ax.text(0.5, 0.97, "too few images to test", transform=ax.transAxes,
                            ha="center", va="top", fontsize=7, color="#777777")
                ax.set_title(marker_name, fontsize=10, fontweight="bold")
                ax.set_ylabel("Median expression per image", fontsize=9)
                vA = in_triad.loc[in_triad["experiment_group"] == gA, col].dropna().astype(float)
                vB = in_triad.loc[in_triad["experiment_group"] == gB, col].dropna().astype(float)
                p_cell = (mannwhitneyu(vA, vB, alternative="two-sided").pvalue
                          if len(vA) >= 3 and len(vB) >= 3 else np.nan)
                all_summary.append({
                    "cell_type": display_name, "comparison": f"{gA}_vs_{gB}_in_triad",
                    "marker": marker_name, "marker_col": col, "group_A": gA, "group_B": gB,
                    "unit": "image", "n_images_A": len(mA), "n_images_B": len(mB),
                    "n_cells_A": len(vA), "n_cells_B": len(vB),
                    "median_A": float(mA.median()) if len(mA) else np.nan,
                    "median_B": float(mB.median()) if len(mB) else np.nan,
                    "p_value": p2, "test": "Mann-Whitney on per-image medians",
                    "p_per_cell_mannwhitney": p_cell,
                })
            fig2.suptitle(f"{display_name}: {gA} vs {gB} — in-triad cells, one dot per image\n"
                          f"(triad radius {func_radius} µm)", fontsize=11, fontweight="bold")
            plt.tight_layout()
            fig2.savefig(os.path.join(func_out, f"functional_{display_name}_experiment_group_compare.png"),
                         dpi=150, bbox_inches="tight")
            plt.close(fig2)

    # ── Summary CSV ───────────────────────────────────────────────────────────
    if all_summary:
        summary_df = pd.DataFrame(all_summary)
        n_tests = int(summary_df["p_value"].notna().sum())
        summary_df["q_bh"] = bh_fdr(summary_df["p_value"].values)
        summary_df["sig_bh_0.05"] = summary_df["q_bh"] < 0.05
        bonf_alpha = 0.05 / max(n_tests, 1)
        summary_df["bonferroni_alpha"] = round(bonf_alpha, 8)
        summary_df["sig_bonferroni"] = summary_df["p_value"] < bonf_alpha
        summary_df["significance"] = [p_to_stars(q) for q in summary_df["q_bh"]]
        out_csv = os.path.join(func_out, "functional_marker_summary.csv")
        summary_df.to_csv(out_csv, index=False)
        print(f"\n[functional] ✓  {n_tests} per-image tests  |  BH q<0.05: "
              f"{int(summary_df['sig_bh_0.05'].sum())}  |  Bonferroni: {int(summary_df['sig_bonferroni'].sum())}")
        hits = summary_df[summary_df["sig_bh_0.05"]]
        for _, r in hits.iterrows():
            print(f"[functional]    hit: {r['cell_type']} {r['marker']} {r['comparison']} "
                  f"median {r['median_A']:.3g} vs {r['median_B']:.3g}  q={r['q_bh']:.2g}")
        print(f"[functional]    Outputs → {func_out}")
    else:
        print("[functional] ⚠️  No results produced.")
