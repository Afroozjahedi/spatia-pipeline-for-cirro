"""
spatia.analysis.survival
========================
Survival analysis step: correlates per-patient triad density with OS and DFS.

Reads
-----
  {output_dir}/triad_summary.csv    — per-image triad counts (from triads step)
  cfg["analysis"]["survival"]["patient_annotation_file"]  — patient-level annotation
      containing OS, OS_Censor, DFS, DFS_Censor, Patient, Group, TMA spot / region

Workflow
--------
  1. Build a complete image-level table (all 140 images, including those with 0 triads).
  2. Aggregate to patient level: total triads / total tissue area → triads per mm².
  3. Merge with survival data on patient ID (via TMA spot/region → reg number mapping).
  4. Split patients into HIGH vs LOW triad density (median cutoff, configurable).
  5. Kaplan-Meier curves for OS and DFS; log-rank test.
  6. Optional: multivariate Cox PH with condition (CLR/DII) as covariate.

Outputs land in {output_dir}/survival/:
  km_OS.png
  km_DFS.png
  patient_triad_density.csv
  survival_logrank_results.csv

Config shape (crc_tma.yaml)
----------------------------
    analysis:
      survival:
        enabled: true
        patient_annotation_file: "/path/to/patient_with_tls_class.csv"
        area_per_image_um2: 1168128.0    # 1920 × 1440 × 0.65² µm²
        split_by: "median"               # or a numeric threshold in triads/mm²
        covariates: ["condition"]        # for Cox PH (optional)
        # Convention: OS_Censor=0 → event (death), 1 → censored (alive)
        event_col_os:  "OS_Censor"
        event_col_dfs: "DFS_Censor"
        censor_is_one: true              # true: Censor=1 means censored; false: inverted
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


# ── Helpers ───────────────────────────────────────────────────────────────────

def _spot_to_reg(spot_str: str) -> list[int]:
    """'1,2' → [1, 2]"""
    try:
        return [int(s.strip()) for s in str(spot_str).split(",") if s.strip()]
    except Exception:
        return []


def _reg_to_image_prefix(reg_num: int) -> str:
    """1 → 'reg001'"""
    return f"reg{reg_num:03d}"


def _build_patient_density(
    triad_summary: pd.DataFrame,
    surv_df: pd.DataFrame,
    all_image_ids: list[str],
    area_per_image_um2: float,
) -> pd.DataFrame:
    """
    Aggregate triad counts to patient level.

    surv_df must have columns: Patient, Group, 'TMA spot / region'
    all_image_ids: every image that was processed (including those with 0 triads).
    """
    # Build a per-image triad count (default 0 for images not in summary)
    img_counts = (
        triad_summary.groupby("image_id")["n_triads"].sum()
        .reindex(all_image_ids, fill_value=0)
        .reset_index()
        .rename(columns={"index": "image_id", 0: "n_triads"})
    )
    img_counts.columns = ["image_id", "n_triads"]

    # Extract reg_num from image_id (e.g. "CLR_reg001_A" → 1)
    img_counts["reg_num"] = (
        img_counts["image_id"]
        .str.extract(r"reg(\d+)", expand=False)
        .astype(float)
        .astype("Int64")
    )

    rows = []
    for _, patient_row in surv_df.dropna(subset=["Patient"]).iterrows():
        patient_id = int(patient_row["Patient"])
        spots      = _spot_to_reg(patient_row["TMA spot / region"])
        group      = patient_row.get("Group", np.nan)

        # Collect all images belonging to this patient's TMA spots
        patient_imgs = img_counts[img_counts["reg_num"].isin(spots)]
        n_images     = len(patient_imgs)
        total_triads = int(patient_imgs["n_triads"].sum())
        total_area_mm2 = n_images * area_per_image_um2 / 1e6

        rows.append({
            "patient_id":          patient_id,
            "group":               int(group) if not pd.isna(group) else np.nan,
            "n_images":            n_images,
            "total_triads":        total_triads,
            "total_area_mm2":      round(total_area_mm2, 4),
            "triads_per_mm2":      round(total_triads / total_area_mm2, 4)
                                   if total_area_mm2 > 0 else 0.0,
        })

    density_df = pd.DataFrame(rows)

    # Merge OS / DFS from surv_df
    surv_cols = [c for c in ["Patient", "Group", "OS", "OS_Censor", "DFS", "DFS_Censor",
                              "MSI_IHC", "pT", "pN", "pM", "Sex", "Age"]
                 if c in surv_df.columns]
    surv_slim = surv_df[surv_cols].copy()
    surv_slim = surv_slim.rename(columns={"Patient": "patient_id"})
    surv_slim["patient_id"] = surv_slim["patient_id"].astype("Int64")
    density_df["patient_id"] = density_df["patient_id"].astype("Int64")

    merged = density_df.merge(surv_slim, on="patient_id", how="left")
    return merged


def _km_plot(
    df: pd.DataFrame,
    duration_col: str,
    event_col: str,
    censor_is_one: bool,
    group_col: str,
    title: str,
    save_path: str,
    split_label: str = "median",
) -> dict:
    """
    Kaplan-Meier curves for two groups defined by group_col (high/low density).
    Returns dict with log-rank p-value and stats.
    """
    df = df.dropna(subset=[duration_col, event_col, group_col]).copy()

    # event_observed: 1 = event happened (death/recurrence)
    if censor_is_one:
        # Censor=1 means censored → event = 1 - Censor
        df["_event"] = 1 - df[event_col].astype(int)
    else:
        # Censor=0 means censored → event = Censor
        df["_event"] = df[event_col].astype(int)

    groups = sorted(df[group_col].unique())
    if len(groups) < 2:
        print(f"  ⚠️  Only one group in {group_col} — cannot plot KM.")
        return {}

    COLORS = {"high": "#d62728", "low": "#1f77b4", "High": "#d62728", "Low": "#1f77b4"}

    fig, ax = plt.subplots(figsize=(9, 6))
    results = {}

    for grp in groups:
        sub = df[df[group_col] == grp]
        kmf = KaplanMeierFitter()
        kmf.fit(sub[duration_col], event_observed=sub["_event"],
                label=f"{grp} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=True,
                                   color=COLORS.get(grp, None), linewidth=2)
        results[grp] = {"n": len(sub), "events": int(sub["_event"].sum())}

    if len(groups) == 2:
        g1 = df[df[group_col] == groups[0]]
        g2 = df[df[group_col] == groups[1]]
        lr = logrank_test(
            g1[duration_col], g2[duration_col],
            event_observed_A=g1["_event"], event_observed_B=g2["_event"],
        )
        pval = lr.p_value
        sig  = ("***" if pval < 0.001 else "**" if pval < 0.01
                else "*" if pval < 0.05 else "ns")
        ax.text(0.02, 0.08, f"Log-rank p = {pval:.3e}  {sig}",
                transform=ax.transAxes, fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        results["logrank_p"] = pval
        results["significance"] = sig

    ax.set_xlabel("Time (months)", fontsize=12)
    ax.set_ylabel("Survival probability", fontsize=12)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return results


# ── Main entry point ──────────────────────────────────────────────────────────

def run_survival_analysis(cfg: dict) -> None:
    """
    Survival analysis step — called by run_pipeline.py.
    Depends on the triads step having written triad_summary.csv.
    """
    exp_name   = cfg["experiment"]["name"]
    output_dir = cfg["paths"]["output_dir"]

    s_cfg = cfg.get("analysis", {}).get("survival", {})
    if not s_cfg.get("enabled", False):
        print("[survival] disabled in config — skipping.")
        return

    annot_file        = s_cfg.get("patient_annotation_file", "")
    area_per_image    = s_cfg.get("area_per_image_um2", 1_168_128.0)
    split_by          = s_cfg.get("split_by", "median")
    censor_is_one     = s_cfg.get("censor_is_one", True)

    if not annot_file or not os.path.exists(annot_file):
        print(f"[survival] ⚠️  Patient annotation file not found: {annot_file!r}")
        print("[survival]    Set analysis.survival.patient_annotation_file in YAML.")
        return

    summary_csv = os.path.join(output_dir, "triad_summary.csv")
    if not os.path.exists(summary_csv):
        print(f"[survival] ⚠️  triad_summary.csv not found in {output_dir}")
        print("[survival]    Run the triads step first.")
        return

    surv_out = os.path.join(output_dir, "survival")
    os.makedirs(surv_out, exist_ok=True)

    print(f"[survival] Experiment       : {exp_name}")
    print(f"[survival] Annotation file  : {annot_file}")
    print(f"[survival] Area/image       : {area_per_image:,.0f} µm²")
    print(f"[survival] Split by         : {split_by}\n")

    # ── Load inputs ───────────────────────────────────────────────────────────
    triad_summary = pd.read_csv(summary_csv)
    surv_df       = pd.read_csv(annot_file)

    # Collect all image IDs that were processed (including 0-triad images)
    input_dir  = cfg["paths"]["input_dir"]
    conditions = cfg["experiment"]["conditions"]
    img_cond_map = cfg["experiment"].get("image_condition_map", {})

    all_image_ids = []
    if os.path.exists(input_dir):
        for f in os.listdir(input_dir):
            if f.endswith("_matched_with_boundaries.csv") and not f.startswith("._"):
                all_image_ids.append(f.replace("_matched_with_boundaries.csv", ""))
    if not all_image_ids:
        # Fall back to images present in the output directory
        all_image_ids = list(triad_summary["image_id"].unique())
        print("[survival]   ⚠️  Could not list input_dir — using only images with triads.")

    # ── Build patient-level density ───────────────────────────────────────────
    patient_df = _build_patient_density(triad_summary, surv_df, all_image_ids, area_per_image)
    print(f"[survival] Patients with survival data : {patient_df['OS'].notna().sum()}")
    print(f"[survival] Patients with triads        : "
          f"{(patient_df['total_triads'] > 0).sum()}")

    # ── Split high / low triad density ───────────────────────────────────────
    if split_by == "median":
        threshold = patient_df["triads_per_mm2"].median()
        print(f"[survival] Median triad density : {threshold:.4f} triads/mm²")
        # When median == 0 (most patients have zero triads), a median split puts
        # everyone in "High" (≥ 0). Fall back to presence vs absence.
        if threshold == 0:
            print("[survival]   Median is 0 — splitting by triad presence (any vs zero).")
            patient_df["density_group"] = np.where(
                patient_df["triads_per_mm2"] > 0, "Triad+", "Triad−"
            )
        else:
            patient_df["density_group"] = np.where(
                patient_df["triads_per_mm2"] >= threshold, "High", "Low"
            )
    else:
        try:
            threshold = float(split_by)
        except ValueError:
            print(f"[survival] ⚠️  Invalid split_by: {split_by!r}. Using median.")
            threshold = patient_df["triads_per_mm2"].median()
        patient_df["density_group"] = np.where(
            patient_df["triads_per_mm2"] >= threshold, "High", "Low"
        )

    # Save patient-level table
    patient_csv = os.path.join(surv_out, "patient_triad_density.csv")
    patient_df.to_csv(patient_csv, index=False)
    print(f"[survival] Patient table saved → {patient_csv}")
    print(patient_df[["patient_id", "total_triads", "triads_per_mm2",
                       "density_group"]].to_string(index=False))

    # ── KM plots ──────────────────────────────────────────────────────────────
    lr_results = []

    for outcome, dur_col, cens_col in [
        ("Overall Survival (OS)", "OS", "OS_Censor"),
        ("Disease-Free Survival (DFS)", "DFS", "DFS_Censor"),
    ]:
        if dur_col not in patient_df.columns:
            print(f"[survival] ⚠️  Column {dur_col!r} not found — skipping {outcome}.")
            continue

        save_path = os.path.join(surv_out, f"km_{dur_col}.png")
        title = (f"{exp_name} — {outcome}\n"
                 f"High (≥{threshold:.3f}) vs Low triad density  |  "
                 f"split by: {split_by}")

        print(f"\n[survival] Plotting {outcome} …")
        res = _km_plot(
            patient_df,
            duration_col   = dur_col,
            event_col      = cens_col,
            censor_is_one  = censor_is_one,
            group_col      = "density_group",
            title          = title,
            save_path      = save_path,
            split_label    = split_by,
        )
        if "logrank_p" in res:
            print(f"  Log-rank p = {res['logrank_p']:.4e}  {res['significance']}")
            lr_results.append({
                "outcome":       outcome,
                "split_by":      split_by,
                "threshold":     threshold,
                "logrank_p":     res["logrank_p"],
                "significance":  res["significance"],
                **{f"n_{k}": v["n"] for k, v in res.items() if isinstance(v, dict)},
            })

    if lr_results:
        lr_df = pd.DataFrame(lr_results)
        lr_df.to_csv(os.path.join(surv_out, "survival_logrank_results.csv"), index=False)

    # ── Also plot KM by condition (CLR vs DII) ────────────────────────────────
    if "Group" in patient_df.columns or "group" in patient_df.columns:
        grp_col = "Group" if "Group" in patient_df.columns else "group"
        group_map = {1: "CLR", 2: "DII"}
        # Group may be stored as float (1.0, 2.0) — convert to int before mapping
        patient_df["condition_label"] = (
            patient_df[grp_col]
            .apply(lambda x: group_map.get(int(x)) if pd.notna(x) else None)
            .fillna("Unknown")
        )

        for outcome, dur_col, cens_col in [
            ("Overall Survival (OS)", "OS", "OS_Censor"),
            ("Disease-Free Survival (DFS)", "DFS", "DFS_Censor"),
        ]:
            if dur_col not in patient_df.columns:
                continue
            save_path = os.path.join(surv_out, f"km_condition_{dur_col}.png")
            print(f"\n[survival] Plotting {outcome} by condition (CLR vs DII) …")
            res = _km_plot(
                patient_df,
                duration_col  = dur_col,
                event_col     = cens_col,
                censor_is_one = censor_is_one,
                group_col     = "condition_label",
                title         = f"{exp_name} — {outcome} by condition (CLR vs DII)",
                save_path     = save_path,
            )
            if "logrank_p" in res:
                print(f"  CLR vs DII log-rank p = {res['logrank_p']:.4e}  {res['significance']}")

    print(f"\n[survival] ✓  All outputs → {surv_out}")
