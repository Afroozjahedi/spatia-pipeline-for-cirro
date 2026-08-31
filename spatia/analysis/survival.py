"""
spatia.analysis.survival
=========================
Cohort-wide survival analysis: correlates per-patient/subject triad density
— and, optionally, per-patient functional-marker exposure — with clinical
outcomes (Kaplan-Meier + log-rank, plus optional multivariate Cox PH).

Rewritten to be dataset-agnostic (previous version hardcoded CRC-TMA-specific
column names, a TMA-spot-number image->patient mapping, and a CLR/DII group
code map — none of that is here anymore; see docs/spatia_analysis_survival.md
for what changed and why).

Image -> patient/subject mapping
---------------------------------
  - Default: image_id IS the patient/subject ID (1:1) — works out of the box
    for any cohort where each processed image is one subject (e.g. a WSI
    cohort like PirB, one image per animal).
  - Optional: analysis.survival.image_patient_map ({image_id: patient_id})
    for cohorts where multiple images belong to one patient (e.g. a TMA with
    several cores per patient) — same override pattern as
    experiment.image_experiment_group_map elsewhere in this pipeline.

Reads
-----
  {output_dir}/triad_summary.csv                     — from triads.py
  {output_dir}/{image_id}_cells_with_triad_flags.csv  — from triads.py, only
      read if analysis.survival.marker_exposures is configured
  {output_dir}/{image_id}_triad_pairs.csv             — from triads.py, only
      read alongside the above (to recompute in-triad membership at
      report_radius_um the same way functional.py does, so a marker-exposure
      number here means the same thing functional.py's own numbers do)
  cfg["analysis"]["survival"]["patient_annotation_file"] — external
      clinical/outcome data. REQUIRED: outcome data (time-to-event, censor
      status) cannot be derived from imaging alone, this always has to come
      from an external file. Column names are config-driven, not assumed —
      see Config shape below.

Why functional_marker_summary.csv itself isn't used for marker exposure
-------------------------------------------------------------------------
That file (written by functional.py) is already pooled across the whole
cohort per (cell_type × marker × comparison) — it has no per-patient
breakdown to correlate against survival. To get a per-patient marker
exposure value, this module re-aggregates the same per-cell files
functional.py reads, the same way functional.py does (same in-triad
recompute logic at report_radius_um), just grouped by patient instead of
pooled across everyone.

Workflow
--------
  1. Enumerate all processed images (from input_dir, including 0-triad ones).
  2. Map each image to a patient/subject (1:1 default, or image_patient_map).
  3. Aggregate triad density to patient level: total triads / total area.
  4. If analysis.survival.marker_exposures is configured, re-aggregate
     per-cell marker intensities to patient level too (cell-count-weighted
     mean across a patient's images).
  5. Merge with the external annotation file on patient_id (config-driven
     column names — no hardcoded CRC-specific columns or clinical fields;
     every column in the annotation file is kept, not a fixed allowlist).
  6. Split into High/Low by triad density (median or numeric threshold, with
     a presence/absence fallback when the median is 0) and, optionally, by
     each configured marker-exposure column (same median-split logic).
  7. Kaplan-Meier + log-rank for every split, for every configured outcome
     (OS/DFS by default, but the outcome list itself is config-driven — not
     hardcoded to exactly two outcomes named OS/DFS).
  8. Optional multivariate Cox PH (if analysis.survival.covariates is set)
     across triad density / marker-exposure / experiment_group covariates —
     this was advertised in the old docstring but never implemented; it's
     real now.

Outputs land in {output_dir}/survival/:
  patient_cohort_summary.csv   (density + marker covariates + outcomes, one row per patient)
  km_<outcome_key>_density.png
  km_<outcome_key>_<marker_exposure_name>.png   (one per configured marker exposure)
  km_<outcome_key>_experiment_group.png
  survival_logrank_results.csv
  cox_ph_<outcome_key>.csv     (one per outcome, only if covariates configured)

Config shape (generic example — no dataset-specific column names implied)
---------------------------------------------------------------------------
    analysis:
      survival:
        enabled: true
        patient_annotation_file: "/path/to/annotation.csv"
        patient_id_col: "patient_id"        # column in the annotation file
        image_patient_map: {}               # optional; {} => image_id is the patient id
        group_col: "group"                  # optional column for experiment_group KM split
        group_code_map: {1: "WT", 2: "KO"}  # optional; falls back to a positional guess
        area_per_image_um2: 1000000.0
        split_by: "median"                  # or a numeric triads/mm² threshold
        censor_is_one: true
        outcomes:
          - name: "Overall Survival (OS)"
            duration_col: "OS"
            event_col: "OS_Censor"
          - name: "Disease-Free Survival (DFS)"
            duration_col: "DFS"
            event_col: "DFS_Censor"
        marker_exposures:                   # optional, per-patient marker covariates
          - name: "GzmB_CD8_intriad"
            cell_type: "CD8 T cells"
            marker_col: "GzmB - cytotoxicity:Cyc_13_ch_2"
            in_triad_only: true             # default true; false = all cells of this type
        report_radius_um: 20.0              # only used if marker_exposures is set —
                                             # should match analysis.functional.report_radius_um
        covariates: []                      # e.g. ["triads_per_mm2", "GzmB_CD8_intriad",
                                             # "experiment_group_label"] — enables Cox PH per outcome
"""

import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

try:
    from lifelines import CoxPHFitter
    _HAS_COX = True
except ImportError:
    _HAS_COX = False


# ── Image -> patient mapping ────────────────────────────────────────────────

def _map_image_to_patient(image_id: str, image_patient_map: dict) -> str:
    """1:1 by default (image_id is the patient id); image_patient_map
    overrides for cohorts where multiple images belong to one patient."""
    if image_patient_map:
        return str(image_patient_map.get(image_id, image_id))
    return str(image_id)


def _normalize_code(x) -> str:
    """Normalizes a group/annotation code so 1, 1.0, '1', '1.0' all compare
    equal — avoids the int()-forcing the previous version used, which broke
    on non-numeric codes."""
    try:
        f = float(x)
        if f.is_integer():
            return str(int(f))
    except (TypeError, ValueError):
        pass
    return str(x)


# ── Triad density aggregation ───────────────────────────────────────────────

def _build_patient_triad_density(
    triad_summary: pd.DataFrame,
    all_image_ids: list,
    image_patient_map: dict,
    area_per_image_um2: float,
) -> pd.DataFrame:
    """Aggregates per-image triad counts (including 0-triad images) up to
    the patient level via the image->patient mapping, and computes
    triads/mm² using a flat per-image area constant.

    Note: area_per_image_um2 is a single constant applied to every image —
    if your images vary a lot in tissue area, this will bias density for
    patients whose images are smaller/larger than average. triads.py has a
    3-tier per-image measured-area fallback (imaging.roi_labels_dir) that
    this module does not yet reuse — a reasonable future enhancement, not
    built here.
    """
    img_counts = (
        triad_summary.groupby("image_id")["n_triads"].sum()
        .reindex(all_image_ids, fill_value=0)
    )

    by_patient: dict = {}
    for image_id in all_image_ids:
        patient_id = _map_image_to_patient(image_id, image_patient_map)
        n_triads = int(img_counts.get(image_id, 0))
        d = by_patient.setdefault(patient_id, {"patient_id": patient_id, "n_images": 0, "total_triads": 0})
        d["n_images"] += 1
        d["total_triads"] += n_triads

    rows = []
    for patient_id, d in by_patient.items():
        total_area_mm2 = d["n_images"] * area_per_image_um2 / 1e6
        rows.append({
            "patient_id":     patient_id,
            "n_images":       d["n_images"],
            "total_triads":   d["total_triads"],
            "total_area_mm2": round(total_area_mm2, 4),
            "triads_per_mm2": round(d["total_triads"] / total_area_mm2, 4)
                               if total_area_mm2 > 0 else 0.0,
        })
    return pd.DataFrame(rows)


# ── Marker exposure aggregation (optional) ──────────────────────────────────

def _compute_in_triad_ids(output_dir: str, image_id: str, report_radius_um: float) -> set:
    """Cell IDs in this image that participate in >=1 triad within
    report_radius_um — mirrors functional.py's own in-triad recompute logic
    (read *_triad_pairs.csv, filter to report_radius_um, union all three
    role's cell IDs) so a cell counted as 'in-triad' here matches what
    functional.py would call in-triad for the same image."""
    pairs_path = os.path.join(output_dir, f"{image_id}_triad_pairs.csv")
    if not os.path.exists(pairs_path):
        return set()
    try:
        pairs = pd.read_csv(pairs_path)
    except Exception as e:
        print(f"[survival]   ⚠️  Could not read {os.path.basename(pairs_path)}: {e} — "
              f"treating {image_id} as having no in-triad cells for marker exposure.")
        return set()
    if pairs.empty or "dist_anchor_p1_um" not in pairs.columns:
        return set()
    pairs_report = pairs[
        pairs[["dist_anchor_p1_um", "dist_anchor_p2_um"]].max(axis=1) <= report_radius_um
    ]
    ids = set()
    for col in ["anchor_cell_id", "partner1_cell_id", "partner2_cell_id"]:
        if col in pairs_report.columns:
            ids |= set(pairs_report[col].astype(str))
    return ids


def _build_patient_marker_exposure(
    output_dir: str,
    all_image_ids: list,
    image_patient_map: dict,
    marker_exposures_cfg: list,
    report_radius_um: float,
) -> "pd.DataFrame | None":
    """For each configured marker exposure, computes a per-image mean marker
    intensity (optionally restricted to in-triad cells of the configured
    cell_type), then aggregates to patient level as a cell-count-weighted
    mean across all of that patient's images — so a patient with more of
    the relevant cell type across their images contributes proportionally,
    rather than every image being weighted equally regardless of n cells.
    Returns None if marker_exposures_cfg is empty (feature not in use)."""
    if not marker_exposures_cfg:
        return None

    per_image_rows = []
    for image_id in all_image_ids:
        cells_path = os.path.join(output_dir, f"{image_id}_cells_with_triad_flags.csv")
        if not os.path.exists(cells_path):
            continue
        try:
            cells = pd.read_csv(cells_path)
        except Exception as e:
            print(f"[survival]   ⚠️  Could not read {os.path.basename(cells_path)}: {e} — "
                  f"skipping {image_id} for marker exposure.")
            continue
        if "cell_id" not in cells.columns:
            cells = cells.copy()
            cells.insert(0, "cell_id", cells.index.astype(str))
        if "cell_type" not in cells.columns:
            print(f"[survival]   ⚠️  {os.path.basename(cells_path)} has no cell_type column — "
                  f"skipping {image_id} for marker exposure.")
            continue

        patient_id = _map_image_to_patient(image_id, image_patient_map)
        in_triad_ids = None  # computed lazily, only if an exposure needs it

        row = {"image_id": image_id, "patient_id": patient_id}
        for exp in marker_exposures_cfg:
            name          = exp["name"]
            cell_type     = exp["cell_type"]
            marker_col    = exp["marker_col"]
            in_triad_only = exp.get("in_triad_only", True)

            sub = cells[cells["cell_type"] == cell_type]
            if in_triad_only:
                if in_triad_ids is None:
                    in_triad_ids = _compute_in_triad_ids(output_dir, image_id, report_radius_um)
                sub = sub[sub["cell_id"].astype(str).isin(in_triad_ids)]

            if marker_col not in cells.columns or sub.empty:
                row[f"{name}__n"], row[f"{name}__mean"] = 0, np.nan
            else:
                vals = sub[marker_col].dropna().astype(float)
                row[f"{name}__n"] = len(vals)
                row[f"{name}__mean"] = float(vals.mean()) if len(vals) else np.nan
        per_image_rows.append(row)

    if not per_image_rows:
        return None
    per_image_df = pd.DataFrame(per_image_rows)

    agg_rows = []
    for patient_id, grp in per_image_df.groupby("patient_id"):
        row = {"patient_id": patient_id}
        for exp in marker_exposures_cfg:
            name = exp["name"]
            n_col, mean_col = f"{name}__n", f"{name}__mean"
            total_n = grp[n_col].sum()
            if total_n > 0:
                weighted = (grp[mean_col].fillna(0) * grp[n_col]).sum() / total_n
                row[name] = round(float(weighted), 4)
                row[f"{name}_n_cells"] = int(total_n)
            else:
                row[name] = np.nan
                row[f"{name}_n_cells"] = 0
        agg_rows.append(row)
    return pd.DataFrame(agg_rows)


# ── Annotation loading ──────────────────────────────────────────────────────

def _load_annotation(
    annot_file: str,
    patient_id_col: str,
    group_col: str,
    group_code_map: dict,
    experiment_groups: list,
) -> pd.DataFrame:
    """Loads the external annotation file. Keeps every column (no hardcoded
    clinical-field allowlist — the previous version only kept a fixed CRC
    TNM-staging-shaped set of columns); just standardizes patient_id to a
    string key for merging, and (if group_col is configured) adds an
    experiment_group_label column."""
    surv_df = pd.read_csv(annot_file)
    if patient_id_col not in surv_df.columns:
        raise ValueError(
            f"patient_id_col {patient_id_col!r} not found in annotation file columns: "
            f"{list(surv_df.columns)}"
        )
    surv_df = surv_df.copy()
    surv_df["patient_id"] = surv_df[patient_id_col].apply(
        lambda x: _normalize_code(x) if pd.notna(x) else None
    )

    if group_col and group_col in surv_df.columns:
        if group_code_map:
            gmap = {_normalize_code(k): v for k, v in group_code_map.items()}
        else:
            codes_present = sorted({_normalize_code(x) for x in surv_df[group_col].dropna().unique()})
            gmap = dict(zip(codes_present, experiment_groups))
            print(f"[survival]   ⚠️  No analysis.survival.group_code_map configured — "
                  f"guessing {group_col!r} code -> experiment_group by position: {gmap}. "
                  f"Set group_code_map explicitly in YAML if this is wrong.")
        surv_df["experiment_group_label"] = surv_df[group_col].apply(
            lambda x: gmap.get(_normalize_code(x)) if pd.notna(x) else None
        ).fillna("Unknown")

    return surv_df


# ── Kaplan-Meier ─────────────────────────────────────────────────────────────

def _km_plot(
    df: pd.DataFrame,
    duration_col: str,
    event_col: str,
    censor_is_one: bool,
    group_col: str,
    title: str,
    save_path: str,
) -> dict:
    """Kaplan-Meier curves for the groups in group_col. Returns a dict with
    per-group n/events and (if exactly 2 groups) a log-rank p-value."""
    df = df.dropna(subset=[duration_col, event_col, group_col]).copy()

    df["_event"] = (1 - df[event_col].astype(int)) if censor_is_one else df[event_col].astype(int)

    groups = sorted(df[group_col].unique(), key=str)
    if len(groups) < 2:
        print(f"  ⚠️  Only one group in {group_col} — cannot plot KM.")
        return {}

    palette = ["#d62728", "#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd"]
    colors = {g: palette[i % len(palette)] for i, g in enumerate(groups)}

    fig, ax = plt.subplots(figsize=(9, 6))
    results = {}

    for grp in groups:
        sub = df[df[group_col] == grp]
        kmf = KaplanMeierFitter()
        kmf.fit(sub[duration_col], event_observed=sub["_event"], label=f"{grp} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[grp], linewidth=2)
        results[grp] = {"n": len(sub), "events": int(sub["_event"].sum())}

    if len(groups) == 2:
        g1 = df[df[group_col] == groups[0]]
        g2 = df[df[group_col] == groups[1]]
        lr = logrank_test(
            g1[duration_col], g2[duration_col],
            event_observed_A=g1["_event"], event_observed_B=g2["_event"],
        )
        pval = lr.p_value
        sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
        ax.text(0.02, 0.08, f"Log-rank p = {pval:.3e}  {sig}",
                transform=ax.transAxes, fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
        results["logrank_p"] = pval
        results["significance"] = sig

    ax.set_xlabel("Time", fontsize=12)
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


# ── Cox PH (optional multivariate) ──────────────────────────────────────────

def _cox_ph(
    patient_df: pd.DataFrame,
    duration_col: str,
    event_col: str,
    censor_is_one: bool,
    covariate_cols: list,
    save_path: str,
):
    """Multivariate Cox proportional-hazards regression across the
    configured covariates. Long promised in this module's docstring, never
    implemented until this rewrite — real now, off by default (only runs if
    analysis.survival.covariates is non-empty)."""
    if not _HAS_COX:
        print("  ⚠️  lifelines.CoxPHFitter not available — skipping Cox PH "
              "(pip install/upgrade lifelines).")
        return None

    missing = [c for c in covariate_cols if c not in patient_df.columns]
    if missing:
        print(f"  ⚠️  Cox PH covariates not found in patient table: {missing} — skipping.")
        return None

    df = patient_df.dropna(subset=[duration_col, event_col] + covariate_cols).copy()
    if len(df) < max(10, 3 * len(covariate_cols)):
        print(f"  ⚠️  Only {len(df)} complete cases for {len(covariate_cols)} covariate(s) — "
              f"skipping Cox PH (too few events/patients relative to covariates).")
        return None

    df["_event"] = (1 - df[event_col].astype(int)) if censor_is_one else df[event_col].astype(int)

    cph_df = df[[duration_col, "_event"] + covariate_cols].copy()
    cat_cols = [c for c in covariate_cols if cph_df[c].dtype == object]
    if cat_cols:
        cph_df = pd.get_dummies(cph_df, columns=cat_cols, drop_first=True)

    cph = CoxPHFitter()
    try:
        cph.fit(cph_df, duration_col=duration_col, event_col="_event")
    except Exception as e:
        print(f"  ⚠️  Cox PH fit failed: {e} — skipping.")
        return None

    summary = cph.summary.reset_index().rename(columns={"index": "covariate"})
    summary.to_csv(save_path, index=False)
    print(f"  Cox PH ({len(df)} patients, {int(df['_event'].sum())} events):")
    print(cph.summary[["coef", "exp(coef)", "p"]].round(4).to_string())
    return summary


# ── Density / marker-exposure split helper ──────────────────────────────────

def _add_split(df: pd.DataFrame, value_col: str, split_col: str, split_by, pos_label: str, neg_label: str):
    """Adds a High/Low (or presence/absence) split column for value_col,
    matching the density-split logic: median split, with a fallback to
    presence/absence when the median is 0 (common when the exposure is rare)."""
    if split_by == "median":
        threshold = df[value_col].median()
        if threshold == 0:
            df[split_col] = np.where(df[value_col] > 0, pos_label, neg_label)
        else:
            df[split_col] = np.where(df[value_col] >= threshold, "High", "Low")
    else:
        try:
            threshold = float(split_by)
        except (TypeError, ValueError):
            threshold = df[value_col].median()
        df[split_col] = np.where(df[value_col] >= threshold, "High", "Low")
    return threshold


# ── Main entry point ─────────────────────────────────────────────────────────

def run_survival_analysis(cfg: dict) -> None:
    """
    Cohort-wide survival analysis — called by run_pipeline.py, after triads
    (and, if marker exposures are configured, after functional — both must
    have already written their per-image output to output_dir).
    """
    exp_name   = cfg["experiment"]["name"]
    output_dir = cfg["paths"]["output_dir"]

    s_cfg = cfg.get("analysis", {}).get("survival", {})
    if not s_cfg.get("enabled", False):
        print("[survival] disabled in config — skipping.")
        return

    annot_file         = s_cfg.get("patient_annotation_file", "")
    patient_id_col      = s_cfg.get("patient_id_col", "patient_id")

    # image_patient_map may be an inline dict (small cohorts) or a path to a
    # JSON file (TMA-scale cohorts where inlining ~100+ entries directly in
    # the YAML would be unwieldy) — see build_image_patient_map.py.
    image_patient_map_cfg = s_cfg.get("image_patient_map", {}) or {}
    if isinstance(image_patient_map_cfg, str):
        with open(image_patient_map_cfg) as f:
            image_patient_map = json.load(f)
        print(f"[survival] Loaded image_patient_map from {image_patient_map_cfg} "
              f"({len(image_patient_map)} entries)")
    else:
        image_patient_map = image_patient_map_cfg

    group_col           = s_cfg.get("group_col", "")
    group_code_map      = s_cfg.get("group_code_map", {})
    area_per_image      = s_cfg.get("area_per_image_um2", 1_000_000.0)
    split_by             = s_cfg.get("split_by", "median")
    censor_is_one        = s_cfg.get("censor_is_one", True)
    outcomes_cfg          = s_cfg.get("outcomes") or [
        {"name": "Overall Survival (OS)", "duration_col": "OS", "event_col": "OS_Censor"},
        {"name": "Disease-Free Survival (DFS)", "duration_col": "DFS", "event_col": "DFS_Censor"},
    ]
    marker_exposures_cfg = s_cfg.get("marker_exposures", []) or []
    report_radius_um      = s_cfg.get("report_radius_um", 20.0)
    covariates            = s_cfg.get("covariates", []) or []

    experiment_groups = cfg["experiment"]["groups"]

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
    print(f"[survival] Image->patient   : "
          f"{'image_patient_map (' + str(len(image_patient_map)) + ' entries)' if image_patient_map else '1:1 (image_id is the patient id)'}")
    print(f"[survival] Area/image       : {area_per_image:,.0f} µm²")
    print(f"[survival] Split by         : {split_by}")
    if marker_exposures_cfg:
        print(f"[survival] Marker exposures : {[m['name'] for m in marker_exposures_cfg]}")
    print()

    # ── Load inputs ──────────────────────────────────────────────────────
    triad_summary = pd.read_csv(summary_csv)
    surv_df = _load_annotation(annot_file, patient_id_col, group_col, group_code_map, experiment_groups)

    input_dir = cfg["paths"]["input_dir"]
    all_image_ids = []
    if os.path.exists(input_dir):
        for f in os.listdir(input_dir):
            if f.endswith("_matched_with_boundaries.csv") and not f.startswith("._"):
                all_image_ids.append(f.replace("_matched_with_boundaries.csv", ""))
    if not all_image_ids:
        all_image_ids = list(triad_summary["image_id"].unique())
        print("[survival]   ⚠️  Could not list input_dir — using only images with triads.")

    # ── Build patient-level table ────────────────────────────────────────
    patient_df = _build_patient_triad_density(triad_summary, all_image_ids, image_patient_map, area_per_image)

    marker_df = _build_patient_marker_exposure(
        output_dir, all_image_ids, image_patient_map, marker_exposures_cfg, report_radius_um
    )
    if marker_df is not None:
        patient_df = patient_df.merge(marker_df, on="patient_id", how="left")

    patient_df = patient_df.merge(surv_df, on="patient_id", how="left")

    outcome_duration_cols = [o["duration_col"] for o in outcomes_cfg]
    n_with_outcome = patient_df[outcome_duration_cols].notna().any(axis=1).sum() if outcome_duration_cols else 0
    print(f"[survival] Patients (subjects) : {len(patient_df)}")
    print(f"[survival] With outcome data   : {n_with_outcome}")
    print(f"[survival] With triads         : {(patient_df['total_triads'] > 0).sum()}")

    # ── Density split ────────────────────────────────────────────────────
    density_threshold = _add_split(patient_df, "triads_per_mm2", "density_group", split_by, "Triad+", "Triad−")
    print(f"[survival] Density split threshold : {density_threshold:.4f} triads/mm²")

    # ── Marker-exposure splits ───────────────────────────────────────────
    marker_thresholds = {}
    for exp in marker_exposures_cfg:
        name = exp["name"]
        if name in patient_df.columns:
            split_col = f"{name}_group"
            marker_thresholds[name] = _add_split(patient_df, name, split_col, split_by, f"{name}+", f"{name}−")

    patient_csv = os.path.join(surv_out, "patient_cohort_summary.csv")
    patient_df.to_csv(patient_csv, index=False)
    print(f"[survival] Patient cohort table saved → {patient_csv}")

    # ── KM + log-rank + optional Cox PH, per outcome ────────────────────
    lr_results = []
    for outcome in outcomes_cfg:
        name, dur_col, cens_col = outcome["name"], outcome["duration_col"], outcome["event_col"]
        if dur_col not in patient_df.columns:
            print(f"[survival] ⚠️  Column {dur_col!r} not found — skipping {name}.")
            continue

        # -- density split --
        save_path = os.path.join(surv_out, f"km_{dur_col}_density.png")
        print(f"\n[survival] Plotting {name} by triad density …")
        res = _km_plot(
            patient_df, duration_col=dur_col, event_col=cens_col, censor_is_one=censor_is_one,
            group_col="density_group",
            title=f"{exp_name} — {name}\nTriad density split (threshold={density_threshold:.3f})",
            save_path=save_path,
        )
        if "logrank_p" in res:
            print(f"  Log-rank p = {res['logrank_p']:.4e}  {res['significance']}")
            lr_results.append({"outcome": name, "split": "triad_density", "split_by": split_by,
                                "threshold": density_threshold, "logrank_p": res["logrank_p"],
                                "significance": res["significance"]})

        # -- marker-exposure splits --
        for exp in marker_exposures_cfg:
            mname = exp["name"]
            split_col = f"{mname}_group"
            if split_col not in patient_df.columns:
                continue
            save_path = os.path.join(surv_out, f"km_{dur_col}_{mname}.png")
            print(f"[survival] Plotting {name} by {mname} …")
            res = _km_plot(
                patient_df, duration_col=dur_col, event_col=cens_col, censor_is_one=censor_is_one,
                group_col=split_col,
                title=f"{exp_name} — {name}\n{mname} split (threshold={marker_thresholds.get(mname, float('nan')):.3f})",
                save_path=save_path,
            )
            if "logrank_p" in res:
                print(f"  Log-rank p = {res['logrank_p']:.4e}  {res['significance']}")
                lr_results.append({"outcome": name, "split": mname, "split_by": split_by,
                                    "threshold": marker_thresholds.get(mname), "logrank_p": res["logrank_p"],
                                    "significance": res["significance"]})

        # -- experiment_group split --
        if "experiment_group_label" in patient_df.columns:
            save_path = os.path.join(surv_out, f"km_{dur_col}_experiment_group.png")
            print(f"[survival] Plotting {name} by experiment_group …")
            res = _km_plot(
                patient_df, duration_col=dur_col, event_col=cens_col, censor_is_one=censor_is_one,
                group_col="experiment_group_label",
                title=f"{exp_name} — {name} by experiment group",
                save_path=save_path,
            )
            if "logrank_p" in res:
                print(f"  Log-rank p = {res['logrank_p']:.4e}  {res['significance']}")
                lr_results.append({"outcome": name, "split": "experiment_group", "split_by": "n/a",
                                    "threshold": np.nan, "logrank_p": res["logrank_p"],
                                    "significance": res["significance"]})

        # -- optional multivariate Cox PH --
        if covariates:
            print(f"[survival] Cox PH for {name} …")
            _cox_ph(patient_df, dur_col, cens_col, censor_is_one, covariates,
                    os.path.join(surv_out, f"cox_ph_{dur_col}.csv"))

    if lr_results:
        pd.DataFrame(lr_results).to_csv(os.path.join(surv_out, "survival_logrank_results.csv"), index=False)

    print(f"\n[survival] ✓  All outputs → {surv_out}")
