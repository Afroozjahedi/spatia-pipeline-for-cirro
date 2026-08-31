"""
spatia/analysis/validation.py
==============================
Post-step output validation for the SPATIA pipeline.

Each validator checks that a step's expected outputs:
  1. Exist on disk
  2. Are non-empty (file size > 0)
  3. Can be opened without error (for h5ad files)
  4. Contain expected content (obs columns, cell counts > 0)

Usage
-----
Called automatically by run_pipeline.py after each step.
Can also be run standalone to check a previous run:

    python -c "
    import yaml
    from spatia.analysis.validation import validate_step
    cfg = yaml.safe_load(open('config_example.yaml'))
    ok, errors = validate_step('preprocessing', cfg)
    print('OK' if ok else errors)
    "
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List, Tuple


# ─────────────────────────────────────────────────────────────────────────────
# CORE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

class ValidationError:
    def __init__(self, step: str, message: str, path: str = ""):
        self.step    = step
        self.message = message
        self.path    = path

    def __str__(self):
        loc = f"\n    Path: {self.path}" if self.path else ""
        return f"[{self.step}] {self.message}{loc}"


def _file_exists_and_nonempty(path: str) -> Tuple[bool, str]:
    """Return (True, '') or (False, reason)."""
    if not os.path.exists(path):
        return False, f"File not found: {path}"
    if os.path.getsize(path) == 0:
        return False, f"File is empty (0 bytes): {path}"
    return True, ""


def _open_h5ad(path: str) -> Tuple[bool, str, object]:
    """
    Try to open an h5ad file with anndata.
    Returns (success, error_message, adata_or_None).
    """
    try:
        import anndata as ad
        adata = ad.read_h5ad(path)
        return True, "", adata
    except Exception as exc:
        return False, f"Cannot open h5ad ({type(exc).__name__}): {exc}", None


def _check_obs_columns(adata, required_cols: List[str], path: str) -> List[str]:
    """Return list of missing columns in adata.obs."""
    missing = [c for c in required_cols if c not in adata.obs.columns]
    return missing


# ─────────────────────────────────────────────────────────────────────────────
# STEP VALIDATORS
# ─────────────────────────────────────────────────────────────────────────────

def validate_segmentation(cfg: dict) -> Tuple[bool, List[ValidationError]]:
    """
    Checks after run_segmentation():
    - segmentation_results_dir exists
    - At least one *_mesmer_result.csv was produced
    - Every masked ROI TIFF found in masked_roi_dir has a corresponding
      *_mesmer_result.csv. A shortfall means some images silently failed
      segmentation or were skipped for a channel mismatch --
      run_cell_segmentation() collects those in its own return value, but
      validate_step() here only receives cfg, not that return value (same as
      every other validator), so this re-derives the same signal from disk:
      TIFF count in vs. CSV count out.
    - Each *_mesmer_result.csv is non-empty and has the columns
      preprocessing.py's run_preprocessing() actually reads (label, x, y, area).

    Note: this does NOT re-run the channel-order check itself -- that already
    runs (and raises) inside run_cell_segmentation() when
    segmentation.channel_check is enabled. This validator's job is to catch
    the *other* failure mode: a generic per-image exception (e.g. a spacec/
    Mesmer crash) that run_cell_segmentation() deliberately does NOT raise on,
    so the pipeline doesn't halt over one bad image -- but also doesn't
    silently continue with fewer cells than expected without this check.
    """
    errors: List[ValidationError] = []
    step = "segmentation"

    paths = cfg.get("paths", {})
    seg_dir = paths.get("segmentation_results_dir")

    if not seg_dir or not os.path.isdir(seg_dir):
        errors.append(ValidationError(
            step,
            "segmentation_results_dir not found — did segmentation run?",
            seg_dir or "",
        ))
        return False, errors

    csv_files = []
    for root, _dirs, files in os.walk(seg_dir):
        for f in files:
            if f.endswith("_mesmer_result.csv"):
                csv_files.append(os.path.join(root, f))

    if not csv_files:
        errors.append(ValidationError(
            step, "No *_mesmer_result.csv files found.", seg_dir,
        ))
        return False, errors

    # Compare against masked ROI TIFFs actually available as input -- a
    # shortfall here is the disk-based signal that some image(s) failed or
    # were skipped, since we can't see run_cell_segmentation's in-memory
    # errors list from here.
    masked_dir = paths.get("masked_roi_dir") or os.path.join(
        paths.get("output_dir", "."), "masked_rois"
    )
    if os.path.isdir(masked_dir):
        tif_count = 0
        for root, _dirs, files in os.walk(masked_dir):
            for f in files:
                if f.endswith(".tif") and "preview" not in f:
                    tif_count += 1
        if tif_count and len(csv_files) < tif_count:
            errors.append(ValidationError(
                step,
                f"{tif_count} masked ROI TIFF(s) found in {masked_dir}, but "
                f"only {len(csv_files)} *_mesmer_result.csv produced "
                f"({tif_count - len(csv_files)} short). Check the segmentation "
                f"step's console log for 'FAILED segmentation' or 'CHANNEL "
                f"MISMATCH' entries for the missing image(s).",
                seg_dir,
            ))

    # Per-CSV content checks
    expected_cols = {"label", "x", "y", "area"}
    for path in csv_files:
        ok, msg = _file_exists_and_nonempty(path)
        if not ok:
            errors.append(ValidationError(step, msg, path))
            continue
        try:
            import pandas as pd
            df = pd.read_csv(path, nrows=5)
            missing = expected_cols - set(df.columns)
            if missing:
                errors.append(ValidationError(
                    step,
                    f"Missing expected columns: {sorted(missing)}",
                    path,
                ))
        except Exception as exc:
            errors.append(ValidationError(step, f"Cannot read CSV: {exc}", path))

    return len(errors) == 0, errors


def validate_preprocessing(cfg: dict) -> Tuple[bool, List[ValidationError]]:
    """
    Checks after run_preprocessing():
    - At least one *_combined_all_experiment_groups.h5ad exists in individual_processed_data/
    - Each combined h5ad opens cleanly
    - Each combined h5ad has a 'experiment_group' column in obs
    - Each combined h5ad has > 0 cells
    """
    errors: List[ValidationError] = []
    step = "preprocessing"

    base_out    = cfg["paths"]["output_dir"]
    tissues_dir = os.path.join(base_out, "combined_processed_data", "individual_processed_data")

    if not os.path.isdir(tissues_dir):
        errors.append(ValidationError(
            step,
            f"Output directory not found — did preprocessing run?",
            tissues_dir,
        ))
        return False, errors

    combined_files = [
        f for f in os.listdir(tissues_dir)
        if f.endswith("_combined_all_experiment_groups.h5ad")
    ]

    if not combined_files:
        errors.append(ValidationError(
            step,
            "No *_combined_all_experiment_groups.h5ad files found.",
            tissues_dir,
        ))
        return False, errors

    for fname in combined_files:
        path = os.path.join(tissues_dir, fname)

        ok, msg = _file_exists_and_nonempty(path)
        if not ok:
            errors.append(ValidationError(step, msg, path))
            continue

        ok, msg, adata = _open_h5ad(path)
        if not ok:
            errors.append(ValidationError(step, msg, path))
            continue

        if adata.n_obs == 0:
            errors.append(ValidationError(
                step, f"h5ad has 0 cells after processing.", path
            ))

        missing = _check_obs_columns(adata, ["experiment_group", "image_ID"], path)
        if missing:
            errors.append(ValidationError(
                step,
                f"Missing expected obs columns: {missing}",
                path,
            ))

        experiment_groups = cfg["experiment"]["groups"]
        if "experiment_group" in adata.obs.columns:
            found = set(adata.obs["experiment_group"].unique())
            expected = set(experiment_groups)
            if not expected.issubset(found) and len(found) < 2:
                errors.append(ValidationError(
                    step,
                    f"Expected experiment_groups {expected} but found {found} in obs['experiment_group'].",
                    path,
                ))

    return len(errors) == 0, errors


def validate_cell_typing(cfg: dict) -> Tuple[bool, List[ValidationError]]:
    """
    Checks after run_cell_typing():

    Automatic mode:
    - {analysis_name}_cell_typed.h5ad exists and opens
    - Has 'cell_type' column in obs
    - At least 2 distinct cell types assigned
    - Cell count > 0

    Semi-automatic, first run (cluster_labels_file is null):
    - {analysis_name}_clustered.h5ad exists (pipeline stopped intentionally)
    - Informs user: this is expected, fill cluster_labels_file and re-run

    Semi-automatic, second run:
    - Same as automatic mode checks
    """
    errors: List[ValidationError] = []
    step = "cell_typing"

    ct_cfg        = cfg.get("cell_typing", {})
    mode          = ct_cfg.get("mode", "automatic")
    exp_name      = cfg["experiment"]["name"]
    analysis_name = ct_cfg.get("analysis_name", exp_name)
    output_dir    = cfg["paths"]["output_dir"]
    data_dir      = os.path.join(output_dir, "cell_typing_data")

    if mode == "semi_automatic" and ct_cfg.get("cluster_labels_file") is None:
        # First semi-auto run — pipeline stops intentionally after clustering.
        # Validate that the clustered h5ad was saved, then stop (this is not an error).
        clustered_path = os.path.join(data_dir, f"{analysis_name}_clustered.h5ad")
        ok, msg = _file_exists_and_nonempty(clustered_path)
        if not ok:
            errors.append(ValidationError(
                step,
                "Semi-auto first run: clustered h5ad not found. "
                "Check that cell typing ran to completion.",
                clustered_path,
            ))
        else:
            # Not an error — inform caller via a special message
            print(
                f"\n[validation] Semi-auto first run complete.\n"
                f"  Clustered h5ad saved: {clustered_path}\n"
                f"  Next step: inspect UMAP plots in cell_typing_plots/, "
                f"fill in cluster_labels_file in your config, then re-run."
            )
        return len(errors) == 0, errors

    # Automatic mode or semi-auto second run — expect cell_typed h5ad
    typed_path = os.path.join(data_dir, f"{analysis_name}_cell_typed.h5ad")

    ok, msg = _file_exists_and_nonempty(typed_path)
    if not ok:
        errors.append(ValidationError(step, msg, typed_path))
        return False, errors

    ok, msg, adata = _open_h5ad(typed_path)
    if not ok:
        errors.append(ValidationError(step, msg, typed_path))
        return False, errors

    if adata.n_obs == 0:
        errors.append(ValidationError(
            step, "Cell-typed h5ad has 0 cells.", typed_path
        ))

    missing = _check_obs_columns(adata, ["cell_type", "experiment_group"], typed_path)
    if missing:
        errors.append(ValidationError(
            step,
            f"Missing expected obs columns: {missing}",
            typed_path,
        ))

    if "cell_type" in adata.obs.columns:
        n_types = adata.obs["cell_type"].nunique()
        if n_types < 2:
            errors.append(ValidationError(
                step,
                f"Only {n_types} distinct cell type(s) assigned — "
                "check GMM thresholds and cell_type_definitions.yaml.",
                typed_path,
            ))
        unassigned_pct = (
            (adata.obs["cell_type"] == "Unassigned").sum() / adata.n_obs * 100
            if "Unassigned" in adata.obs["cell_type"].values
            else 0
        )
        if unassigned_pct > 50:
            errors.append(ValidationError(
                step,
                f"{unassigned_pct:.1f}% of cells are Unassigned — "
                "GMM thresholds may be too strict or marker names may not match.",
                typed_path,
            ))

    return len(errors) == 0, errors


def validate_triads(cfg: dict) -> Tuple[bool, List[ValidationError]]:
    """
    Checks after run_triad_analysis():
    - triad_summary.csv exists and is non-empty
    - experiment_group_comparison_counts.csv exists
    - At least one *_triad_pairs.csv exists (at least one image had triads)
    - triad_summary.csv has expected columns
    """
    errors: List[ValidationError] = []
    step = "triads"

    output_dir = cfg["paths"]["output_dir"]

    # Required summary files
    required_files = [
        "triad_summary.csv",
        "experiment_group_comparison_counts.csv",
    ]
    for fname in required_files:
        path = os.path.join(output_dir, fname)
        ok, msg = _file_exists_and_nonempty(path)
        if not ok:
            errors.append(ValidationError(step, msg, path))

    # Check triad_summary columns
    summary_path = os.path.join(output_dir, "triad_summary.csv")
    if os.path.exists(summary_path) and os.path.getsize(summary_path) > 0:
        try:
            import pandas as pd
            df = pd.read_csv(summary_path)
            expected_cols = ["image_id", "experiment_group", "n_triads"]
            missing = [c for c in expected_cols if c not in df.columns]
            if missing:
                errors.append(ValidationError(
                    step,
                    f"triad_summary.csv missing columns: {missing}",
                    summary_path,
                ))
            if len(df) == 0:
                errors.append(ValidationError(
                    step,
                    "triad_summary.csv has no rows — no images were processed.",
                    summary_path,
                ))
        except Exception as exc:
            errors.append(ValidationError(
                step, f"Cannot read triad_summary.csv: {exc}", summary_path
            ))

    # At least one per-image triad pairs file
    pair_files = [
        f for f in os.listdir(output_dir)
        if f.endswith("_triad_pairs.csv")
    ] if os.path.isdir(output_dir) else []

    if not pair_files:
        errors.append(ValidationError(
            step,
            "No *_triad_pairs.csv files found — no triads were detected in any image. "
            "Check anchor_type / partner_type names match your cell_typed h5ad.",
            output_dir,
        ))

    return len(errors) == 0, errors


# ─────────────────────────────────────────────────────────────────────────────
# DISPATCH
# ─────────────────────────────────────────────────────────────────────────────

def _validate_passthrough(cfg: dict) -> Tuple[bool, List]:
    """Stub validator for steps that produce their own internal checks."""
    return True, []


_VALIDATORS = {
    "segmentation":  validate_segmentation,
    "preprocessing": validate_preprocessing,
    "cell_typing":   validate_cell_typing,
    "triads":        validate_triads,
    "functional":    _validate_passthrough,
    "survival":      _validate_passthrough,
}


def validate_step(step: str, cfg: dict) -> Tuple[bool, List[ValidationError]]:
    """
    Run the validator for a given pipeline step.

    Parameters
    ----------
    step : one of 'segmentation', 'preprocessing', 'cell_typing', 'triads',
           'functional', 'survival'
    cfg  : parsed YAML config dict

    Returns
    -------
    (passed: bool, errors: List[ValidationError])
    """
    if step not in _VALIDATORS:
        raise ValueError(f"No validator registered for step '{step}'. "
                         f"Known steps: {list(_VALIDATORS.keys())}")
    return _VALIDATORS[step](cfg)


def validate_all(cfg: dict, steps: List[str] = None) -> Tuple[bool, dict]:
    """
    Run validators for all (or a subset of) steps.

    Returns
    -------
    (all_passed: bool, results: {step: (passed, errors)})
    """
    steps = steps or list(_VALIDATORS.keys())
    results = {}
    for step in steps:
        if step in _VALIDATORS:
            results[step] = validate_step(step, cfg)
    all_passed = all(passed for passed, _ in results.values())
    return all_passed, results
