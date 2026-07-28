"""
spatia.analysis.cell_typing
============================
Generalizable cell typing for spatial proteomics data.
Supports two modes controlled by config:

  mode: "automatic"
    - GMM thresholds each marker
    - Maps marker combos → cell types via cell_type_definitions.yaml
    - No human input required

  mode: "semi_automatic"
    - GMM thresholds each marker
    - Runs PCA → KNN → Leiden clustering
    - STOPS and saves plots if cluster_labels_file is null/missing
    - Resumes and assigns labels once cluster_labels_file is filled in

Usage
-----
    from spatia.analysis.cell_typing import run_cell_typing
    run_cell_typing(cfg)   # cfg loaded from config YAML via yaml.safe_load

CD45 gating note
----------------
Semi-auto mode historically used cd45_std_multiplier=8 (very tight gate).
Auto mode used cd45_std_multiplier=3 (more permissive).
Both are now explicit config parameters so the choice is intentional and documented.
"""

import os
import yaml
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from sklearn.mixture import GaussianMixture

try:
    import scanpy as sc
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False
    print("[cell_typing] WARNING: scanpy not installed — semi_automatic mode unavailable.")


# ── Config helpers ────────────────────────────────────────────────────────────

def _load_cell_type_definitions(path: str) -> dict:
    with open(path) as f:
        raw = yaml.safe_load(f)
    return raw.get("cell_type_definitions", raw)


def _load_cluster_labels(path: str) -> dict:
    """Returns {cluster_id_str: cell_type_label} or None if file missing/incomplete."""
    if not path or not os.path.exists(path):
        return None
    with open(path) as f:
        data = yaml.safe_load(f)
    labels = data.get("cluster_labels", {})
    # Treat any "Unknown" / null as incomplete — warn but continue
    if not labels:
        return None
    return {str(k): v for k, v in labels.items()}


# ── GMM thresholding ──────────────────────────────────────────────────────────

def _gmm_threshold(values: np.ndarray, std_multiplier: float, n_components: int = 2,
                   random_state: int = 42, n_init: int = 1,
                   max_cells: int = 50_000) -> float:
    """
    Fit a 2-component GMM and return threshold = mean_low + std_multiplier * std_low.
    Returns the global 95th percentile as fallback if GMM fails.

    max_cells: if more cells than this, subsample for GMM fitting (threshold
               is still applied to all cells). Keeps runtime O(max_cells).
    """
    clean = values[np.isfinite(values)].reshape(-1, 1)
    if len(clean) < 10:
        return float(np.percentile(clean, 95))
    # Subsample for speed on large datasets
    if len(clean) > max_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(len(clean), size=max_cells, replace=False)
        fit_data = clean[idx]
    else:
        fit_data = clean
    try:
        gmm = GaussianMixture(n_components=n_components, random_state=random_state,
                              max_iter=200, n_init=n_init)
        gmm.fit(fit_data)
        means = gmm.means_.flatten()
        stds  = np.sqrt(gmm.covariances_.flatten())
        low_idx = int(np.argmin(means))
        return float(means[low_idx] + std_multiplier * stds[low_idx])
    except Exception as e:
        print(f"    [GMM] WARNING: {e} — using 95th-percentile fallback")
        return float(np.percentile(clean, 95))


def compute_marker_thresholds(adata, markers: list, std_multipliers: dict,
                               default_std: float = 2.0,
                               n_components: int = 2, random_state: int = 42,
                               n_init: int = 1, max_cells: int = 50_000) -> dict:
    """
    Returns {marker: threshold_value} for all markers present in adata.
    Uses per-marker std multiplier from std_multipliers dict, falling back to default_std.
    """
    thresholds = {}
    for m in markers:
        if m not in adata.var_names:
            print(f"    [GMM] marker '{m}' not in data — skipping")
            continue
        vals = adata[:, m].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray().flatten()
        else:
            vals = np.array(vals).flatten()
        mult = std_multipliers.get(m, default_std)
        t = _gmm_threshold(vals, mult, n_components=n_components,
                           random_state=random_state, n_init=n_init, max_cells=max_cells)
        thresholds[m] = t
        print(f"    {m:<20} std_mult={mult:.1f}  threshold={t:.4f}")
    return thresholds


def add_positivity_columns(adata, thresholds: dict) -> None:
    """
    Adds boolean columns '<marker>_pos' and intensity columns '<marker>_intensity'
    (0=negative, 1=+, 2=++, 3=+++) to adata.obs in-place.

    Uses a single pd.concat at the end to avoid DataFrame fragmentation.
    """
    new_cols: dict[str, np.ndarray] = {}

    for marker, threshold in thresholds.items():
        if marker not in adata.var_names:
            continue
        vals = adata[:, marker].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray().flatten()
        else:
            vals = np.array(vals).flatten()

        pos_mask = vals > threshold
        new_cols[f"{marker}_pos"] = pos_mask

        intensity = np.zeros(len(vals), dtype=np.int8)
        if pos_mask.sum() > 3:
            pos_vals = vals[pos_mask]
            q1 = np.percentile(pos_vals, 33)
            q2 = np.percentile(pos_vals, 66)
            intensity[pos_mask & (vals <= q1)] = 1
            intensity[pos_mask & (vals > q1) & (vals <= q2)] = 2
            intensity[pos_mask & (vals > q2)] = 3
        new_cols[f"{marker}_intensity"] = intensity

    # Assign all new columns in one shot — avoids DataFrame fragmentation
    new_df = pd.DataFrame(new_cols, index=adata.obs.index)
    adata.obs = pd.concat([adata.obs, new_df], axis=1)


# ── CD45 gating ───────────────────────────────────────────────────────────────

def gate_cd45_positive(adata, cd45_std_multiplier: float, n_components: int = 2,
                        random_state: int = 42, plot_dir: str = None):
    """
    Returns adata filtered to CD45+ cells.
    Saves a threshold histogram to plot_dir if provided.
    """
    if "CD45" not in adata.var_names:
        print("  [CD45 gate] WARNING: CD45 not found — returning all cells")
        return adata

    vals = adata[:, "CD45"].X
    if hasattr(vals, "toarray"):
        vals = vals.toarray().flatten()
    else:
        vals = np.array(vals).flatten()

    threshold = _gmm_threshold(vals, cd45_std_multiplier, n_components, random_state)
    print(f"  [CD45 gate] threshold={threshold:.4f}  (std_mult={cd45_std_multiplier})")

    if plot_dir:
        os.makedirs(plot_dir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(vals, bins=50, color="steelblue", alpha=0.7, label="All cells")
        ax.axvline(threshold, color="tomato", linewidth=2, label=f"Threshold={threshold:.3f}")
        ax.set_xlabel("CD45 expression")
        ax.set_ylabel("Cell count")
        ax.set_title(f"CD45 gate  (std_mult={cd45_std_multiplier})")
        ax.legend()
        plt.tight_layout()
        fig.savefig(os.path.join(plot_dir, "cd45_gate.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)

    mask = vals > threshold
    print(f"  [CD45 gate] {mask.sum():,} / {len(mask):,} cells pass ({mask.mean()*100:.1f}%)")
    return adata[mask].copy()


# ── Automatic cell typing ─────────────────────────────────────────────────────

def _score_cell_type(obs_row: pd.Series, rule: dict) -> int:
    """
    Score how well a cell matches one rule.
    Returns: 2 if required+preferred match, 1 if only required match, 0 if fails.
    """
    for req in rule.get("required", []):
        if not obs_row.get(req, False):
            return 0
    req_any = rule.get("required_any", [])
    if req_any and not any(obs_row.get(r, False) for r in req_any):
        return 0
    for exc in rule.get("excluded", []):
        if obs_row.get(exc, False):
            return 0
    preferred = rule.get("preferred", [])
    bonus = sum(1 for p in preferred if obs_row.get(p, False))
    return 2 + bonus if preferred else 1


def assign_cell_types_automatic(adata, cell_type_definitions: dict) -> pd.Series:
    """
    For each cell in adata.obs, score against every rule in cell_type_definitions
    and assign the best-matching cell type.

    Fully vectorised — no row-by-row Python loops. Matches the original
    scoring: base=1 (no preferred) or base=2 (preferred defined), +1 per
    preferred marker that is True.
    """
    obs    = adata.obs
    labels = pd.Series("Unassigned", index=obs.index, dtype=object)
    scores = pd.Series(0,            index=obs.index, dtype=np.int32)

    def _col(name: str) -> pd.Series:
        """Return boolean obs column, or all-False if missing."""
        if name in obs.columns:
            return obs[name].astype(bool)
        return pd.Series(False, index=obs.index)

    for cell_type, rules in cell_type_definitions.items():
        if isinstance(rules, dict):
            rules = [rules]
        for rule in rules:
            passes = pd.Series(True, index=obs.index)

            for col in rule.get("required", []):
                passes &= _col(col)

            req_any = rule.get("required_any", [])
            if req_any:
                any_true = pd.Series(False, index=obs.index)
                for col in req_any:
                    any_true |= _col(col)
                passes &= any_true

            for col in rule.get("excluded", []):
                passes &= ~_col(col)

            # Score: matches original _score_cell_type logic exactly
            preferred_list = rule.get("preferred", [])
            cell_scores = passes.astype(np.int32)   # 1 if passes, else 0
            if preferred_list:
                cell_scores += passes.astype(np.int32)  # bump base to 2 when preferred defined
                for col in preferred_list:
                    cell_scores += (passes & _col(col)).astype(np.int32)

            update_mask = cell_scores > scores
            labels[update_mask] = cell_type
            scores[update_mask] = cell_scores[update_mask]

    return labels


# ── Semi-automatic cell typing ────────────────────────────────────────────────

def run_clustering(adata, n_neighbors: int = 15, leiden_resolution: float = 0.5,
                   random_state: int = 42, plot_dir: str = None):
    """
    Runs PCA → KNN → UMAP → Leiden on adata.
    Saves UMAP plots to plot_dir.
    Returns adata with leiden cluster column added.
    """
    if not HAS_SCANPY:
        raise ImportError("scanpy is required for semi_automatic mode")

    print("  [clustering] Normalizing and log-transforming...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    sc.pp.pca(adata, random_state=random_state)

    print(f"  [clustering] KNN graph (k={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state)
    sc.tl.umap(adata, random_state=random_state)

    print(f"  [clustering] Leiden (resolution={leiden_resolution})...")
    sc.tl.leiden(adata, resolution=leiden_resolution, random_state=random_state,
                 key_added="leiden")

    n_clusters = adata.obs["leiden"].nunique()
    print(f"  [clustering] Found {n_clusters} clusters")

    if plot_dir:
        os.makedirs(plot_dir, exist_ok=True)
        # Leiden cluster UMAP
        fig, ax = plt.subplots(figsize=(8, 7))
        sc.pl.umap(adata, color="leiden", ax=ax, show=False,
                   title=f"Leiden clusters (res={leiden_resolution})")
        fig.savefig(os.path.join(plot_dir, "umap_leiden_clusters.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)

        # Condition UMAP
        if "condition" in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(8, 7))
            sc.pl.umap(adata, color="condition", ax=ax, show=False, title="Condition")
            fig.savefig(os.path.join(plot_dir, "umap_condition.png"),
                        dpi=150, bbox_inches="tight")
            plt.close(fig)

        # Cluster composition table
        comp = adata.obs.groupby("leiden").size().reset_index(name="n_cells")
        comp.to_csv(os.path.join(plot_dir, "cluster_sizes.csv"), index=False)

        # Print cluster sizes for reference
        print("\n  Cluster sizes (use these to fill cluster_labels.yaml):")
        print(comp.to_string(index=False))

    return adata


def assign_cluster_labels(adata, cluster_labels: dict) -> pd.Series:
    """Maps leiden cluster IDs to cell type labels via cluster_labels dict."""
    mapped = adata.obs["leiden"].map(cluster_labels)
    n_unmapped = mapped.isna().sum()
    if n_unmapped > 0:
        print(f"  [semi_auto] WARNING: {n_unmapped} cells in clusters not in cluster_labels — labeled 'Unassigned'")
        mapped = mapped.fillna("Unassigned")
    return mapped


# ── Main entry point ──────────────────────────────────────────────────────────

def run_cell_typing(cfg: dict) -> None:
    """
    Run cell typing for an experiment defined by cfg.

    Parameters
    ----------
    cfg : dict
        Loaded from config YAML (yaml.safe_load). Must have a 'cell_typing' key.

    Behavior by mode
    ----------------
    automatic:
        Runs GMM → positivity → cell type assignment → saves h5ad + plots.

    semi_automatic (cluster_labels_file = null):
        Runs GMM → positivity → PCA → Leiden → saves UMAP plots.
        STOPS here. Fill in cluster_labels_file, then re-run.

    semi_automatic (cluster_labels_file = path to filled yaml):
        Loads saved clustered h5ad → assigns labels → saves final h5ad + plots.
    """
    if not HAS_SCANPY:
        raise ImportError("scanpy is required for cell typing")

    ct_cfg   = cfg["cell_typing"]
    mode     = ct_cfg["mode"]           # "automatic" or "semi_automatic"
    exp_name = cfg["experiment"]["name"]

    input_file    = ct_cfg["input_file"]
    analysis_name = ct_cfg.get("analysis_name", exp_name)
    output_dir    = cfg["paths"]["output_dir"]
    plot_dir      = os.path.join(output_dir, "cell_typing_plots")
    data_dir      = os.path.join(output_dir, "cell_typing_data")
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)

    markers_cfg   = ct_cfg["markers"]
    panel         = markers_cfg["panel"]
    gating_only   = markers_cfg.get("gating_only", ["CD45"])

    gmm_cfg       = ct_cfg["gmm"]
    cd45_mult     = gmm_cfg["cd45_std_multiplier"]
    default_std   = gmm_cfg["default_std_multiplier"]
    per_marker    = gmm_cfg.get("per_marker_overrides", {})
    n_components  = gmm_cfg.get("n_components", 2)
    random_state  = gmm_cfg.get("random_state", 42)
    n_init        = gmm_cfg.get("n_init", 1)
    max_cells_gmm = gmm_cfg.get("max_cells_gmm", 50_000)

    # Merge per-marker overrides with default
    std_multipliers = {m: default_std for m in panel}
    std_multipliers.update(per_marker)
    for g in gating_only:
        std_multipliers[g] = cd45_mult

    clust_cfg     = ct_cfg.get("clustering", {})
    n_neighbors   = clust_cfg.get("n_neighbors", 15)
    leiden_res    = clust_cfg.get("leiden_resolution", 0.5)

    print(f"[cell_typing] Experiment : {exp_name}")
    print(f"[cell_typing] Mode       : {mode}")
    print(f"[cell_typing] Input      : {input_file}")

    # ── Load data ─────────────────────────────────────────────
    adata = sc.read(input_file)
    print(f"[cell_typing] Loaded {adata.n_obs:,} cells × {adata.n_vars} markers")

    # ── CD45 gate ─────────────────────────────────────────────
    # skip_cd45_gate: true  → run on all cells (needed when panel includes
    # both immune and non-immune cell types, e.g. CRC TMA).
    # Also auto-skipped when gating_only is empty.
    skip_gate = gmm_cfg.get("skip_cd45_gate", False) or not gating_only
    if skip_gate:
        print(f"[cell_typing] CD45 gate  : SKIPPED — running on all {adata.n_obs:,} cells")
        adata_cd45 = adata
    else:
        adata_cd45 = gate_cd45_positive(
            adata, cd45_mult, n_components=n_components,
            random_state=random_state, plot_dir=plot_dir
        )

    # ── Markers for analysis (exclude gating-only markers) ───
    markers_for_analysis = [m for m in panel if m not in gating_only and m in adata_cd45.var_names]
    missing = [m for m in panel if m not in gating_only and m not in adata_cd45.var_names]
    if missing:
        print(f"  [cell_typing] WARNING: markers not in data: {missing}")

    # ── GMM thresholds ────────────────────────────────────────
    print(f"\n  Computing GMM thresholds (n_init={n_init}, max_cells={max_cells_gmm:,})...")
    thresholds = compute_marker_thresholds(
        adata_cd45, markers_for_analysis, std_multipliers,
        default_std=default_std, n_components=n_components, random_state=random_state,
        n_init=n_init, max_cells=max_cells_gmm,
    )
    pd.Series(thresholds, name="threshold").to_csv(
        os.path.join(data_dir, "marker_thresholds.csv"), header=True
    )

    # ── Add positivity columns ────────────────────────────────
    add_positivity_columns(adata_cd45, thresholds)

    # ══════════════════════════════════════════════════════════
    # AUTOMATIC MODE
    # ══════════════════════════════════════════════════════════
    if mode == "automatic":
        defs_file = ct_cfg.get("cell_type_definitions_file")
        if not defs_file or not os.path.exists(defs_file):
            raise FileNotFoundError(
                f"cell_type_definitions_file not found: {defs_file}\n"
                "Set 'cell_typing.cell_type_definitions_file' in your config."
            )
        cell_type_defs = _load_cell_type_definitions(defs_file)
        print(f"\n  Assigning cell types from {len(cell_type_defs)} definitions...")
        adata_cd45.obs["cell_type"] = assign_cell_types_automatic(adata_cd45, cell_type_defs)

        _save_and_plot(adata_cd45, data_dir, plot_dir, analysis_name, mode)

    # ══════════════════════════════════════════════════════════
    # SEMI-AUTOMATIC MODE
    # ══════════════════════════════════════════════════════════
    elif mode == "semi_automatic":
        cluster_labels_file = ct_cfg.get("cluster_labels_file")
        clustered_h5ad = os.path.join(data_dir, f"{analysis_name}_clustered.h5ad")

        # ── Phase 1: clustering (run if cluster_labels not yet provided) ──
        if not cluster_labels_file or not os.path.exists(cluster_labels_file):
            print("\n  [semi_auto] Phase 1: Running clustering...")
            adata_cd45 = run_clustering(
                adata_cd45, n_neighbors=n_neighbors,
                leiden_resolution=leiden_res,
                random_state=random_state,
                plot_dir=plot_dir
            )
            adata_cd45.write(clustered_h5ad)
            print(f"\n  [semi_auto] ✓ Clustered data saved: {clustered_h5ad}")
            print(f"  [semi_auto] ✓ UMAP plots saved:     {plot_dir}/umap_leiden_clusters.png")
            print("\n" + "="*65)
            print("  NEXT STEP: Inspect the UMAP plots, then fill in:")
            print(f"  {cluster_labels_file or '<path/to/cluster_labels.yaml>'}")
            print("  Then set 'cell_typing.cluster_labels_file' in your config and re-run.")
            print("="*65)
            return  # ← intentional stop

        # ── Phase 2: label assignment ──
        print(f"\n  [semi_auto] Phase 2: Loading cluster labels from {cluster_labels_file}...")
        cluster_labels = _load_cluster_labels(cluster_labels_file)
        if cluster_labels is None:
            raise ValueError(
                f"cluster_labels_file exists but has no labels: {cluster_labels_file}\n"
                "Fill in the cluster_labels section and re-run."
            )

        # Load saved clustered h5ad if available, else re-cluster
        if os.path.exists(clustered_h5ad):
            print(f"  [semi_auto] Loading saved clustered data: {clustered_h5ad}")
            adata_cd45 = sc.read(clustered_h5ad)
        else:
            print("  [semi_auto] No saved clustered data found — re-running clustering...")
            adata_cd45 = run_clustering(
                adata_cd45, n_neighbors=n_neighbors,
                leiden_resolution=leiden_res,
                random_state=random_state,
                plot_dir=plot_dir
            )

        adata_cd45.obs["cell_type"] = assign_cluster_labels(adata_cd45, cluster_labels)
        _save_and_plot(adata_cd45, data_dir, plot_dir, analysis_name, mode)

    else:
        raise ValueError(f"Unknown cell_typing.mode: '{mode}'. Use 'automatic' or 'semi_automatic'.")


# ── Shared output logic ───────────────────────────────────────────────────────

def _save_and_plot(adata, data_dir, plot_dir, analysis_name, mode):
    out_h5ad = os.path.join(data_dir, f"{analysis_name}_cell_typed.h5ad")
    adata.write(out_h5ad)
    print(f"\n  ✓ Cell-typed h5ad: {out_h5ad}")

    # Cell type counts
    ct_counts = adata.obs["cell_type"].value_counts()
    ct_counts.to_csv(os.path.join(data_dir, "cell_type_counts.csv"), header=True)
    print("\n  Cell type distribution:")
    print(ct_counts.to_string())

    # Bar chart
    fig, ax = plt.subplots(figsize=(10, 5))
    ct_counts.plot(kind="bar", ax=ax, color="steelblue", edgecolor="white")
    ax.set_title(f"Cell Type Counts — {analysis_name} ({mode})", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("Cells")
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()
    fig.savefig(os.path.join(plot_dir, "cell_type_counts.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # UMAP colored by cell type (if UMAP exists)
    if "X_umap" in adata.obsm:
        import scanpy as sc
        fig, ax = plt.subplots(figsize=(9, 8))
        sc.pl.umap(adata, color="cell_type", ax=ax, show=False,
                   title=f"Cell Types — {analysis_name}")
        fig.savefig(os.path.join(plot_dir, "umap_cell_types.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)

    # Condition comparison if available
    if "condition" in adata.obs.columns:
        comp = (adata.obs.groupby(["condition", "cell_type"])
                .size().unstack(fill_value=0))
        comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100
        comp_pct.to_csv(os.path.join(data_dir, "cell_type_by_condition_pct.csv"))

        fig, ax = plt.subplots(figsize=(12, 5))
        comp_pct.T.plot(kind="bar", ax=ax, edgecolor="white")
        ax.set_title(f"Cell Type % by Condition — {analysis_name}", fontsize=12)
        ax.set_ylabel("% of cells")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(title="Condition", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        fig.savefig(os.path.join(plot_dir, "cell_type_by_condition.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"\n[cell_typing] Complete. Outputs in: {data_dir}")
