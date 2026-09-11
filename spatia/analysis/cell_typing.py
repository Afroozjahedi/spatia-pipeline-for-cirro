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

Threshold mode note (2026-08-19)
---------------------------------
Two ways to turn a fitted 2-component GMM into a positive/negative call,
selected via gmm.threshold_mode:

  "std_multiplier" (default -- unchanged prior behavior)
    threshold = mean_low + std_multiplier * std_low
    A cell is positive if its value is above this single scalar cutoff.
    Simple, but only looks at the negative component's mean/std -- the
    positive component's own shape and the relative size of the two
    populations never factor in.

  "posterior"
    A cell is positive if its GMM posterior probability of belonging to
    the higher-mean component exceeds gmm.confidence_level. This uses
    both components' full distributions (mean, std, and mixture weight),
    which is the more statistically complete way to ask "how sure are we
    this cell is positive." An equivalent scalar "effective threshold"
    (the value where the posterior crosses confidence_level) is still
    computed and reported in marker_thresholds.csv for continuity with
    existing plots/reports, but the actual per-cell positivity call in
    "posterior" mode uses the real posterior probability, not that
    scalar re-derived cutoff.

Both modes depend on the GMM's two components actually being separated.
A degenerate/collapsed component (e.g. a zero-inflated marker where the
"low" component collapses onto the zero spike) makes a posterior just as
unreliable as a std-multiplier threshold -- see gmm.transform below.
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
from typing import Optional, Tuple, Dict
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# GROUP_COLORS / _color_for: shared with preprocessing_summary.py so
# experiment_group coloring is consistent across the QC summary report and
# the cell-typing diagnostic plots. Falls back to matplotlib's default
# cycle if that module can't be imported for some reason.
try:
    from spatia.analysis.preprocessing_summary import GROUP_COLORS, _color_for
except Exception:
    GROUP_COLORS = {}
    def _color_for(group, groups_in_order):
        cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        idx = groups_in_order.index(group) if group in groups_in_order else 0
        return cycle[idx % len(cycle)]

try:
    import scanpy as sc
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False
    print("[cell_typing] WARNING: scanpy not installed — semi_automatic mode unavailable.")

# Diagnostic/reporting plots (GMM per-marker distributions, positivity/
# expression heatmaps, restyled counts chart, stacked-by-group chart, UMAP)
# -- separate module, added 2026-09-10 at Afrouz's request. Imported lazily
# inside run_cell_typing() (not here at module load) so a bug or missing
# dependency in the plotting module can never block import of cell_typing.py
# itself -- see the try/except around its call site below.


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


# ── Transforms (Q16 remediation + Yeo-Johnson) ─────────────────────────────────

def _yeojohnson_inverse(y: np.ndarray, lmbda: float) -> np.ndarray:
    """
    Inverse of scipy.stats.yeojohnson for a known/fitted lambda. scipy
    provides the forward transform (and can fit lambda for you) but has no
    public inverse function the way it does for Box-Cox (inv_boxcox), so
    this implements it directly from the Yeo-Johnson definition. Handles
    the y>=0 and y<0 branches (corresponding to x>=0 / x<0 in the forward
    transform) and the lambda==0 / lambda==2 special cases separately.
    """
    y = np.asarray(y, dtype=float)
    x = np.empty_like(y)

    pos = y >= 0
    if lmbda != 0:
        x[pos] = np.power(y[pos] * lmbda + 1.0, 1.0 / lmbda) - 1.0
    else:
        x[pos] = np.expm1(y[pos])

    neg = ~pos
    if lmbda != 2:
        x[neg] = 1.0 - np.power(1.0 - (2.0 - lmbda) * y[neg], 1.0 / (2.0 - lmbda))
    else:
        x[neg] = -np.expm1(-y[neg])

    return x


def _transform_values(values: np.ndarray, transform: str = "none",
                       arcsinh_cofactor: float = 5.0,
                       yeojohnson_lambda: Optional[float] = None) -> Tuple[np.ndarray, dict]:
    """
    Q16 remediation (decided by Afrouz, 2026-08-18): optionally transform
    marker intensities before GMM fitting, to fix threshold degeneracy on
    zero-inflated channels (see Day 6 log entry — 14/53 CRC markers have
    30-50% exact-zero raw intensity, causing the GMM's "low" component to
    collapse onto the zero spike and produce threshold ~= 0). This is
    standard cytometry practice -- transform to a scale where the two
    populations are closer to Gaussian, fit the mixture there, transform
    the result back -- not a workaround.

    transform: "none" (legacy behavior, unchanged) | "arcsinh" | "log1p"
               | "yeojohnson".

    Returns (transformed_values, transform_params). transform_params is
    empty for "none"/"arcsinh"/"log1p" (nothing needs to be remembered to
    invert them). For "yeojohnson" it contains {"lambda": float} -- the
    fitted (or passed-through, if yeojohnson_lambda was given) power
    parameter, needed to invert the transform and worth reporting
    per-marker for transparency (see compute_marker_thresholds).

    arcsinh_cofactor: standard cytometry practice divides by a cofactor before
    taking asinh (5.0 is the common CyTOF default; empirically checked against
    this project's real CRC intensities on 2026-08-18 — median non-zero value
    for the most zero-inflated markers (CD138, CK) is ~6-7, so cofactor=5 keeps
    the dim-positive population in a reasonably spread part of the transformed
    range rather than compressing it against zero). This default has NOT been
    tuned per-marker or per-panel — treat as a starting point, not a validated
    optimum, and re-check if applying this to a new panel/platform.
    log1p is only defined for values >= -1; only use it on non-negative,
    non-z-scored intensities (true of the CRC h5ad, not guaranteed for other
    inputs — arcsinh is the safer default since it's defined everywhere).

    yeojohnson (new): unlike arcsinh/log1p (fixed-shape transforms), this
    fits a lambda parameter per marker from the data itself via
    scipy.stats.yeojohnson, so the transform's shape adapts to how skewed
    that specific marker's distribution actually is (lambda=0 reduces to
    log1p-like behavior; lambda=1 is close to identity; other values
    interpolate). Also the only option here that's well-defined for
    negative values (arcsinh handles them fine too, but log1p does not),
    which matters if this is ever run on z-scored rather than raw
    intensities.
    """
    if transform == "none":
        return values, {}
    elif transform == "arcsinh":
        return np.arcsinh(values / arcsinh_cofactor), {}
    elif transform == "log1p":
        return np.log1p(values), {}
    elif transform == "yeojohnson":
        from scipy.stats import yeojohnson
        flat = np.asarray(values, dtype=float).flatten()
        if yeojohnson_lambda is None:
            transformed, lam = yeojohnson(flat)
        else:
            transformed = yeojohnson(flat, lmbda=yeojohnson_lambda)
            lam = yeojohnson_lambda
        return transformed.reshape(np.asarray(values).shape), {"lambda": float(lam)}
    else:
        raise ValueError(
            f"Unknown gmm.transform: '{transform}'. Use 'none', 'arcsinh', 'log1p', or 'yeojohnson'."
        )


def _inverse_transform_scalar(value_t: float, transform: str, arcsinh_cofactor: float,
                               transform_params: dict) -> float:
    """Inverse of _transform_values for a single scalar (e.g. a fitted threshold)."""
    if transform == "none":
        return float(value_t)
    elif transform == "arcsinh":
        return float(np.sinh(value_t) * arcsinh_cofactor)
    elif transform == "log1p":
        return float(np.expm1(value_t))
    elif transform == "yeojohnson":
        lam = transform_params.get("lambda")
        if lam is None:
            raise ValueError("yeojohnson inverse requires transform_params['lambda']")
        return float(_yeojohnson_inverse(np.array([value_t]), lam)[0])
    else:
        raise ValueError(f"Unknown transform: '{transform}'")


# ── GMM fitting (shared by std_multiplier and posterior modes) ────────────────

def _detect_point_mass(values: np.ndarray, min_frac: float = 0.05) -> Optional[float]:
    """
    Detect a single repeated constant value that makes up at least
    `min_frac` of `values` -- e.g. a hardware/detector floor sitting at
    exactly 0 in raw fluorescence intensity, OR that same floor after a
    linear normalization (z-score, etc.) has shifted it to some other
    constant (mean and std are per-marker, so the floor value differs by
    marker; it is emphatically not 0 anymore).

    Added 2026-09-11 to generalize _fit_gmm's zero-exclusion (previously
    hardcoded to `values == 0`) so it still finds the point mass on data
    that already went through spatia/analysis/preprocessing.py's per-image
    z-score normalization (sp.pp.format(..., method="zscore")) before
    reaching cell_typing -- which is what the REAL pipeline run does,
    unlike the crc_tma_celltyping.yaml validation dataset (raw CSV
    intensities, no normalization, where the floor genuinely is 0 and the
    old check happened to be correct by coincidence). A linear transform
    moves a point mass to a new constant, it does not remove it -- the
    same GMM-collapse failure mode this function's zero-exclusion exists
    to prevent (see docstring below) would otherwise resurface silently
    on z-scored real data because `== 0` would simply stop matching it.

    Returns the detected value, or None if no single value reaches
    min_frac (i.e. no meaningful point mass -- most markers on already
    zero-excluded/normalized data). Values are rounded to 6 decimals first
    so floating-point noise doesn't split one truly-repeated value (e.g.
    an exact 0.0 written by multiple upstream tools) into near-duplicate
    bins that individually fall under min_frac.
    """
    if len(values) == 0:
        return None
    rounded = np.round(values, 6)
    vals, counts = np.unique(rounded, return_counts=True)
    idx = np.argmax(counts)
    frac = counts[idx] / len(values)
    return float(vals[idx]) if frac >= min_frac else None


def _fit_gmm(values: np.ndarray, n_components: int = 2, random_state: int = 42,
             n_init: int = 1, max_cells: int = 50_000, transform: str = "none",
             arcsinh_cofactor: float = 5.0, yeojohnson_lambda: Optional[float] = None) -> dict:
    """
    Fit a GMM once on (optionally transformed) marker values. Shared by
    both threshold_mode paths so a marker's GMM is only ever fit a single
    time, regardless of which mode consumes it.

    Returns a dict:
        "gmm"              : fitted GaussianMixture, or None if fitting
                              wasn't possible (too few finite cells) or failed
        "transform_params" : see _transform_values
        "low_idx"           : component index with the lower mean ("negative")
        "high_idx"          : component index with the higher mean ("positive")
        "clean"             : raw (untransformed) finite values, for the
                               percentile fallback when gmm is None

    Zero-inflation note: when transform != "none", the dominant point-mass
    value is excluded from the GMM fit (still included in "clean" for the
    fallback) -- see _detect_point_mass(). A monotonic transform
    (arcsinh/log1p/yeojohnson) reshapes smooth right-skew, but it cannot
    fix a genuine point mass: every one of these transforms maps a
    constant to a constant, so a real "30-50% of cells sit at the same
    floor value" spike (documented for 14/53 CRC markers, using raw
    fluorescence where that floor is exactly 0) survives the transform
    unchanged and can still make one GMM component collapse onto it,
    exactly the degeneracy the transform was meant to fix. Excluding the
    point mass and fitting the 2-component GMM on the remaining continuum
    is the standard complement to transforming.

    Generalized 2026-09-11 from a hardcoded `values == 0` check to
    _detect_point_mass(), which finds whatever constant is actually
    repeated -- because on data that already went through
    preprocessing.py's per-image z-score normalization before reaching
    cell_typing (the real pipeline's path, as opposed to the
    crc_tma_celltyping.yaml validation dataset's raw un-normalized CSV),
    a floor that was exactly 0 pre-normalization is shifted to some other
    marker-specific constant post-normalization -- a linear transform
    moves a point mass, it does not remove it. The old `== 0` check would
    silently stop finding it in that case, reintroducing this exact
    degeneracy on real data even with this exclusion "on".
    transform="none" is left untouched (no point-mass exclusion) to keep
    that path exactly backward compatible with pre-2026-08-19 behavior.
    Confirmed with a standalone test: on a 33%-exact-zero, right-skewed
    synthetic marker, yeojohnson WITHOUT this exclusion scored 52.6%
    agreement with ground truth (worse than doing nothing) because the
    fitted "low" component collapsed onto the zero spike and the
    resulting threshold landed near zero, calling almost everything
    positive; WITH zero exclusion it should recover discrimination
    between the real negative and positive populations.
    """
    clean = values[np.isfinite(values)].reshape(-1, 1)
    result = {"gmm": None, "transform_params": {}, "low_idx": None,
              "high_idx": None, "clean": clean}
    if len(clean) < 10:
        return result

    if transform != "none":
        point_mass = _detect_point_mass(clean.flatten())
        if point_mass is not None:
            fit_input = clean[np.abs(clean.flatten() - point_mass) > 1e-9].reshape(-1, 1)
        else:
            fit_input = clean
        if len(fit_input) < 10:
            print(f"    [GMM] WARNING: fewer than 10 values after excluding the point mass "
                  f"({point_mass}) for transform='{transform}' — falling back to percentile threshold")
            return result
    else:
        fit_input = clean

    clean_t, tparams = _transform_values(fit_input, transform, arcsinh_cofactor, yeojohnson_lambda)
    result["transform_params"] = tparams

    if len(clean_t) > max_cells:
        rng = np.random.default_rng(random_state)
        idx = rng.choice(len(clean_t), size=max_cells, replace=False)
        fit_data = clean_t[idx]
    else:
        fit_data = clean_t

    try:
        gmm = GaussianMixture(n_components=n_components, random_state=random_state,
                              max_iter=200, n_init=n_init)
        gmm.fit(fit_data)
        means = gmm.means_.flatten()
        result.update({
            "gmm": gmm,
            "low_idx": int(np.argmin(means)),
            "high_idx": int(np.argmax(means)),
        })
    except Exception as e:
        print(f"    [GMM] WARNING: {e} — GMM fit failed")

    return result


def _gmm_threshold_from_fit(fit: dict, std_multiplier: float, transform: str,
                             arcsinh_cofactor: float) -> float:
    """std_multiplier-mode threshold: mean_low + std_multiplier * std_low, on
    the raw (inverse-transformed) scale. Falls back to the 95th percentile
    of raw values if the GMM couldn't be fit."""
    if fit["gmm"] is None:
        return float(np.percentile(fit["clean"], 95)) if len(fit["clean"]) else float("nan")
    gmm = fit["gmm"]
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_.flatten())
    low_idx = fit["low_idx"]
    threshold_t = means[low_idx] + std_multiplier * stds[low_idx]
    return _inverse_transform_scalar(threshold_t, transform, arcsinh_cofactor, fit["transform_params"])


def _gmm_posterior_effective_threshold(fit: dict, confidence_level: float, transform: str,
                                        arcsinh_cofactor: float, n_grid: int = 2000) -> float:
    """
    posterior-mode "effective threshold": the raw-scale value at which the
    posterior probability of positive-component membership first crosses
    confidence_level, scanning from low to high. Found by a numerical grid
    search over the fit data's own transformed range rather than solving
    the two-Gaussian crossing point analytically -- that equation can have
    0, 1, or 2 real roots depending on the two components' relative
    variances, and a grid search sidesteps picking the "right" root.

    The grid search starts at the LOW component's own mean, not the
    minimum observed value. Two unequal-variance Gaussians can cross twice
    -- if the high (positive) component also happens to be wider, its
    slower-decaying tail can make the posterior spuriously high again far
    below the low component's mean, which is never a region anyone would
    consider "positive." Starting the search there avoids ever reporting
    an effective threshold from that spurious low-tail crossing.
    Confirmed with a standalone test: a mixture with means (2, 8) and
    stds (0.5, 1.5) has 280 grid points below x=0.32 where posterior
    non-monotonically flips back toward "positive" -- none of which are
    a plausible real cell value, but all of which would corrupt a naive
    grid search starting at the data minimum.

    This value is NOT what actually decides positivity in "posterior"
    mode (predict_proba is, per-cell -- see add_positivity_columns, which
    applies the same low-component-mean floor for the same reason) -- it
    exists so marker_thresholds.csv and existing threshold-based plots
    still get a single reportable number, for continuity with
    "std_multiplier" mode's output shape.

    Falls back to the 95th percentile of raw values if the GMM couldn't
    be fit, same as the std_multiplier path.
    """
    if fit["gmm"] is None:
        return float(np.percentile(fit["clean"], 95)) if len(fit["clean"]) else float("nan")

    gmm = fit["gmm"]
    high_idx = fit["high_idx"]
    low_idx = fit["low_idx"]
    clean = fit["clean"]
    clean_t, _ = _transform_values(clean, transform, arcsinh_cofactor,
                                    fit["transform_params"].get("lambda"))
    low_mean_t = float(gmm.means_.flatten()[low_idx])
    hi = float(np.max(clean_t))
    lo = min(low_mean_t, hi)
    if lo >= hi:
        return _inverse_transform_scalar(hi, transform, arcsinh_cofactor, fit["transform_params"])

    grid_t = np.linspace(lo, hi, n_grid).reshape(-1, 1)
    proba_high = gmm.predict_proba(grid_t)[:, high_idx]
    above = np.where(proba_high >= confidence_level)[0]
    if len(above) == 0:
        # Posterior never reaches confidence_level anywhere in the observed
        # range -- report the max observed value as a (very conservative)
        # stand-in, rather than raising, so batch runs don't halt on one
        # under-separated marker. This is exactly the kind of case worth a
        # manual look (see Notes/risks in the doc).
        crossing_t = hi
    else:
        crossing_t = float(grid_t[above[0], 0])

    return _inverse_transform_scalar(crossing_t, transform, arcsinh_cofactor, fit["transform_params"])


def compute_marker_thresholds(adata, markers: list, std_multipliers: dict,
                               default_std: float = 2.0,
                               n_components: int = 2, random_state: int = 42,
                               n_init: int = 1, max_cells: int = 50_000,
                               transform: str = "none",
                               arcsinh_cofactor: float = 5.0,
                               threshold_mode: str = "std_multiplier",
                               confidence_level: float = 0.8,
                               confidence_overrides: Optional[dict] = None,
                               column_map: Optional[dict] = None,
                               transform_overrides: Optional[dict] = None) -> Tuple[dict, dict]:
    """
    Fits a GMM per marker and returns (thresholds, fit_info).

    thresholds : {marker: threshold_value} -- always populated regardless
        of threshold_mode, for backward-compatible reporting
        (marker_thresholds.csv, threshold-line plots, etc). In
        "std_multiplier" mode this is the actual decision boundary. In
        "posterior" mode it's the equivalent "effective threshold" (see
        _gmm_posterior_effective_threshold) -- informative, but NOT what
        add_positivity_columns uses to decide positivity for that marker.

    fit_info : {marker: {"gmm", "transform_params", "low_idx", "high_idx",
        "threshold_mode", "confidence_level", "transform"}} -- everything
        add_positivity_columns needs to make the actual per-cell call,
        including the fitted GMM itself for "posterior" mode. "transform"
        (added 2026-09-11) is the transform actually used to fit THIS
        marker -- transform_overrides.get(m, transform) -- so downstream
        consumers (add_positivity_columns' posterior-mode path,
        plot_marker_distributions, plot_marker_umaps) never need the
        global `transform` value passed in separately to stay correct
        when different markers use different transforms.

    confidence_overrides: per-marker gmm.confidence_level overrides,
        mirroring how std_multipliers/per_marker_overrides already work --
        only used when threshold_mode == "posterior".

    transform_overrides: {marker: transform_name}, optional. Added
        2026-09-11 so an individual marker can be fit in a different
        transform space than the global `transform` default -- e.g. a
        zero-inflated marker (see docs/spatia_analysis_cell_typing.md's
        marker-shape-diagnosis notes) can use "arcsinh" while the rest of
        the panel stays "none", without forcing one transform on all 55
        markers. A marker not listed here falls back to `transform`,
        identical to prior behavior for any config that doesn't set this.

    column_map: {short_marker_name: actual_adata_column_name}, optional.
        Added 2026-09-10 for panels where the config's marker names (and
        cell_type_definitions.yaml's *_pos rules, and gmm.per_marker_overrides)
        use clean short names like "CD45", but the real column in this
        dataset's adata is something like "CD45 - hematopoietic cells (C14)"
        (CODEX-style export with descriptions + cycle suffixes -- see
        docs/spatia_analysis_cell_typing.md Notes/risks). Every dict key
        below (thresholds, fit_info, and later new_cols in
        add_positivity_columns) stays keyed by the SHORT name `m` throughout
        -- only the adata lookup itself is redirected to the real column via
        column_map.get(m, m). A config with no column_map (the previous,
        still-default behavior) resolves every marker to itself, identical
        to before this parameter existed.
    """
    confidence_overrides = confidence_overrides or {}
    column_map = column_map or {}
    transform_overrides = transform_overrides or {}
    thresholds: Dict[str, float] = {}
    fit_info: Dict[str, dict] = {}

    for m in markers:
        real_col = column_map.get(m, m)
        if real_col not in adata.var_names:
            print(f"    [GMM] marker '{m}' (column '{real_col}') not in data — skipping")
            continue
        vals = adata[:, real_col].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray().flatten()
        else:
            vals = np.array(vals).flatten()

        marker_transform = transform_overrides.get(m, transform)
        fit = _fit_gmm(vals, n_components=n_components, random_state=random_state,
                       n_init=n_init, max_cells=max_cells, transform=marker_transform,
                       arcsinh_cofactor=arcsinh_cofactor)

        mult = std_multipliers.get(m, default_std)
        conf = confidence_overrides.get(m, confidence_level)
        override_note = "  [override]" if m in transform_overrides else ""

        if threshold_mode == "posterior":
            t = _gmm_posterior_effective_threshold(fit, conf, marker_transform, arcsinh_cofactor)
            lam_note = f"  lambda={fit['transform_params']['lambda']:.3f}" if "lambda" in fit["transform_params"] else ""
            print(f"    {m:<20} mode=posterior  confidence={conf:.2f}  transform={marker_transform}"
                  f"  effective_threshold={t:.4f}{lam_note}{override_note}")
        else:
            t = _gmm_threshold_from_fit(fit, mult, marker_transform, arcsinh_cofactor)
            lam_note = f"  lambda={fit['transform_params']['lambda']:.3f}" if "lambda" in fit["transform_params"] else ""
            print(f"    {m:<20} mode=std_multiplier  std_mult={mult:.1f}  transform={marker_transform}"
                  f"  threshold={t:.4f}{lam_note}{override_note}")

        thresholds[m] = t
        fit_info[m] = {
            **fit,
            "threshold_mode": threshold_mode,
            "confidence_level": conf,
            "transform": marker_transform,
        }

    return thresholds, fit_info


def add_positivity_columns(adata, thresholds: dict, fit_info: Optional[dict] = None,
                            transform: str = "none", arcsinh_cofactor: float = 5.0,
                            column_map: Optional[dict] = None) -> None:
    """
    Adds boolean columns '<marker>_pos' and intensity columns '<marker>_intensity'
    (0=negative, 1=+, 2=++, 3=+++) to adata.obs in-place. In "posterior" mode
    also adds '<marker>_posterior' with the raw per-cell posterior probability,
    for transparency/debugging (e.g. spotting cells that are borderline).

    Uses a single pd.concat at the end to avoid DataFrame fragmentation.

    fit_info: from compute_marker_thresholds. Required for "posterior" mode
        (needs the fitted GMM); ignored for "std_multiplier" mode, which
        only needs the scalar thresholds dict (kept as a required, simple
        argument so any external caller with just thresholds still works
        exactly as before this feature was added).

    column_map: same {short_name: real_column_name} resolution as
        compute_marker_thresholds() (added 2026-09-10) -- new_cols is still
        keyed as "<marker>_pos"/"<marker>_intensity" using the SHORT name
        `marker` (from `thresholds`, which is itself keyed by short names),
        so cell_type_definitions.yaml's rules (which reference short-name
        "*_pos" columns) keep working unchanged. Only the adata lookup is
        redirected.
    """
    fit_info = fit_info or {}
    column_map = column_map or {}
    new_cols: dict[str, np.ndarray] = {}

    for marker, threshold in thresholds.items():
        real_col = column_map.get(marker, marker)
        if real_col not in adata.var_names:
            continue
        vals = adata[:, real_col].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray().flatten()
        else:
            vals = np.array(vals).flatten()

        info = fit_info.get(marker, {})
        mode = info.get("threshold_mode", "std_multiplier")

        if mode == "posterior" and info.get("gmm") is not None:
            gmm = info["gmm"]
            high_idx = info["high_idx"]
            low_idx = info["low_idx"]
            tparams = info.get("transform_params", {})
            # Use THIS marker's actual fitted transform (added 2026-09-11
            # alongside per-marker transform_overrides in
            # compute_marker_thresholds) rather than the function's global
            # `transform` default -- otherwise a marker fit with an
            # override would have its GMM evaluated against values
            # transformed the wrong way here.
            marker_transform = info.get("transform", transform)
            vals_t, _ = _transform_values(vals.reshape(-1, 1), marker_transform, arcsinh_cofactor,
                                           tparams.get("lambda"))
            posterior = gmm.predict_proba(vals_t)[:, high_idx]
            # Floor: never call a cell positive below the negative
            # component's own mean, even if the posterior says so -- see
            # _gmm_posterior_effective_threshold's docstring for why this
            # can happen with unequal-variance components (the positive
            # component's wider tail can spuriously outweigh the negative
            # component's narrower tail far below where any real negative
            # cell would sit).
            low_mean_t = float(gmm.means_.flatten()[low_idx])
            above_low_mean = vals_t.flatten() >= low_mean_t
            pos_mask = (posterior >= info.get("confidence_level", 0.8)) & above_low_mean
            new_cols[f"{marker}_posterior"] = posterior.astype(np.float32)
        else:
            # std_multiplier mode, or posterior mode with a GMM that failed
            # to fit (fit_info["gmm"] is None) -- same scalar-threshold
            # fallback compute_marker_thresholds itself uses in that case,
            # so behavior stays consistent instead of silently producing
            # an all-False positivity column.
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
                        random_state: int = 42, plot_dir: str = None,
                        transform: str = "none", arcsinh_cofactor: float = 5.0,
                        threshold_mode: str = "std_multiplier",
                        confidence_level: float = 0.8,
                        cd45_col: str = "CD45"):
    """
    Returns adata filtered to CD45+ cells.
    Saves a threshold histogram to plot_dir if provided.

    threshold_mode/confidence_level: same std_multiplier/posterior choice
    as the per-marker thresholds (see module docstring), applied here for
    consistency rather than leaving CD45 gating permanently on the older
    std_multiplier-only path while other markers can use posteriors.

    transform/arcsinh_cofactor: Q16 remediation, see `_transform_values()`.

    cd45_col: the actual adata column name for CD45 (added 2026-09-10).
        Defaults to the literal "CD45" -- unchanged prior behavior for any
        panel where that's already the real column name. Pass
        markers.column_map.get("CD45", "CD45") from run_cell_typing() for a
        panel (like this CRC one) whose real CD45 column is something like
        "CD45 - hematopoietic cells (C14)".
    """
    if cd45_col not in adata.var_names:
        print(f"  [CD45 gate] WARNING: '{cd45_col}' not found — returning all cells")
        return adata

    vals = adata[:, cd45_col].X
    if hasattr(vals, "toarray"):
        vals = vals.toarray().flatten()
    else:
        vals = np.array(vals).flatten()

    fit = _fit_gmm(vals, n_components=n_components, random_state=random_state,
                   transform=transform, arcsinh_cofactor=arcsinh_cofactor)

    if threshold_mode == "posterior":
        threshold = _gmm_posterior_effective_threshold(fit, confidence_level, transform, arcsinh_cofactor)
        print(f"  [CD45 gate] mode=posterior  confidence={confidence_level:.2f}  "
              f"effective_threshold={threshold:.4f}")
        if fit["gmm"] is not None:
            tparams = fit["transform_params"]
            vals_t, _ = _transform_values(vals.reshape(-1, 1), transform, arcsinh_cofactor,
                                           tparams.get("lambda"))
            posterior = fit["gmm"].predict_proba(vals_t)[:, fit["high_idx"]]
            # Same low-component-mean floor as add_positivity_columns --
            # see _gmm_posterior_effective_threshold's docstring.
            low_mean_t = float(fit["gmm"].means_.flatten()[fit["low_idx"]])
            above_low_mean = vals_t.flatten() >= low_mean_t
            mask = (posterior >= confidence_level) & above_low_mean
        else:
            mask = vals > threshold
    else:
        threshold = _gmm_threshold_from_fit(fit, cd45_std_multiplier, transform, arcsinh_cofactor)
        print(f"  [CD45 gate] threshold={threshold:.4f}  (std_mult={cd45_std_multiplier})")
        mask = vals > threshold

    if plot_dir:
        os.makedirs(plot_dir, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.hist(vals, bins=50, color="steelblue", alpha=0.7, label="All cells")
        ax.axvline(threshold, color="tomato", linewidth=2, label=f"Threshold={threshold:.3f}")
        ax.set_xlabel("CD45 expression")
        ax.set_ylabel("Cell count")
        title_suffix = (f"confidence={confidence_level}" if threshold_mode == "posterior"
                         else f"std_mult={cd45_std_multiplier}")
        ax.set_title(f"CD45 gate  ({threshold_mode}, {title_suffix})")
        ax.legend()
        plt.tight_layout()
        fig.savefig(os.path.join(plot_dir, "cd45_gate.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)

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

        # Experiment group UMAP
        if "experiment_group" in adata.obs.columns:
            fig, ax = plt.subplots(figsize=(8, 7))
            sc.pl.umap(adata, color="experiment_group", ax=ax, show=False, title="Experiment Group")
            fig.savefig(os.path.join(plot_dir, "umap_experiment_group.png"),
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
    # {short_marker_name: actual_adata_column_name} -- added 2026-09-10 for
    # panels (like this CRC one) whose real h5ad columns are the CODEX
    # export's full "<name> - <description> (C<n>)" strings, not the clean
    # short names used in `panel`/gating_only/gmm.per_marker_overrides/
    # cell_type_definitions.yaml. Empty by default -- identical prior
    # behavior (exact match) for any config that doesn't set it.
    column_map    = markers_cfg.get("column_map", {})

    gmm_cfg       = ct_cfg["gmm"]
    cd45_mult     = gmm_cfg["cd45_std_multiplier"]
    default_std   = gmm_cfg["default_std_multiplier"]
    per_marker    = gmm_cfg.get("per_marker_overrides", {})
    n_components  = gmm_cfg.get("n_components", 2)
    random_state  = gmm_cfg.get("random_state", 42)
    n_init        = gmm_cfg.get("n_init", 1)
    max_cells_gmm = gmm_cfg.get("max_cells_gmm", 50_000)
    # Q16 remediation (decided by Afrouz, 2026-08-18): "none" preserves the
    # exact prior behavior for any config that doesn't opt in. "yeojohnson"
    # (2026-08-19) adds a data-driven alternative to the fixed-shape
    # arcsinh/log1p transforms -- see _transform_values().
    gmm_transform = gmm_cfg.get("transform", "none")
    arcsinh_cofactor = gmm_cfg.get("arcsinh_cofactor", 5.0)

    # threshold_mode (2026-08-19): "std_multiplier" (default, unchanged
    # prior behavior) or "posterior" (GMM posterior probability of
    # positive-component membership >= confidence_level). See module
    # docstring for the statistical rationale.
    threshold_mode = gmm_cfg.get("threshold_mode", "std_multiplier")
    if threshold_mode not in ("std_multiplier", "posterior"):
        raise ValueError(
            f"Unknown gmm.threshold_mode: '{threshold_mode}'. "
            "Use 'std_multiplier' or 'posterior'."
        )
    confidence_level = gmm_cfg.get("confidence_level", 0.8)
    confidence_overrides = gmm_cfg.get("per_marker_confidence_overrides", {})
    # Per-marker transform overrides (added 2026-09-11): a marker can be
    # fit in a different transform space than gmm.transform's global
    # default -- e.g. zero-inflated markers get "arcsinh" while the rest
    # of the panel stays "none". See compute_marker_thresholds docstring.
    transform_overrides = gmm_cfg.get("per_marker_transform_overrides", {})

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
    print(f"[cell_typing] Threshold  : {threshold_mode}"
          + (f"  (confidence={confidence_level})" if threshold_mode == "posterior" else ""))
    print(f"[cell_typing] Transform  : {gmm_transform}")

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
            random_state=random_state, plot_dir=plot_dir,
            transform=gmm_transform, arcsinh_cofactor=arcsinh_cofactor,
            threshold_mode=threshold_mode, confidence_level=confidence_level,
            cd45_col=column_map.get("CD45", "CD45"),
        )

    # ── Markers for analysis (exclude gating-only markers) ───
    # Resolved through column_map so a marker like "CD45" whose real column
    # is "CD45 - hematopoietic cells (C14)" is correctly found here instead
    # of being reported missing (added 2026-09-10).
    markers_for_analysis = [m for m in panel
                             if m not in gating_only and column_map.get(m, m) in adata_cd45.var_names]
    missing = [m for m in panel
               if m not in gating_only and column_map.get(m, m) not in adata_cd45.var_names]
    if missing:
        print(f"  [cell_typing] WARNING: markers not in data: {missing}")

    # ── GMM thresholds ────────────────────────────────────────
    print(f"\n  Computing GMM thresholds (n_init={n_init}, max_cells={max_cells_gmm:,}, "
          f"transform={gmm_transform}, threshold_mode={threshold_mode})...")
    thresholds, fit_info = compute_marker_thresholds(
        adata_cd45, markers_for_analysis, std_multipliers,
        default_std=default_std, n_components=n_components, random_state=random_state,
        n_init=n_init, max_cells=max_cells_gmm,
        transform=gmm_transform, arcsinh_cofactor=arcsinh_cofactor,
        threshold_mode=threshold_mode, confidence_level=confidence_level,
        confidence_overrides=confidence_overrides,
        column_map=column_map,
        transform_overrides=transform_overrides,
    )
    threshold_report = pd.DataFrame({
        "threshold": thresholds,
        "threshold_mode": {m: threshold_mode for m in thresholds},
        # Per-marker actual transform (added 2026-09-11) -- reads
        # fit_info[m]["transform"], not the global gmm_transform, so a
        # marker with a per_marker_transform_overrides entry is reported
        # correctly instead of showing the panel-wide default.
        "transform": {m: fit_info[m]["transform"] for m in thresholds},
        "transform_lambda": {
            m: fit_info[m]["transform_params"].get("lambda") for m in thresholds
        },
    })
    threshold_report.index.name = "marker"
    threshold_report.to_csv(os.path.join(data_dir, "marker_thresholds.csv"))

    # ── Add positivity columns ────────────────────────────────
    add_positivity_columns(adata_cd45, thresholds, fit_info=fit_info,
                            transform=gmm_transform, arcsinh_cofactor=arcsinh_cofactor,
                            column_map=column_map)

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

    # -- Diagnostic plot suite (GMM per-marker distributions, positivity/
    # expression heatmaps, restyled counts chart, stacked-by-group chart,
    # UMAP by cell_type/group/core) -----------------------------------
    # Only reached by a completed automatic run or semi_automatic Phase 2
    # (Phase 1 returns early above, before cell_type exists). Local import
    # + try/except, matching preprocessing.py's summary-report hook: a bug
    # here must never fail the cell_typing step itself, since triads/
    # functional/survival depend on this function's actual output, not on
    # these plots.
    if "cell_type" in adata_cd45.obs.columns:
        try:
            generate_cell_typing_diagnostics(
                adata_cd45, thresholds, fit_info, column_map, data_dir, plot_dir,
                transform=gmm_transform, arcsinh_cofactor=arcsinh_cofactor,
                random_state=random_state,
            )
        except Exception as e:
            print(f"WARNING: cell-typing diagnostic plots failed (non-fatal -- "
                  f"cell_typing itself succeeded): {e}")


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

    # Experiment group comparison if available
    if "experiment_group" in adata.obs.columns:
        comp = (adata.obs.groupby(["experiment_group", "cell_type"])
                .size().unstack(fill_value=0))
        comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100
        comp_pct.to_csv(os.path.join(data_dir, "cell_type_by_experiment_group_pct.csv"))

        fig, ax = plt.subplots(figsize=(12, 5))
        comp_pct.T.plot(kind="bar", ax=ax, edgecolor="white")
        ax.set_title(f"Cell Type % by Experiment Group — {analysis_name}", fontsize=12)
        ax.set_ylabel("% of cells")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(title="Experiment Group", bbox_to_anchor=(1, 1))
        plt.tight_layout()
        fig.savefig(os.path.join(plot_dir, "cell_type_by_experiment_group.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)

    print(f"\n[cell_typing] Complete. Outputs in: {data_dir}")

# ══════════════════════════════════════════════════════════════════════════
# Diagnostic / reporting plots for a completed cell-typing run
# ══════════════════════════════════════════════════════════════════════════
# Moved here from the formerly-separate spatia/analysis/cell_typing_plots.py
# on 2026-09-10, at Afrouz's request, so the whole cell-typing step -- GMM
# fitting, cell-type assignment, AND its diagnostic plots -- lives in one
# file rather than being split across two. Still dataset-agnostic (no
# marker/cell-type/group names hardcoded -- everything below is driven by
# the thresholds/fit_info/column_map/adata passed in from run_cell_typing()
# above); see docs/spatia_analysis_cell_typing.md for the verification that
# established that.

# ── Small shared helpers ────────────────────────────────────────────────────

def _savefig(fig, path_no_ext: str, formats=("png", "pdf")):
    """Save one figure in every requested format. 'png' always gets dpi=150."""
    for fmt in formats:
        kwargs = {"bbox_inches": "tight"}
        if fmt == "png":
            kwargs["dpi"] = 150
        fig.savefig(f"{path_no_ext}.{fmt}", **kwargs)


def _inverse_transform_scalar_local(value_t, transform: str, arcsinh_cofactor: float,
                                     transform_params: dict):
    """
    Local, minimal re-implementation of cell_typing.py's
    _inverse_transform_scalar -- duplicated (not imported) to keep this
    plotting module decoupled from cell_typing.py's private internals and
    avoid a circular import (cell_typing.py imports FROM this module).
    Covers "none"/"arcsinh"/"log1p" exactly. "yeojohnson" is intentionally
    NOT reproduced here (it needs the fitted lambda's full inverse formula)
    -- not needed by plot_marker_distributions() any more (2026-09-11): the
    transformed-space panel now uses the FORWARD transform (see
    _yeojohnson_forward_local below), which has a closed form for every
    supported transform including yeojohnson, so no inverse is required.
    """
    if transform == "none":
        return float(value_t)
    elif transform == "arcsinh":
        return float(np.sinh(value_t) * arcsinh_cofactor)
    elif transform == "log1p":
        return float(np.expm1(value_t))
    else:
        raise ValueError(f"transform '{transform}' not supported for local inverse")


def _yeojohnson_forward_local(x: np.ndarray, lmbda: float) -> np.ndarray:
    """
    Forward Yeo-Johnson transform, vectorized. Local re-implementation
    (not imported, same decoupling reasoning as _inverse_transform_scalar_local
    above) -- standard closed form, e.g. Yeo & Johnson (2000):
        x >= 0, lmbda != 0 :  ((x + 1)**lmbda - 1) / lmbda
        x >= 0, lmbda == 0 :  log(x + 1)
        x <  0, lmbda != 2 :  -((-x + 1)**(2 - lmbda) - 1) / (2 - lmbda)
        x <  0, lmbda == 2 :  -log(-x + 1)
    Added 2026-09-11 so plot_marker_distributions() can draw a correct
    transformed-space panel (and native-space GMM component curves, no
    Jacobian correction needed since the mixture is drawn directly in the
    space it was fit in) for yeojohnson markers too -- previously that
    transform fell back to histogram + threshold line only, with no
    curves, because only the (harder) inverse was implemented and only
    for arcsinh/log1p/none.
    """
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    pos = x >= 0
    if lmbda != 0:
        out[pos] = ((x[pos] + 1) ** lmbda - 1) / lmbda
    else:
        out[pos] = np.log1p(x[pos])
    neg = ~pos
    if lmbda != 2:
        out[neg] = -((-x[neg] + 1) ** (2 - lmbda) - 1) / (2 - lmbda)
    else:
        out[neg] = -np.log1p(-x[neg])
    return out


def _get_column(adata, marker: str, column_map: dict):
    real_col = column_map.get(marker, marker)
    if real_col not in adata.var_names:
        return None, real_col
    vals = adata[:, real_col].X
    if hasattr(vals, "toarray"):
        vals = vals.toarray().flatten()
    else:
        vals = np.asarray(vals).flatten()
    return vals, real_col


# ── 1. Per-marker GMM distribution plots ────────────────────────────────────

def plot_marker_distributions(adata, thresholds: dict, fit_info: dict, column_map: dict,
                               plot_dir: str, transform: str = "none",
                               arcsinh_cofactor: float = 5.0, formats=("png", "pdf")):
    """
    Three panels per marker, left to right:
      1. Raw expression -- histogram + threshold line (raw-scale, always
         meaningful regardless of transform).
      2. This marker's own transformed space -- histogram of the actual
         transformed values PLUS the fitted GMM component curves drawn
         natively in that space (no Jacobian correction needed, since this
         is exactly the space compute_marker_thresholds() fit the GMM in --
         more correct than overlaying curves on the raw axis, which the
         single-panel version before 2026-09-11 did approximately). Reads
         fit_info[marker]["transform"] -- the transform THIS marker
         actually used (per_marker_transform_overrides-aware), not the
         `transform` argument, which is only a fallback for markers with
         no fit_info entry. If that marker's transform is "none", this
         panel is identical in shape to panel 1 (labeled as such) rather
         than omitted, so every marker gets the same 3-panel layout.
      3. Positive vs. negative split -- the raw-expression histogram
         colored by this marker's actual `{marker}_pos` call, so the plot
         shows what the threshold decision actually did to the real cell
         population, not just where the line sits.

    Panel 2's curves REUSE fit_info[marker]["gmm"] -- the exact fit that
    produced the real threshold, not a fresh re-fit -- and respect this
    panel's documented zero-inflation handling (14/53 CRC markers have
    30-50% exact-zero/point-mass intensity -- see cell_typing.py's
    _fit_gmm / _detect_point_mass docstrings), since fit_info's
    "clean"/"gmm" already reflect that exclusion. yeojohnson markers now
    get real component curves too (2026-09-11, via _yeojohnson_forward_local)
    -- previously skipped because only the harder inverse transform was
    implemented, and inverting isn't needed when panel 2 draws natively in
    transformed space instead of back on the raw axis.
    """
    column_map = column_map or {}
    fit_info = fit_info or {}
    sub_dir = os.path.join(plot_dir, "marker_distributions")
    os.makedirs(sub_dir, exist_ok=True)

    NEG_COLOR, POS_COLOR = "#4C72B0", "#DD8452"  # colorblind-safe blue/orange pair

    n_done, skipped = 0, []
    for marker, threshold in thresholds.items():
        info = fit_info.get(marker)
        vals, real_col = _get_column(adata, marker, column_map)
        if vals is None or info is None or info.get("gmm") is None:
            skipped.append(marker)
            continue

        finite_mask = np.isfinite(vals)
        clean = vals[finite_mask]
        if len(clean) < 10:
            skipped.append(marker)
            continue

        marker_transform = info.get("transform", transform)
        gmm = info["gmm"]
        low_idx, high_idx = info["low_idx"], info["high_idx"]
        means_t = gmm.means_.flatten()
        stds_t = np.sqrt(gmm.covariances_.flatten())
        weights = gmm.weights_
        tparams = info.get("transform_params", {})

        pos_col = f"{marker}_pos"
        pos_bool = (adata.obs[pos_col].to_numpy()[finite_mask]
                    if pos_col in adata.obs.columns else None)

        try:
            lo, hi = np.percentile(clean, [0.1, 99.9])
            if lo >= hi:
                lo, hi = float(np.min(clean)), float(np.max(clean))
            x_raw = np.linspace(lo, hi, 500)

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(19, 5.5))

            # ── Panel 1: raw expression ──────────────────────────────
            ax1.hist(clean, bins=60, alpha=0.7, density=True, color="#8c8c8c",
                     edgecolor="white", linewidth=0.3, label=f"All cells (n={len(clean):,})")
            ax1.axvline(threshold, color="#C44E52", linestyle="--", linewidth=1.8,
                        label=f"Threshold: {threshold:.3f}")
            ax1.set_title("Raw expression", fontsize=11)
            ax1.set_xlabel("Expression level (raw)")
            ax1.set_ylabel("Density")
            ax1.legend(fontsize="small", frameon=False)

            # ── Panel 2: this marker's own transformed space ────────
            if marker_transform == "none":
                ax2.hist(clean, bins=60, alpha=0.7, density=True, color="#8c8c8c",
                         edgecolor="white", linewidth=0.3)
                ax2.axvline(threshold, color="#C44E52", linestyle="--", linewidth=1.8)
                ax2.set_title("Transformed space (transform: none)", fontsize=11)
                ax2.set_xlabel("Expression level (raw)")
            else:
                if marker_transform == "arcsinh":
                    clean_t = np.arcsinh(clean / arcsinh_cofactor)
                    x_t = np.arcsinh(x_raw / arcsinh_cofactor)
                elif marker_transform == "log1p":
                    clean_t = np.log1p(clean)
                    x_t = np.log1p(x_raw)
                elif marker_transform == "yeojohnson":
                    lam = tparams.get("lambda", 1.0)
                    clean_t = _yeojohnson_forward_local(clean, lam)
                    x_t = _yeojohnson_forward_local(x_raw, lam)
                else:
                    clean_t, x_t = clean, x_raw

                ax2.hist(clean_t, bins=60, alpha=0.7, density=True, color="#8c8c8c",
                         edgecolor="white", linewidth=0.3, label=f"All cells (n={len(clean_t):,})")
                comp_names = {low_idx: "negative", high_idx: "positive"}
                mixture = np.zeros_like(x_t)
                for i in (low_idx, high_idx):
                    comp_pdf = weights[i] * norm.pdf(x_t, means_t[i], stds_t[i])
                    mixture += comp_pdf
                    ax2.plot(x_t, comp_pdf, linewidth=2,
                             label=f"{comp_names[i]} (μ={means_t[i]:.2f}, w={weights[i]:.2f})")
                ax2.plot(x_t, mixture, "k--", linewidth=1.3, label="GMM mixture")
                # Threshold in this same transformed space: forward-transform
                # the already-known raw threshold (exact, since it was
                # produced by inverse-transforming this same fit).
                if marker_transform == "arcsinh":
                    threshold_t = np.arcsinh(threshold / arcsinh_cofactor)
                elif marker_transform == "log1p":
                    threshold_t = np.log1p(threshold)
                else:
                    threshold_t = _yeojohnson_forward_local(np.array([threshold]), tparams.get("lambda", 1.0))[0]
                ax2.axvline(threshold_t, color="#C44E52", linestyle="--", linewidth=1.8,
                            label=f"Threshold: {threshold_t:.3f}")
                ax2.set_title(f"Transformed space (transform: {marker_transform})", fontsize=11)
                ax2.set_xlabel("Expression level (transformed)")
                ax2.legend(fontsize="small", frameon=False)
            ax2.set_ylabel("Density")

            # ── Panel 3: positive vs. negative split ─────────────────
            if pos_bool is not None and pos_bool.any() and (~pos_bool).any():
                ax3.hist([clean[~pos_bool], clean[pos_bool]], bins=60, stacked=True,
                         density=False, color=[NEG_COLOR, POS_COLOR],
                         label=[f"negative (n={int((~pos_bool).sum()):,})",
                                f"positive (n={int(pos_bool.sum()):,})"],
                         edgecolor="white", linewidth=0.2)
                ax3.legend(fontsize="small", frameon=False)
            else:
                ax3.hist(clean, bins=60, color="#8c8c8c", edgecolor="white", linewidth=0.3)
            ax3.axvline(threshold, color="#C44E52", linestyle="--", linewidth=1.8)
            ax3.set_title("Positive / negative split (raw scale)", fontsize=11)
            ax3.set_xlabel("Expression level (raw)")
            ax3.set_ylabel("Cells")

            if pos_col in adata.obs.columns:
                pct_pos = 100 * adata.obs[pos_col].mean()
                n_pos = int(adata.obs[pos_col].sum())
                suptitle = f"{marker}  ({real_col})\nPositive: {n_pos:,} cells ({pct_pos:.1f}%)"
            else:
                suptitle = f"{marker}  ({real_col})"
            fig.suptitle(suptitle, fontsize=12)
            for ax in (ax1, ax2, ax3):
                ax.spines[["top", "right"]].set_visible(False)
            plt.tight_layout(rect=[0, 0, 1, 0.93])

            safe_marker = marker.replace("/", "_").replace(" ", "_")
            _savefig(fig, os.path.join(sub_dir, f"{safe_marker}_distribution"), formats)
            plt.close(fig)
            n_done += 1
        except Exception as e:
            plt.close("all")
            skipped.append(f"{marker} ({e})")

    print(f"  [diagnostics] Marker distribution plots: {n_done} written to {sub_dir}/"
          + (f"  -- skipped: {skipped}" if skipped else ""))


# ── 2. Marker-positivity-% heatmap by cell type ─────────────────────────────

def plot_positivity_heatmap(adata, markers: list, plot_dir: str, data_dir: str,
                             formats=("png", "pdf")):
    """
    Rows = cell type, columns = marker, value = % of that cell type's cells
    positive for that marker (from the existing {marker}_pos columns --
    same convention plot_marker_distributions() reads). Sequential single-hue
    colormap (0-100%), per the project's dataviz convention for magnitude.
    """
    if "cell_type" not in adata.obs.columns:
        print("  [diagnostics] Skipping positivity heatmap -- no 'cell_type' column.")
        return

    cell_types = adata.obs["cell_type"].value_counts().index.tolist()
    pos_cols = [m for m in markers if f"{m}_pos" in adata.obs.columns]
    if not pos_cols:
        print("  [diagnostics] Skipping positivity heatmap -- no *_pos columns found.")
        return

    table = pd.DataFrame(index=cell_types, columns=pos_cols, dtype=float)
    for ct in cell_types:
        sub = adata.obs.loc[adata.obs["cell_type"] == ct, [f"{m}_pos" for m in pos_cols]]
        table.loc[ct] = 100 * sub.mean().values

    table.to_csv(os.path.join(data_dir, "cell_type_marker_positivity_pct.csv"))

    fig_w = max(10, 0.35 * len(pos_cols) + 3)
    fig_h = max(6, 0.35 * len(cell_types) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(table, cmap="rocket_r", vmin=0, vmax=100, cbar_kws={"label": "% positive"}, ax=ax)
    ax.set_title("Marker Positivity (%) by Cell Type", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, "cell_type_marker_positivity_heatmap"), formats)
    plt.close(fig)
    print(f"  [diagnostics] Positivity heatmap: {len(cell_types)} cell types x {len(pos_cols)} markers")


# ── 3. Mean-expression heatmap by cell type ─────────────────────────────────

def plot_mean_expression_heatmap(adata, markers: list, column_map: dict, plot_dir: str,
                                  data_dir: str, formats=("png", "pdf")):
    """Rows = cell type, columns = marker, value = mean raw expression."""
    if "cell_type" not in adata.obs.columns:
        return
    column_map = column_map or {}
    cell_types = adata.obs["cell_type"].value_counts().index.tolist()
    usable = [(m, column_map.get(m, m)) for m in markers if column_map.get(m, m) in adata.var_names]
    if not usable:
        print("  [diagnostics] Skipping mean-expression heatmap -- no usable markers.")
        return

    table = pd.DataFrame(index=cell_types, columns=[m for m, _ in usable], dtype=float)
    for ct in cell_types:
        mask = (adata.obs["cell_type"] == ct).values
        sub = adata[mask]
        for m, real_col in usable:
            vals = sub[:, real_col].X
            vals = vals.toarray().flatten() if hasattr(vals, "toarray") else np.asarray(vals).flatten()
            table.loc[ct, m] = float(np.mean(vals)) if len(vals) else np.nan

    table.to_csv(os.path.join(data_dir, "cell_type_marker_means.csv"))

    fig_w = max(10, 0.35 * len(usable) + 3)
    fig_h = max(6, 0.35 * len(cell_types) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(table, cmap="mako_r", cbar_kws={"label": "Mean expression (raw)"}, ax=ax,
                center=0)
    ax.set_title("Mean Marker Expression by Cell Type", fontsize=13)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, "cell_type_marker_expression_heatmap"), formats)
    plt.close(fig)
    print(f"  [diagnostics] Mean-expression heatmap: {len(cell_types)} cell types x {len(usable)} markers")


# ── 4. Restyled cell-type-counts bar chart ──────────────────────────────────

def plot_cell_type_counts_pretty(adata, plot_dir: str, formats=("png", "pdf")):
    """
    Replaces automatic mode's original cell_type_counts.png/.pdf (written
    moments earlier by _save_and_plot()) with a horizontal, sorted,
    value-labeled version -- the vertical/rotated-label version was flagged
    as unreadable with this many cell-type categories (long names collide
    at 45deg rotation). Horizontal bars sidestep that entirely.
    """
    if "cell_type" not in adata.obs.columns:
        return
    counts = adata.obs["cell_type"].value_counts().sort_values(ascending=True)
    total = int(counts.sum())
    pct = counts / total * 100

    fig_h = max(5, 0.32 * len(counts) + 1.5)
    fig, ax = plt.subplots(figsize=(9.5, fig_h))
    colors = ["#9a9a9a" if ct == "Unassigned" else "#4C72B0" for ct in counts.index]
    bars = ax.barh(counts.index.astype(str), counts.values, color=colors,
                    edgecolor="white", height=0.72, zorder=3)
    for bar, c, p in zip(bars, counts.values, pct.values):
        ax.text(bar.get_width() + total * 0.006, bar.get_y() + bar.get_height() / 2,
                f"{c:,}  ({p:.1f}%)", va="center", fontsize=8.5, color="#333333")

    ax.set_xlabel("Number of cells")
    ax.set_title(f"Cell Type Distribution  —  {total:,} cells total", fontsize=13, pad=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.grid(True, color="#e3e3e3", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.margins(x=0.14)
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, "cell_type_counts"), formats)
    plt.close(fig)
    print(f"  [diagnostics] Restyled cell_type_counts.png ({len(counts)} cell types)")


# ── 5. Stacked cell-type-by-group bar chart ─────────────────────────────────

def plot_stacked_by_group(adata, group_col: str, plot_dir: str, data_dir: str,
                           formats=("png", "pdf")):
    """
    100%-stacked bar chart: one bar per experiment_group (e.g. CLR, DII),
    each fully stacked into its cell-type composition -- segments sum to
    exactly 100% per bar. This is the standard reading of "stacked barplot
    of cell types" and is DIFFERENT from a first draft of this function
    that put cell_type on the x-axis and stacked CLR+DII on top of each
    other per type -- caught in verification (2026-09-10, synthetic test):
    that version's bars summed to (%-within-CLR + %-within-DII), which
    isn't bounded or meaningful (e.g. a bar could total ~96%, another
    ~28%), since each group's percentages are independently normalized to
    100% across cell types, not across groups. Group-on-x-axis, cell-type
    stacked is the only arrangement where the stack total is actually 100%.

    Adds this alongside (not replacing -- Afrouz confirmed the existing
    grouped/dodged version is fine as-is) _save_and_plot()'s existing
    cell_type_by_experiment_group.png, which shows the same per-group
    percentages as side-by-side bars instead of stacked ones.
    """
    if "cell_type" not in adata.obs.columns or group_col not in adata.obs.columns:
        return
    comp = pd.crosstab(adata.obs[group_col], adata.obs["cell_type"], normalize="index") * 100
    # Canonical, reproducible column (segment) order: most-common cell type
    # first, so colors/order are stable across re-runs of this cohort.
    order = adata.obs["cell_type"].value_counts().index.tolist()
    comp = comp[[c for c in order if c in comp.columns]]
    comp.to_csv(os.path.join(data_dir, f"cell_type_by_{group_col}_stacked100_pct.csv"))

    n_types = comp.shape[1]
    palette = sns.color_palette("husl", n_types)
    if "Unassigned" in comp.columns:
        palette[list(comp.columns).index("Unassigned")] = (0.6, 0.6, 0.6)

    fig, ax = plt.subplots(figsize=(max(6, 1.6 * comp.shape[0] + 2), 7))
    comp.plot(kind="bar", stacked=True, ax=ax, color=palette, edgecolor="white",
              linewidth=0.5, width=0.55)
    ax.set_title(f"Cell Type Composition by {group_col} (100% stacked)", fontsize=13)
    ax.set_ylabel("% of cells")
    ax.set_xlabel("")
    ax.set_ylim(0, 100)
    plt.setp(ax.get_xticklabels(), rotation=0)
    ax.legend(title="Cell type", bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=7, ncol=1 if n_types <= 15 else 2)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, f"cell_type_by_{group_col}_stacked100"), formats)
    plt.close(fig)
    print(f"  [diagnostics] 100%-stacked cell-type-by-{group_col} chart written ({n_types} cell types)")


# ── 6. UMAP: cell type / group / core ───────────────────────────────────────

def compute_and_plot_umaps(adata, plot_dir: str, data_dir: str, group_col: str, tissue_col: str,
                            random_state: int = 42, formats=("png", "pdf")):
    """
    Computes PCA -> neighbors -> UMAP if not already present, then saves
    three colorings: cell_type (the "whole-TMA" map), group_col, tissue_col
    (one combined UMAP colored by core -- Afrouz's 2026-09-10 choice over a
    per-core facet grid). Also exports umap_coordinates.csv (UMAP1/UMAP2 +
    cell_type/group/tissue per cell) since the h5ad _save_and_plot() already
    wrote does NOT include these embeddings (they're computed after that
    write) -- re-writing the full 300k x 80 h5ad just to add 2 columns
    wasn't worth the extra I/O; the CSV is the cheap way to keep this
    reusable outside this script.

    IMPORTANT, and different from run_clustering() (semi_automatic mode):
    this does NOT call sc.pp.normalize_total / sc.pp.log1p / sc.pp.scale
    before neighbors. This pipeline's preprocessing already per-cell
    z-scores every marker (confirmed: real thresholds on this exact
    dataset are legitimately negative, e.g. CD7=-0.0039, CD3=-0.0393,
    CD163=-0.0552 -- log1p on any z-scored value below -1 produces
    NaN/-inf and would silently corrupt the embedding). NOTE:
    run_clustering() (used by semi_automatic mode) still has this exact
    log1p-on-already-z-scored-data issue -- not fixed here since automatic
    mode never calls it, but worth knowing before switching this config to
    semi_automatic mode.

    No PCA step either -- matches Afrouz's own automatic cell-typing
    reference notebook (05_1_a, LILRB2 mouse study), which calls
    sc.pp.neighbors directly on the marker matrix with no PCA reduction
    first (2026-09-10 change; the panel here is already low-dimensional --
    dozens of markers, not thousands of genes -- so PCA isn't needed the
    way it is for scRNA-seq). No cell subsampling either, matching
    Afrouz's explicit choice to run on the full pooled cohort.

    use_rep="X" is REQUIRED below, not cosmetic (bug found 2026-09-10 on
    the first real 80-marker/300k-cell run): scanpy's sc.pp.neighbors
    silently ignores "no PCA" intent once .X has more dimensions than its
    own internal default PC count (50) -- without use_rep="X" it prints
    "Falling back to preprocessing with sc.pp.pca and default params" and
    runs PCA anyway. The earlier synthetic test (6 markers) never
    triggered this, since 6 < 50. Confirmed via the real run's own log.
    """
    if "X_umap" not in adata.obsm:
        print(f"  [diagnostics] Computing neighbors -> UMAP directly on the marker matrix "
              f"({adata.n_vars} markers, {adata.n_obs:,} cells, no PCA/subsampling -- "
              f"matches the reference notebook's approach)...")
        sc.pp.neighbors(adata, n_neighbors=15, random_state=random_state, use_rep="X")
        sc.tl.umap(adata, random_state=random_state)
    else:
        print("  [diagnostics] X_umap already present (semi_automatic clustering ran earlier) -- reusing it.")

    coords = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs.index)
    for col in ("cell_type", group_col, tissue_col):
        if col in adata.obs.columns:
            coords[col] = adata.obs[col].values
    coords.to_csv(os.path.join(data_dir, "umap_coordinates.csv"))

    colorings = [("cell_type", "UMAP — Cell Types (whole TMA)", "umap_cell_types_all")]
    if group_col in adata.obs.columns:
        colorings.append((group_col, f"UMAP — {group_col}", f"umap_by_{group_col}"))
    if tissue_col in adata.obs.columns:
        n_cores = adata.obs[tissue_col].nunique()
        colorings.append((tissue_col, f"UMAP — {tissue_col} ({n_cores} cores)", f"umap_by_{tissue_col}"))

    for color_key, title, fname in colorings:
        try:
            fig, ax = plt.subplots(figsize=(9, 8))
            legend_loc = "right margin" if adata.obs[color_key].nunique() <= 30 else None
            sc.pl.umap(adata, color=color_key, ax=ax, show=False, title=title,
                       legend_loc=legend_loc, legend_fontsize=7)
            _savefig(fig, os.path.join(plot_dir, fname), formats)
            plt.close(fig)
        except Exception as e:
            plt.close("all")
            print(f"  [diagnostics] UMAP coloring by '{color_key}' failed (non-fatal): {e}")

    print(f"  [diagnostics] UMAP: {len(colorings)} coloring(s) written, coordinates saved to "
          f"{data_dir}/umap_coordinates.csv")


# ── 7. Dotplot: marker expression across cell types ─────────────────────────
# Added 2026-09-10, at Afrouz's request to mirror her reference notebook's
# outputs one-for-one. sc.pl.dotplot with dendrogram=True clusters cell
# types by expression similarity and standard_scale="var" puts every marker
# on the same 0-1 color scale, exactly matching the notebook's call.

def plot_dotplot(adata, thresholds: dict, column_map: dict, plot_dir: str, formats=("png", "pdf")):
    if "cell_type" not in adata.obs.columns or not HAS_SCANPY:
        return
    column_map = column_map or {}
    var_names = [column_map.get(m, m) for m in thresholds if column_map.get(m, m) in adata.var_names]
    n_types = adata.obs["cell_type"].nunique()
    if len(var_names) < 2 or n_types < 2:
        print(f"  [diagnostics] Skipping dotplot -- need >=2 markers and >=2 cell types "
              f"(have {len(var_names)} markers, {n_types} cell types).")
        return
    try:
        plt.figure(figsize=(max(10, 0.4 * len(var_names) + 4), max(6, 0.35 * n_types + 3)))
        sc.pl.dotplot(adata, var_names=var_names, groupby="cell_type",
                      standard_scale="var", dendrogram=True, color_map="viridis", show=False)
        fig = plt.gcf()
        fig.suptitle("Marker Expression Across Cell Types", fontsize=14, y=1.02)
        plt.tight_layout()
        _savefig(fig, os.path.join(plot_dir, "cell_type_marker_dotplot"), formats)
        plt.close(fig)
        print(f"  [diagnostics] Dotplot: {len(var_names)} markers x {n_types} cell types")
    except Exception as e:
        plt.close("all")
        raise e


# ── 8. Per-marker UMAP: expression + positivity, side by side ───────────────
# Added 2026-09-10. Matches the notebook's per-marker UMAP panels, one file
# per marker (raw expression on the left, binary {marker}_pos on the
# right). Requires compute_and_plot_umaps() to have already run (reuses
# its X_umap embedding -- never recomputes it).

def plot_marker_umaps(adata, thresholds: dict, column_map: dict, plot_dir: str,
                       fit_info: Optional[dict] = None, arcsinh_cofactor: float = 5.0,
                       formats=("png", "pdf")):
    """
    Three panels per marker, left to right, all coloring the SAME fixed 2D
    UMAP embedding (computed once on the untransformed marker matrix, per
    this project's no-PCA convention -- only the point colors differ
    across panels, never the embedding):
      1. Raw expression.
      2. This marker's own transformed-space expression -- reads
         fit_info[marker]["transform"] (per_marker_transform_overrides-aware,
         added 2026-09-11), i.e. the same space actually used to threshold
         this marker. Falls back to "none" (identical to panel 1, just a
         different colormap so it's visually distinguishable) if fit_info
         has no entry for this marker.
      3. Positive/negative call (unchanged from before 2026-09-11).
    """
    if "X_umap" not in adata.obsm or not HAS_SCANPY:
        print("  [diagnostics] Skipping per-marker UMAPs -- no UMAP embedding yet "
              "(compute_and_plot_umaps() must run first).")
        return
    column_map = column_map or {}
    fit_info = fit_info or {}
    sub_dir = os.path.join(plot_dir, "marker_umaps")
    os.makedirs(sub_dir, exist_ok=True)

    n_done, skipped = 0, []
    for marker, threshold in thresholds.items():
        real_col = column_map.get(marker, marker)
        pos_col = f"{marker}_pos"
        if real_col not in adata.var_names or pos_col not in adata.obs.columns:
            skipped.append(marker)
            continue
        tmp_col = f"__{marker}_transformed_tmp"
        try:
            info = fit_info.get(marker, {})
            marker_transform = info.get("transform", "none")

            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
            sc.pl.umap(adata, color=real_col, cmap="plasma", ax=ax1, show=False)
            ax1.set_title(f"{marker} expression (raw)")

            if marker_transform == "none":
                sc.pl.umap(adata, color=real_col, cmap="viridis", ax=ax2, show=False)
                ax2.set_title(f"{marker} expression (transform: none)")
            else:
                raw_vals, _ = _get_column(adata, marker, column_map)
                if marker_transform == "arcsinh":
                    t_vals = np.arcsinh(raw_vals / arcsinh_cofactor)
                elif marker_transform == "log1p":
                    t_vals = np.log1p(raw_vals)
                elif marker_transform == "yeojohnson":
                    lam = info.get("transform_params", {}).get("lambda", 1.0)
                    t_vals = _yeojohnson_forward_local(raw_vals, lam)
                else:
                    t_vals = raw_vals
                adata.obs[tmp_col] = t_vals
                sc.pl.umap(adata, color=tmp_col, cmap="viridis", ax=ax2, show=False)
                ax2.set_title(f"{marker} expression (transform: {marker_transform})")

            sc.pl.umap(adata, color=pos_col, cmap="coolwarm", ax=ax3, show=False)
            ax3.set_title(f"{marker}+ cells (threshold: {threshold:.3f})")
            plt.tight_layout()
            safe_marker = marker.replace("/", "_").replace(" ", "_")
            _savefig(fig, os.path.join(sub_dir, f"{safe_marker}_umap"), formats)
            plt.close(fig)
            n_done += 1
        except Exception as e:
            plt.close("all")
            skipped.append(f"{marker} ({e})")
        finally:
            # Never leave the scratch column behind, success or failure.
            if tmp_col in adata.obs.columns:
                del adata.obs[tmp_col]

    print(f"  [diagnostics] Per-marker UMAPs: {n_done} written to {sub_dir}/"
          + (f"  -- skipped: {skipped}" if skipped else ""))


# ── 9. Group statistics: chi-square test + significant-difference plots ─────
# Added 2026-09-10. The reference notebook's chi2 test only ever handled
# exactly 2 conditions (a fixed 2x2 contingency table). Generalized here to
# any number of groups: for each cell type, a [this type vs everything
# else] x [group1..groupN] table. The signed diverging "difference" bar
# chart is inherently pairwise, so it's only produced when there are
# exactly 2 groups (matching the notebook exactly in that specific case).

def plot_group_statistics(adata, group_col: str, plot_dir: str, data_dir: str, formats=("png", "pdf")):
    if "cell_type" not in adata.obs.columns or group_col not in adata.obs.columns:
        return
    from scipy.stats import chi2_contingency

    groups = sorted(adata.obs[group_col].dropna().unique().tolist())
    if len(groups) < 2:
        print(f"  [diagnostics] Skipping group statistics -- only {len(groups)} group(s) in '{group_col}'.")
        return

    cell_types = adata.obs["cell_type"].value_counts().index.tolist()
    pct = (pd.crosstab(adata.obs["cell_type"], adata.obs[group_col], normalize="columns") * 100)
    pct = pct.reindex(index=cell_types, columns=groups)

    results = {}
    for ct in cell_types:
        is_type = (adata.obs["cell_type"] == ct)
        table = np.array([
            [int((is_type & (adata.obs[group_col] == g)).sum()) for g in groups],
            [int((~is_type & (adata.obs[group_col] == g)).sum()) for g in groups],
        ])
        try:
            chi2, p, _, _ = chi2_contingency(table)
        except ValueError:
            chi2, p = np.nan, np.nan
        results[ct] = {"chi2": chi2, "p_value": p,
                       "significant": bool(p < 0.05) if pd.notna(p) else False}

    chi2_df = pd.DataFrame(results).T
    chi2_df.to_csv(os.path.join(data_dir, f"cell_type_{group_col}_chi2_statistics.csv"))

    sig = chi2_df[chi2_df["significant"]]
    if not sig.empty:
        sig_pct = pct.loc[sig.index]
        fig, ax = plt.subplots(figsize=(11, max(4, 0.5 * len(sig) + 2)))
        sig_pct.plot(kind="barh", ax=ax)
        ax.set_title(f"Cell Types Significantly Different Across {group_col} (p<0.05)", fontsize=12)
        ax.set_xlabel("% of cells")
        ax.legend(title=group_col, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        _savefig(fig, os.path.join(plot_dir, f"cell_type_{group_col}_significant_differences"), formats)
        plt.close(fig)

    if len(groups) == 2:
        diff = (pct[groups[1]] - pct[groups[0]]).sort_values()
        fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(diff) + 2)))
        colors = ["indianred" if v > 0 else "steelblue" for v in diff.values]
        ax.barh(diff.index.astype(str), diff.values, color=colors)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(f"Cell Type % Difference: {groups[1]} vs {groups[0]}", fontsize=12)
        ax.set_xlabel("Percentage-point difference")
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        _savefig(fig, os.path.join(plot_dir, f"cell_type_{group_col}_difference"), formats)
        plt.close(fig)

    n_sig = int(chi2_df["significant"].sum())
    print(f"  [diagnostics] Group statistics ({group_col}, {len(groups)} groups): "
          f"{n_sig}/{len(cell_types)} cell types significant (p<0.05) -- "
          f"{data_dir}/cell_type_{group_col}_chi2_statistics.csv")


# ── 10. Unassigned-cell diagnostics ──────────────────────────────────────────
# Added 2026-09-10. This pipeline's rule engine (assign_cell_types_automatic)
# labels unmatched cells "Unassigned" -- the equivalent of the reference
# notebook's "Unknown" catch-all. Reproduces the notebook's per-marker
# known-vs-unknown histogram grid and impact-ranked threshold-adjustment
# recommendations. Skipped entirely if there are no Unassigned cells.

def plot_unassigned_diagnostics(adata, thresholds: dict, column_map: dict, plot_dir: str,
                                 data_dir: str, formats=("png", "pdf")):
    if "cell_type" not in adata.obs.columns:
        return
    unassigned_mask = (adata.obs["cell_type"] == "Unassigned").values
    n_unassigned = int(unassigned_mask.sum())
    if n_unassigned == 0:
        print("  [diagnostics] No 'Unassigned' cells -- skipping unassigned-cell diagnostics.")
        return
    column_map = column_map or {}
    known_mask = ~unassigned_mask

    pos_cols = [m for m in thresholds if f"{m}_pos" in adata.obs.columns]
    if pos_cols:
        rows = []
        for m in pos_cols:
            col = f"{m}_pos"
            rows.append({
                "marker": m,
                "Unassigned": 100 * adata.obs.loc[unassigned_mask, col].mean(),
                "Known": 100 * adata.obs.loc[known_mask, col].mean(),
            })
        comp = pd.DataFrame(rows).set_index("marker")
        comp.to_csv(os.path.join(data_dir, "unassigned_vs_known_positivity_pct.csv"))
        fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(pos_cols)), 6))
        comp.plot(kind="bar", ax=ax, color=["#C44E52", "#4C72B0"])
        ax.set_title(f"Marker Positivity: Unassigned ({n_unassigned:,} cells) vs Known", fontsize=12)
        ax.set_ylabel("% positive")
        ax.tick_params(axis="x", rotation=90)
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        _savefig(fig, os.path.join(plot_dir, "unassigned_vs_known_positivity"), formats)
        plt.close(fig)

    markers = list(thresholds.keys())
    n_cols = 3
    n_rows = max(1, -(-len(markers) // n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows))
    axes = np.atleast_1d(axes).flatten()
    recommendations = {}

    for i, marker in enumerate(markers):
        vals, real_col = _get_column(adata, marker, column_map)
        ax = axes[i]
        if vals is None:
            ax.axis("off")
            continue
        unk_expr = vals[unassigned_mask]
        kn_expr = vals[known_mask]
        threshold = thresholds[marker]
        unk_mean = float(np.mean(unk_expr)) if len(unk_expr) else float("nan")
        kn_mean = float(np.mean(kn_expr)) if len(kn_expr) else float("nan")
        ax.hist(kn_expr, bins=40, alpha=0.5, color="#4C72B0", label=f"Known (μ={kn_mean:.2f})")
        ax.hist(unk_expr, bins=40, alpha=0.5, color="#C44E52", label=f"Unassigned (μ={unk_mean:.2f})")
        ax.axvline(threshold, color="black", linestyle="--", linewidth=1.2,
                    label=f"Threshold: {threshold:.2f}")
        pct_unk_pos = 100 * float(np.mean(unk_expr > threshold)) if len(unk_expr) else 0.0
        if unk_mean > kn_mean:
            suggested = float(np.percentile(unk_expr, 90)) if len(unk_expr) else threshold
            direction = "up" if suggested > threshold else "down"
        else:
            suggested, direction = threshold, "none"
        if direction != "none":
            ax.axvline(suggested, color="#C44E52", linestyle=":", linewidth=1.2,
                        label=f"Suggested: {suggested:.2f}")
        recommendations[marker] = {
            "current_threshold": threshold, "unassigned_mean": unk_mean, "known_mean": kn_mean,
            "unassigned_pct_positive": pct_unk_pos, "suggested_threshold": suggested,
            "direction": direction, "impact_score": pct_unk_pos * abs(suggested - threshold),
        }
        ax.set_title(marker, fontsize=9)
        ax.legend(fontsize=6, loc="upper right")
    for j in range(len(markers), len(axes)):
        axes[j].axis("off")
    fig.suptitle(f"Marker Expression: Unassigned ({n_unassigned:,}) vs Known Cells", fontsize=14, y=1.01)
    plt.tight_layout()
    _savefig(fig, os.path.join(plot_dir, "unassigned_marker_analysis_grid"), formats)
    plt.close(fig)

    rec_df = pd.DataFrame(recommendations).T.sort_values("impact_score", ascending=False)
    rec_df.to_csv(os.path.join(data_dir, "unassigned_threshold_recommendations.csv"))

    top = rec_df.head(min(10, len(rec_df)))
    if not top.empty and top["impact_score"].max() > 0:
        fig, ax = plt.subplots(figsize=(10, max(3, 0.4 * len(top) + 1.5)))
        ax.barh(top.index.astype(str)[::-1], top["impact_score"].values[::-1], color="#DD8452")
        ax.set_title("Top Markers by Potential Impact on Unassigned-Cell Reclassification", fontsize=11)
        ax.set_xlabel("Impact score (Unassigned+% x |suggested threshold shift|)")
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        _savefig(fig, os.path.join(plot_dir, "unassigned_threshold_recommendations_top10"), formats)
        plt.close(fig)

    print(f"  [diagnostics] Unassigned-cell diagnostics: {n_unassigned:,} / {len(adata):,} cells "
          f"({100 * n_unassigned / len(adata):.1f}%) -- recommendations in "
          f"{data_dir}/unassigned_threshold_recommendations.csv")


# ── 11. Report bundle: HTML + PDF + text summary ─────────────────────────────
# Added 2026-09-10, to mirror the reference notebook's report/ output.
# Composed from the PNGs this module already wrote to plot_dir (rather than
# a third duplicate of every plotting call, which is how the notebook does
# it -- once for the main script, once as a base64 PNG for its Jinja2 HTML,
# once again inside PdfPages -- three copies of the same plotting logic
# that can silently drift apart). Built with plain string formatting, not
# jinja2/markdown: this must also run unattended on JupyterHub/HPC, where
# those two packages are not guaranteed to be installed, so pulling them in
# just for report text would be a fragile, avoidable new dependency.

def generate_report(adata, thresholds: dict, fit_info: dict, column_map: dict,
                     group_col: str, data_dir: str, plot_dir: str):
    import base64
    from matplotlib.backends.backend_pdf import PdfPages

    report_dir = os.path.join(os.path.dirname(os.path.normpath(plot_dir)), "cell_typing_report")
    os.makedirs(report_dir, exist_ok=True)

    total_cells = int(len(adata))
    n_markers = len(thresholds)
    ct_counts = adata.obs["cell_type"].value_counts() if "cell_type" in adata.obs.columns else pd.Series(dtype=int)
    n_cell_types = int(len(ct_counts))
    n_unassigned = int(ct_counts.get("Unassigned", 0))
    pct_unassigned = 100 * n_unassigned / total_cells if total_cells else 0.0
    now_str = datetime.now().strftime("%Y-%m-%d %H:%M")

    figures = [("cell_type_counts.png", "Cell Type Distribution"),
               ("umap_cell_types_all.png", "UMAP -- Cell Types")]
    if group_col:
        figures.append((f"umap_by_{group_col}.png", f"UMAP -- {group_col}"))
    figures += [
        ("cell_type_marker_positivity_heatmap.png", "Marker Positivity (%) by Cell Type"),
        ("cell_type_marker_expression_heatmap.png", "Mean Marker Expression by Cell Type"),
        ("cell_type_marker_dotplot.png", "Marker Expression Across Cell Types"),
    ]
    if group_col:
        figures += [
            (f"cell_type_by_{group_col}_stacked100.png", f"Cell Type Composition by {group_col}"),
            (f"cell_type_{group_col}_significant_differences.png", f"Significant Differences Across {group_col}"),
            (f"cell_type_{group_col}_difference.png", f"{group_col} Difference"),
        ]
    if n_unassigned:
        figures += [
            ("unassigned_vs_known_positivity.png", "Unassigned vs Known: Marker Positivity"),
            ("unassigned_threshold_recommendations_top10.png", "Top Threshold-Adjustment Candidates"),
        ]
    figures = [(fname, title) for fname, title in figures
               if os.path.exists(os.path.join(plot_dir, fname))]

    threshold_rows = []
    for marker, thr in thresholds.items():
        pos_col = f"{marker}_pos"
        pct_pos = 100 * adata.obs[pos_col].mean() if pos_col in adata.obs.columns else float("nan")
        threshold_rows.append((marker, thr, pct_pos))

    chi2_df = None
    if group_col:
        chi2_csv = os.path.join(data_dir, f"cell_type_{group_col}_chi2_statistics.csv")
        if os.path.exists(chi2_csv):
            chi2_df = pd.read_csv(chi2_csv, index_col=0)
            if "significant" in chi2_df.columns:
                chi2_df["significant"] = chi2_df["significant"].astype(bool)

    def _b64(fname):
        with open(os.path.join(plot_dir, fname), "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")

    html = [f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Cell Typing Analysis Report</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 24px; color: #222; }}
h1, h2, h3 {{ color: #333366; }}
.container {{ max-width: 1100px; margin: 0 auto; }}
.section {{ margin-bottom: 32px; border-bottom: 1px solid #eee; padding-bottom: 20px; }}
.figure {{ margin: 15px 0; text-align: center; }}
.figure img {{ max-width: 100%; border: 1px solid #ddd; }}
table {{ border-collapse: collapse; width: 100%; margin: 12px 0; font-size: 13px; }}
th, td {{ border: 1px solid #ddd; padding: 6px 8px; text-align: left; }}
th {{ background-color: #f2f2f2; }}
tr:nth-child(even) {{ background-color: #fafafa; }}
.highlight {{ background-color: #fff3cd; }}
.footer {{ font-size: 11px; color: #777; margin-top: 30px; text-align: center; }}
</style></head><body><div class="container">
<div class="section">
<h1>Cell Typing Analysis Report</h1>
<p><strong>Generated:</strong> {now_str}</p>
<p><strong>Total cells:</strong> {total_cells:,}</p>
<p><strong>Markers used:</strong> {n_markers}</p>
<p><strong>Cell types identified:</strong> {n_cell_types}</p>
<p><strong>Unassigned:</strong> {n_unassigned:,} cells ({pct_unassigned:.2f}%)</p>
</div>
"""]

    if figures:
        html.append('<div class="section"><h2>Overview</h2>')
        for fname, title in figures:
            html.append(f'<div class="figure"><h3>{title}</h3>'
                        f'<img src="data:image/png;base64,{_b64(fname)}" alt="{title}"></div>')
        html.append("</div>")

    html.append('<div class="section"><h2>Cell Type Counts</h2><table>'
                '<tr><th>Cell Type</th><th>Count</th><th>%</th></tr>')
    for ct, count in ct_counts.items():
        pct = 100 * count / total_cells if total_cells else 0
        cls = ' class="highlight"' if ct == "Unassigned" else ""
        html.append(f"<tr{cls}><td>{ct}</td><td>{count:,}</td><td>{pct:.2f}%</td></tr>")
    html.append("</table></div>")

    html.append('<div class="section"><h2>GMM-Detected Marker Thresholds</h2><table>'
                '<tr><th>Marker</th><th>Threshold</th><th>% Positive</th></tr>')
    for marker, thr, pct_pos in threshold_rows:
        html.append(f"<tr><td>{marker}</td><td>{thr:.3f}</td><td>{pct_pos:.2f}%</td></tr>")
    html.append("</table></div>")

    if chi2_df is not None:
        n_sig = int(chi2_df["significant"].sum()) if "significant" in chi2_df.columns else 0
        html.append(f'<div class="section"><h2>{group_col} Comparison</h2>'
                    f'<p>{n_sig} of {len(chi2_df)} cell types differ significantly '
                    f'across {group_col} (p&lt;0.05, chi-square).</p><table>'
                    '<tr><th>Cell Type</th><th>Chi-square</th><th>P-value</th><th>Significant</th></tr>')
        for ct, row in chi2_df.iterrows():
            cls = ' class="highlight"' if row.get("significant") else ""
            html.append(f"<tr{cls}><td>{ct}</td><td>{row.get('chi2', float('nan')):.3f}</td>"
                        f"<td>{row.get('p_value', float('nan')):.4f}</td>"
                        f"<td>{'Yes' if row.get('significant') else 'No'}</td></tr>")
        html.append("</table></div>")

    html.append(f'<div class="footer">Generated automatically by SPATIA cell_typing.py on {now_str}. '
                f'Full plot set (per-marker distributions/UMAPs, CSVs) in the sibling '
                f'cell_typing_plots/ and cell_typing_data/ directories.</div>')
    html.append("</div></body></html>")

    with open(os.path.join(report_dir, "cell_typing_report.html"), "w") as f:
        f.write("".join(html))

    pdf_path = os.path.join(report_dir, "cell_typing_report.pdf")
    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.88, "Cell Typing Analysis Report", ha="center", fontsize=22, weight="bold")
        fig.text(0.5, 0.82, f"Generated: {now_str}", ha="center", fontsize=12)
        fig.text(0.5, 0.78, f"Total cells: {total_cells:,}", ha="center", fontsize=12)
        fig.text(0.5, 0.75, f"Markers used: {n_markers}", ha="center", fontsize=12)
        fig.text(0.5, 0.72, f"Cell types identified: {n_cell_types}", ha="center", fontsize=12)
        fig.text(0.5, 0.69, f"Unassigned: {n_unassigned:,} ({pct_unassigned:.2f}%)", ha="center", fontsize=12)
        plt.axis("off")
        pdf.savefig(fig)
        plt.close(fig)

        for fname, title in figures:
            try:
                img = plt.imread(os.path.join(plot_dir, fname))
                fig, ax = plt.subplots(figsize=(11, 8.5))
                ax.imshow(img)
                ax.axis("off")
                ax.set_title(title, fontsize=13)
                pdf.savefig(fig)
                plt.close(fig)
            except Exception:
                plt.close("all")

    lines = [
        "CELL TYPING ANALYSIS SUMMARY", "=" * 40,
        f"Generated: {now_str}", f"Total cells: {total_cells:,}",
        f"Markers used: {n_markers}", f"Cell types identified: {n_cell_types}",
        f"Unassigned: {n_unassigned:,} ({pct_unassigned:.2f}%)",
        "", "CELL TYPE DISTRIBUTION", "-" * 40,
    ]
    for ct, count in ct_counts.items():
        pct = 100 * count / total_cells if total_cells else 0
        lines.append(f"{ct}: {count:,} cells ({pct:.2f}%)")
    if chi2_df is not None:
        lines += ["", f"{str(group_col).upper()} COMPARISON", "-" * 40]
        sig = chi2_df[chi2_df["significant"]] if "significant" in chi2_df.columns else chi2_df.iloc[0:0]
        lines.append(f"Significantly different cell types: {len(sig)}")
        for ct, row in sig.iterrows():
            lines.append(f"  * {ct}: p={row.get('p_value', float('nan')):.4f}")
    with open(os.path.join(report_dir, "cell_typing_report_summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"  [diagnostics] Report bundle written to {report_dir}/ "
          f"(cell_typing_report.html, .pdf, _summary.txt)")


# ── Entry point ──────────────────────────────────────────────────────────────

def generate_cell_typing_diagnostics(adata, thresholds: dict, fit_info: dict, column_map: dict,
                                      data_dir: str, plot_dir: str,
                                      transform: str = "none", arcsinh_cofactor: float = 5.0,
                                      group_col: str = "experiment_group",
                                      tissue_col: str = "tissue_id",
                                      random_state: int = 42, formats=("png", "pdf")):
    """
    Full diagnostic plot suite for a completed cell_typing run. Called once
    from run_cell_typing() after _save_and_plot(), for both automatic and
    semi_automatic mode. See module docstring for what's included and what
    was deliberately left out of this first pass.
    """
    print(f"\n{'=' * 80}\nGenerating cell-typing diagnostic plots...\n{'=' * 80}")

    try:
        plot_marker_distributions(adata, thresholds, fit_info, column_map, plot_dir,
                                   transform, arcsinh_cofactor, formats)
    except Exception as e:
        print(f"  WARNING: marker distribution plots failed (non-fatal): {e}")

    try:
        plot_positivity_heatmap(adata, list(thresholds.keys()), plot_dir, data_dir, formats)
    except Exception as e:
        print(f"  WARNING: positivity heatmap failed (non-fatal): {e}")

    try:
        plot_mean_expression_heatmap(adata, list(thresholds.keys()), column_map, plot_dir, data_dir, formats)
    except Exception as e:
        print(f"  WARNING: mean-expression heatmap failed (non-fatal): {e}")

    try:
        plot_cell_type_counts_pretty(adata, plot_dir, formats)
    except Exception as e:
        print(f"  WARNING: restyled cell_type_counts plot failed (non-fatal): {e}")

    try:
        plot_stacked_by_group(adata, group_col, plot_dir, data_dir, formats)
    except Exception as e:
        print(f"  WARNING: stacked cell-type-by-group plot failed (non-fatal): {e}")

    try:
        compute_and_plot_umaps(adata, plot_dir, data_dir, group_col, tissue_col, random_state, formats)
    except Exception as e:
        print(f"  WARNING: UMAP diagnostics failed (non-fatal): {e}")

    try:
        plot_dotplot(adata, thresholds, column_map, plot_dir, formats)
    except Exception as e:
        print(f"  WARNING: dotplot failed (non-fatal): {e}")

    try:
        plot_marker_umaps(adata, thresholds, column_map, plot_dir,
                           fit_info=fit_info, arcsinh_cofactor=arcsinh_cofactor, formats=formats)
    except Exception as e:
        print(f"  WARNING: per-marker UMAPs failed (non-fatal): {e}")

    try:
        plot_group_statistics(adata, group_col, plot_dir, data_dir, formats)
    except Exception as e:
        print(f"  WARNING: group statistics failed (non-fatal): {e}")

    try:
        plot_unassigned_diagnostics(adata, thresholds, column_map, plot_dir, data_dir, formats)
    except Exception as e:
        print(f"  WARNING: unassigned-cell diagnostics failed (non-fatal): {e}")

    try:
        generate_report(adata, thresholds, fit_info, column_map, group_col, data_dir, plot_dir)
    except Exception as e:
        print(f"  WARNING: report bundle generation failed (non-fatal): {e}")

    print(f"\n[diagnostics] Complete. Outputs in: {plot_dir} (figures) / {data_dir} (CSVs) / "
          f"{os.path.join(os.path.dirname(os.path.normpath(plot_dir)), 'cell_typing_report')} (report)")
