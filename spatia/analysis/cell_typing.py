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

    Zero-inflation note: when transform != "none", exact-zero raw values
    are excluded from the GMM fit (still included in "clean" for the
    fallback). A monotonic transform (arcsinh/log1p/yeojohnson) reshapes
    smooth right-skew, but it cannot fix a genuine point mass at exactly
    zero -- every one of these transforms maps 0 -> 0, so a real "30-50%
    of cells have exact-zero intensity" spike (documented for 14/53 CRC
    markers) survives the transform unchanged and can still make one GMM
    component collapse onto it, exactly the degeneracy the transform was
    meant to fix. Excluding true zeros and fitting the 2-component GMM on
    the remaining continuum is the standard complement to transforming.
    transform="none" is left untouched (no zero exclusion) to keep that
    path exactly backward compatible with pre-2026-08-19 behavior.
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
        fit_input = clean[clean.flatten() != 0].reshape(-1, 1)
        if len(fit_input) < 10:
            print(f"    [GMM] WARNING: fewer than 10 non-zero values after excluding exact "
                  f"zeros for transform='{transform}' — falling back to percentile threshold")
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
                               confidence_overrides: Optional[dict] = None) -> Tuple[dict, dict]:
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
        "threshold_mode", "confidence_level"}} -- everything
        add_positivity_columns needs to make the actual per-cell call,
        including the fitted GMM itself for "posterior" mode.

    confidence_overrides: per-marker gmm.confidence_level overrides,
        mirroring how std_multipliers/per_marker_overrides already work --
        only used when threshold_mode == "posterior".
    """
    confidence_overrides = confidence_overrides or {}
    thresholds: Dict[str, float] = {}
    fit_info: Dict[str, dict] = {}

    for m in markers:
        if m not in adata.var_names:
            print(f"    [GMM] marker '{m}' not in data — skipping")
            continue
        vals = adata[:, m].X
        if hasattr(vals, "toarray"):
            vals = vals.toarray().flatten()
        else:
            vals = np.array(vals).flatten()

        fit = _fit_gmm(vals, n_components=n_components, random_state=random_state,
                       n_init=n_init, max_cells=max_cells, transform=transform,
                       arcsinh_cofactor=arcsinh_cofactor)

        mult = std_multipliers.get(m, default_std)
        conf = confidence_overrides.get(m, confidence_level)

        if threshold_mode == "posterior":
            t = _gmm_posterior_effective_threshold(fit, conf, transform, arcsinh_cofactor)
            lam_note = f"  lambda={fit['transform_params']['lambda']:.3f}" if "lambda" in fit["transform_params"] else ""
            print(f"    {m:<20} mode=posterior  confidence={conf:.2f}  "
                  f"effective_threshold={t:.4f}{lam_note}")
        else:
            t = _gmm_threshold_from_fit(fit, mult, transform, arcsinh_cofactor)
            lam_note = f"  lambda={fit['transform_params']['lambda']:.3f}" if "lambda" in fit["transform_params"] else ""
            print(f"    {m:<20} mode=std_multiplier  std_mult={mult:.1f}  threshold={t:.4f}{lam_note}")

        thresholds[m] = t
        fit_info[m] = {
            **fit,
            "threshold_mode": threshold_mode,
            "confidence_level": conf,
        }

    return thresholds, fit_info


def add_positivity_columns(adata, thresholds: dict, fit_info: Optional[dict] = None,
                            transform: str = "none", arcsinh_cofactor: float = 5.0) -> None:
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
    """
    fit_info = fit_info or {}
    new_cols: dict[str, np.ndarray] = {}

    for marker, threshold in thresholds.items():
        if marker not in adata.var_names:
            continue
        vals = adata[:, marker].X
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
            vals_t, _ = _transform_values(vals.reshape(-1, 1), transform, arcsinh_cofactor,
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
                        confidence_level: float = 0.8):
    """
    Returns adata filtered to CD45+ cells.
    Saves a threshold histogram to plot_dir if provided.

    threshold_mode/confidence_level: same std_multiplier/posterior choice
    as the per-marker thresholds (see module docstring), applied here for
    consistency rather than leaving CD45 gating permanently on the older
    std_multiplier-only path while other markers can use posteriors.

    transform/arcsinh_cofactor: Q16 remediation, see `_transform_values()`.
    """
    if "CD45" not in adata.var_names:
        print("  [CD45 gate] WARNING: CD45 not found — returning all cells")
        return adata

    vals = adata[:, "CD45"].X
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
        )

    # ── Markers for analysis (exclude gating-only markers) ───
    markers_for_analysis = [m for m in panel if m not in gating_only and m in adata_cd45.var_names]
    missing = [m for m in panel if m not in gating_only and m not in adata_cd45.var_names]
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
    )
    threshold_report = pd.DataFrame({
        "threshold": thresholds,
        "threshold_mode": {m: threshold_mode for m in thresholds},
        "transform": {m: gmm_transform for m in thresholds},
        "transform_lambda": {
            m: fit_info[m]["transform_params"].get("lambda") for m in thresholds
        },
    })
    threshold_report.index.name = "marker"
    threshold_report.to_csv(os.path.join(data_dir, "marker_thresholds.csv"))

    # ── Add positivity columns ────────────────────────────────
    add_positivity_columns(adata_cd45, thresholds, fit_info=fit_info,
                            transform=gmm_transform, arcsinh_cofactor=arcsinh_cofactor)

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
