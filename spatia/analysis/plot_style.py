"""
Shared plotting style for every SPATIA step (added 2026-11-05).

One place decides the colour of each experiment group, so CLR / DII (or any
other groups) look the same in preprocessing, cell typing, triads,
functional and survival plots.

Config (optional, any step's YAML):

    plotting:
      experiment_group_colors:
        CLR: "#2166AC"    # blue
        DII: "#B2182B"    # red

Groups not listed get colours from a fixed fallback list in the order they
are first seen, and keep that colour for the rest of the run.
"""
from typing import Iterable, Optional

import numpy as np

DEFAULT_GROUP_COLORS = {"CLR": "#2166AC", "DII": "#B2182B"}
_FALLBACK = ["#2166AC", "#B2182B", "#1B7837", "#E08214", "#762A83", "#4D4D4D"]

# Mutable on purpose: modules import this dict once and see later updates.
GROUP_COLORS: dict = dict(DEFAULT_GROUP_COLORS)


def set_group_colors(cfg: Optional[dict]) -> dict:
    """Load plotting.experiment_group_colors from a run config (if any)."""
    user = ((cfg or {}).get("plotting") or {}).get("experiment_group_colors") or {}
    GROUP_COLORS.update({str(k): str(v) for k, v in user.items()})
    return GROUP_COLORS


def group_color(group, groups_in_order: Optional[list] = None) -> str:
    g = str(group)
    if g in GROUP_COLORS:
        return GROUP_COLORS[g]
    used = set(GROUP_COLORS.values())
    free = [c for c in _FALLBACK if c not in used] or _FALLBACK
    GROUP_COLORS[g] = free[0]
    return GROUP_COLORS[g]


def group_palette(groups: Iterable) -> dict:
    return {str(g): group_color(g) for g in groups}


def p_to_stars(p) -> str:
    if p is None or not np.isfinite(p):
        return "n/a"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def add_sig_bracket(ax, x1, x2, y, p, h=None, label_fmt="{stars} p={p:.2g}", fontsize=9):
    """Draw a bracket between x1 and x2 at height y with the p-value."""
    if h is None:
        lo, hi = ax.get_ylim()
        h = (hi - lo) * 0.03
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", linewidth=1)
    ax.text((x1 + x2) / 2, y + h, label_fmt.format(stars=p_to_stars(p), p=p),
            ha="center", va="bottom", fontsize=fontsize)
    lo, hi = ax.get_ylim()
    if y + 4 * h > hi:
        ax.set_ylim(lo, y + 4 * h)


def bh_fdr(pvals) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-safe)."""
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out
