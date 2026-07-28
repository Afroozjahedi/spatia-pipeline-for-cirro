"""
spatia.analysis.agentic
=======================
Agentic (LLM-driven) decision layer for the SPATIA pipeline.

This module closes the human-in-the-loop gaps in `cell_typing.py`. In
`semi_automatic` mode the pipeline currently STOPS after Leiden clustering and
waits for a human to inspect UMAPs and hand-fill `cluster_labels.yaml`. Here an
LLM reads each cluster's quantitative marker profile and proposes that mapping —
with a confidence and a one-line rationale for every call, so the decision stays
auditable rather than opaque.

Design principles
-----------------
1. **Propose, don't overwrite.** The agent writes a *candidate*
   `cluster_labels.yaml` (with rationale/confidence as comments). The existing
   Phase-2 code path consumes it unchanged, and you can edit any call before the
   re-run. Nothing downstream is bypassed.
2. **Grounded, not from memory.** The prompt contains only numbers computed from
   the data (z-scored per-cluster marker means, cluster sizes, % positive). The
   LLM is constrained to a vocabulary you supply.
3. **Provider seam.** Runs against the `anthropic` SDK on Cirro/Seadragon
   (reads ANTHROPIC_API_KEY), OR against a `host.llm`-style callable when run
   inside Claude Science. Inject any callable via `llm_fn=` for testing.

Usage
-----
    # As a pipeline step (after semi_automatic Phase 1 clustering):
    from spatia.analysis.agentic import agentic_label_clusters
    agentic_label_clusters(cfg)     # writes candidate cluster_labels.yaml
    # then re-run:  python run_pipeline.py -c <cfg> --steps cell_typing

    # Standalone on any clustered h5ad:
    agentic_label_clusters(cfg, clustered_h5ad="…/X_clustered.h5ad")
"""
from __future__ import annotations

import json
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd
import yaml


# ── LLM provider seam ─────────────────────────────────────────────────────────

def _default_llm_fn(prompt: str, model: str | None = None) -> str:
    """Call Anthropic's API and return the text. Reads ANTHROPIC_API_KEY.

    Swappable: pass your own `llm_fn(prompt, model)->str` to agentic_* functions
    (e.g. a `host.llm` wrapper inside Claude Science, or a mock in tests).
    """
    try:
        import anthropic
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "The agentic layer needs the 'anthropic' package "
            "(`pip install anthropic`) or an injected llm_fn=."
        ) from e
    client = anthropic.Anthropic()  # ANTHROPIC_API_KEY from env
    msg = client.messages.create(
        model=model or "claude-sonnet-4-5",
        max_tokens=2000,
        messages=[{"role": "user", "content": prompt}],
    )
    return "".join(block.text for block in msg.content if block.type == "text")


def _extract_json(text: str) -> dict:
    """Pull the first JSON object out of an LLM response."""
    m = re.search(r"\{.*\}", text, re.DOTALL)
    if not m:
        raise ValueError(f"No JSON object found in LLM response:\n{text[:500]}")
    return json.loads(m.group(0))


# ── Cluster marker profiles → prompt ──────────────────────────────────────────

def cluster_marker_profiles(adata, cluster_key="leiden", markers=None,
                            n_top=8, n_bottom=3) -> dict:
    """Per-cluster z-scored marker profile the LLM can read.

    z-scoring is ACROSS clusters, so the LLM sees each cluster's most
    DISTINCTIVE markers rather than globally abundant ones. Also reports
    cluster size and, when '<marker>_pos' columns exist, the % positive.

    Returns {cluster_id: {n, pct, high:[(m,z)], low:[(m,z)]}}.
    """
    if markers is None:
        # Use the marker panel present in the data (exclude engineered columns).
        markers = [m for m in adata.var_names]
    markers = [m for m in markers if m in adata.var_names]

    X = adata[:, markers].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    df = pd.DataFrame(X, columns=markers)
    df[cluster_key] = adata.obs[cluster_key].values

    cl_mean = df.groupby(cluster_key, observed=True).mean()
    cl_z = (cl_mean - cl_mean.mean()) / (cl_mean.std() + 1e-9)
    sizes = adata.obs[cluster_key].value_counts()

    # Optional: % positive per cluster from existing *_pos columns
    pos_cols = {m: f"{m}_pos" for m in markers if f"{m}_pos" in adata.obs.columns}

    out = {}
    for cl in cl_z.index:
        row = cl_z.loc[cl].sort_values(ascending=False)
        pct = {}
        if pos_cols:
            mask = adata.obs[cluster_key].values == cl
            for m, col in pos_cols.items():
                pct[m] = round(float(adata.obs.loc[mask, col].mean()) * 100, 1)
        out[str(cl)] = {
            "n": int(sizes[cl]),
            "high": [(m, round(float(row[m]), 2)) for m in row.index[:n_top]],
            "low": [(m, round(float(row[m]), 2)) for m in row.index[-n_bottom:]],
            "pct_pos": pct,
        }
    return out


def build_cluster_prompt(profiles: dict, markers: list, vocabulary: list,
                         context: str = "") -> str:
    """Construct the LLM prompt for cluster → cell-type labeling."""
    lines = []
    for cl, r in profiles.items():
        hi = ", ".join(f"{m}={v}" for m, v in r["high"])
        lo = ", ".join(f"{m}={v}" for m, v in r["low"])
        extra = ""
        if r["pct_pos"]:
            top_pos = sorted(r["pct_pos"].items(), key=lambda kv: -kv[1])[:6]
            extra = "  %pos[" + ", ".join(f"{m}={v}%" for m, v in top_pos) + "]"
        lines.append(f"Cluster {cl} (n={r['n']}): HIGH[{hi}] LOW[{lo}]{extra}")

    vocab = ", ".join(f"'{v}'" for v in vocabulary)
    ctx = f"\nStudy context: {context}\n" if context else ""
    return (
        "You are an expert in imaging spatial proteomics (CODEX / IMC / CyCIF) "
        "annotating cell clusters from marker-expression profiles.\n"
        f"{ctx}"
        f"Marker panel: {', '.join(markers)}\n"
        "Below are Leiden clusters with their most DISTINCTIVELY high/low markers "
        "(z-scored across clusters); %pos is the fraction positive per marker.\n\n"
        f"Assign each cluster ONE label. Prefer labels from this list: {vocab}. "
        "Use 'Unassigned' only if no marker signal is interpretable.\n"
        "Weigh canonical lineage markers over ambiguous ones; note when a cluster "
        "looks like a doublet or low-quality (few distinctive markers).\n\n"
        "Return ONLY JSON mapping each cluster id (string) to an object "
        '{"label": <str>, "confidence": "high|medium|low", "reason": "<=15 words"}.\n\n'
        "Clusters:\n" + "\n".join(lines)
    )


# ── Main entry point ──────────────────────────────────────────────────────────

def agentic_label_clusters(cfg: dict, clustered_h5ad: str | None = None,
                           llm_fn=None, model: str | None = None,
                           overwrite: bool = False) -> str:
    """Read a clustered h5ad, ask an LLM to label each Leiden cluster, and write
    a CANDIDATE cluster_labels.yaml that Phase 2 of cell_typing consumes.

    Parameters
    ----------
    cfg : dict                  loaded experiment config
    clustered_h5ad : str        path to <analysis>_clustered.h5ad
                                (defaults to the semi_auto Phase-1 output path)
    llm_fn : callable           llm_fn(prompt, model)->str; defaults to Anthropic
    model : str                 model id passed to llm_fn
    overwrite : bool            refuse to clobber an existing labels file unless True

    Returns the path to the written cluster_labels.yaml.
    """
    import scanpy as sc

    ct_cfg = cfg["cell_typing"]
    analysis = ct_cfg.get("analysis_name", cfg["experiment"]["name"])
    out_dir = cfg["paths"]["output_dir"]
    data_dir = os.path.join(out_dir, "cell_typing_data")

    if clustered_h5ad is None:
        clustered_h5ad = os.path.join(data_dir, f"{analysis}_clustered.h5ad")
    if not os.path.exists(clustered_h5ad):
        raise FileNotFoundError(
            f"Clustered h5ad not found: {clustered_h5ad}\n"
            "Run semi_automatic Phase 1 (clustering) first."
        )

    labels_path = ct_cfg.get("cluster_labels_file") or \
        os.path.join(data_dir, f"{analysis}_cluster_labels.yaml")
    if os.path.exists(labels_path) and not overwrite:
        raise FileExistsError(
            f"{labels_path} already exists. Pass overwrite=True to replace it "
            "(existing human edits would be lost)."
        )

    print(f"[agentic] Loading clustered data: {clustered_h5ad}")
    adata = sc.read(clustered_h5ad)
    n_clusters = adata.obs["leiden"].nunique()
    print(f"[agentic] {adata.n_obs:,} cells in {n_clusters} Leiden clusters")

    # Vocabulary: cell types the pipeline knows about, from the definitions file
    # if present, else a sensible immune/stromal/tumor default.
    vocab = _vocabulary_from_cfg(cfg)
    markers = [m for m in ct_cfg.get("markers", {}).get("panel", [])
               if m in adata.var_names] or list(adata.var_names)

    profiles = cluster_marker_profiles(adata, markers=markers)
    context = cfg.get("experiment", {}).get("description", "")
    prompt = build_cluster_prompt(profiles, markers, vocab, context=context)

    print(f"[agentic] Querying LLM to label {n_clusters} clusters...")
    llm_fn = llm_fn or _default_llm_fn
    calls = _extract_json(llm_fn(prompt, model))

    # Assemble YAML with rationale/confidence preserved as comments.
    cluster_labels = {str(k): v["label"] for k, v in calls.items()}
    n_low = sum(1 for v in calls.values() if v.get("confidence") == "low")

    header = [
        "# ============================================================",
        "# cluster_labels.yaml  —  AGENT-PROPOSED (review before use)",
        f"# experiment : {cfg['experiment']['name']}",
        f"# generated  : {datetime.now().isoformat(timespec='seconds')}",
        f"# model      : {model or 'default'}",
        f"# clusters   : {n_clusters}   low-confidence: {n_low}",
        "# Each line below shows the agent's confidence + reason as a comment.",
        "# Edit any label before re-running cell_typing.",
        "# ============================================================",
        "cluster_labels:",
    ]
    body = []
    for cl in sorted(calls, key=lambda x: int(x) if x.isdigit() else x):
        v = calls[cl]
        body.append(
            f'  "{cl}": "{v["label"]}"'
            f'   # [{v.get("confidence","?")}] {v.get("reason","")}'
        )

    os.makedirs(os.path.dirname(labels_path), exist_ok=True)
    with open(labels_path, "w") as f:
        f.write("\n".join(header + body) + "\n")

    # Also dump the raw calls as JSON for provenance / downstream audit.
    audit_path = os.path.join(data_dir, f"{analysis}_agentic_calls.json")
    with open(audit_path, "w") as f:
        json.dump({"profiles": profiles, "calls": calls}, f, indent=2)

    print(f"[agentic] ✓ Candidate labels written: {labels_path}")
    print(f"[agentic] ✓ Audit trail written:      {audit_path}")
    if n_low:
        print(f"[agentic] ⚠ {n_low} cluster(s) flagged low-confidence — review these first.")
    print("\n[agentic] NEXT: review the file, set "
          "cell_typing.cluster_labels_file in your config, then re-run:")
    print(f"  python run_pipeline.py -c <config> --steps cell_typing")
    return labels_path


def _vocabulary_from_cfg(cfg: dict) -> list:
    """Cell-type vocabulary for the LLM: keys of the definitions file if present,
    plus the triad anchor/partner types, else a generic default."""
    vocab = []
    defs_file = cfg.get("cell_typing", {}).get("cell_type_definitions_file")
    if defs_file and os.path.exists(defs_file):
        with open(defs_file) as f:
            raw = yaml.safe_load(f) or {}
        vocab = list((raw.get("cell_type_definitions", raw) or {}).keys())
    tri = cfg.get("analysis", {}).get("triad", {})
    for k in ("anchor_type", "partner_type_1", "partner_type_2"):
        if tri.get(k) and tri[k] not in vocab:
            vocab.append(tri[k])
    if not vocab:
        vocab = ["tumor cells", "CD4+ T cells", "CD8+ T cells", "B cells",
                 "CD11c+ DCs", "macrophages", "NK cells", "endothelial",
                 "fibroblasts", "Treg", "granulocytes", "Unassigned"]
    return vocab


# ── Agentic threshold QC (insertion point #3) ─────────────────────────────────

def agentic_threshold_qc(cfg: dict, llm_fn=None, model: str | None = None) -> dict:
    """Read the GMM marker_thresholds.csv and ask the LLM to flag implausible
    thresholds (e.g. the CD11c drift that inverts KO/WT signal). Advisory only —
    returns {marker: {verdict, note}} and writes a QC report. Does not modify
    thresholds.
    """
    out_dir = cfg["paths"]["output_dir"]
    data_dir = os.path.join(out_dir, "cell_typing_data")
    thr_path = os.path.join(data_dir, "marker_thresholds.csv")
    if not os.path.exists(thr_path):
        raise FileNotFoundError(f"{thr_path} not found — run cell_typing first.")

    thr = pd.read_csv(thr_path, index_col=0)
    lines = [f"{m}: threshold={row['threshold']:.4f}"
             for m, row in thr.iterrows()]
    prompt = (
        "You are QC-reviewing GMM-derived positivity thresholds for a spatial "
        "proteomics marker panel. Flag any threshold that looks implausibly "
        "high (would call almost nothing positive) or low (would call almost "
        "everything positive) for its marker's known biology.\n"
        "Return ONLY JSON {marker: {\"verdict\":\"ok|suspect\","
        "\"note\":\"<=15 words\"}}.\n\nThresholds:\n" + "\n".join(lines)
    )
    llm_fn = llm_fn or _default_llm_fn
    verdicts = _extract_json(llm_fn(prompt, model))
    report_path = os.path.join(data_dir, "agentic_threshold_qc.json")
    with open(report_path, "w") as f:
        json.dump(verdicts, f, indent=2)
    n_suspect = sum(1 for v in verdicts.values() if v.get("verdict") == "suspect")
    print(f"[agentic] Threshold QC written: {report_path} "
          f"({n_suspect} suspect marker(s))")
    return verdicts
