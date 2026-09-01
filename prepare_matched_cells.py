#!/usr/bin/env python3
"""
prepare_matched_cells.py
=========================
Generic, dataset-agnostic prep step: converts a raw per-cell data file
(CSV/TSV, or an AnnData .h5ad) into the "{experiment_group}_{sample}_matched_
with_boundaries.csv" files that spatia/analysis/triads.py (and the rest of
the pipeline) reads from paths.input_dir.

Why this script exists
-----------------------
triads.py never creates its own input — it only reads whatever's already in
input_dir. Two dataset-specific scripts already did this conversion by hand
(prepare_crc_data.py for the CRC TMA CSV, prepare_matusiak_data.py for the
Matusiak h5ad), each with its own hardcoded column names and relabeling
rules. This script generalizes both into one reusable tool, driven entirely
by CLI arguments, so a *new* dataset never needs a bespoke prep script of
its own — only its real column names (and, optionally, label-remapping
rules) passed as arguments. Deterministic: same input + same arguments
always produces the same output files, no guessing, no silent defaults for
anything that affects which cell goes in which file.

prepare_crc_data.py and prepare_matusiak_data.py are left as-is (already
validated against real runs) — this is the tool to reach for on any dataset
that isn't one of those two.

Phase 1 — INSPECT (always run first, on the real file — do not guess column
names):
    python prepare_matched_cells.py --inspect --input <path to .csv/.tsv/.h5ad>

    Prints columns/dtypes (tabular) or obs.columns/obsm keys/var_names
    (h5ad), plus value previews for anything that looks like a cell-type,
    experiment-group, sample, or coordinate column. Writes nothing.

Phase 2 — CONVERT (only after Phase 1 confirms real column names):
    python prepare_matched_cells.py \\
        --input  <path to .csv/.tsv/.h5ad> \\
        --output <matched_cells output dir> \\
        --cell-type-col <col> \\
        --x-col <col> --y-col <col> \\
        --experiment-group-col <col> \\
        --sample-col <col> \\
        [--experiment-group-map '{"1": "CLR", "2": "DII"}'] \\
        [--cell-type-merge-map '{"^CD4\\\\+ T cells.*": "CD4+ T cells"}'] \\
        [--force sample1,sample2]

    Writes one CSV per (experiment_group, sample):
        {experiment_group}_{sample}_matched_with_boundaries.csv
    with centroid_x, centroid_y, cell_type, experiment_group columns (all
    original columns, including marker intensities, are preserved for
    downstream functional-marker analysis).

--experiment-group-map / --cell-type-merge-map are both optional. Omit them
to use the raw column values as-is. They exist so dataset-specific
relabeling (e.g. CRC's numeric group codes -> CLR/DII, or merging CD4+ T
cell subtypes into one label) can be expressed as a declarative argument
instead of forking a new script.
"""

import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

REQUIRED_OUTPUT_COLS = ["centroid_x", "centroid_y", "cell_type", "experiment_group"]

CANDIDATE_COLS = {
    "cell type":        ["cell_type", "Cell_Type", "CellType", "cellType", "cell_type_name", "ClusterName"],
    "experiment_group": ["experiment_group", "cancer_type", "tissue_type", "Tissue", "groups", "condition"],
    "sample":           ["sample_id", "SampleID", "image_id", "Sample", "unique_region", "File Name", "Region"],
    "x coord":          ["centroid_x", "x", "X", "x_coordinate", "X_centroid", "X:X"],
    "y coord":          ["centroid_y", "y", "Y", "y_coordinate", "Y_centroid", "Y:Y"],
}


# ── Format detection ────────────────────────────────────────────────────────

def _detect_format(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    if ext in (".h5ad", ".h5"):
        return "h5ad"
    if ext in (".csv", ".tsv", ".txt"):
        return "tabular"
    print(f"ERROR: could not infer format from extension '{ext}'. "
          f"Expected .csv/.tsv/.txt or .h5ad/.h5 — pass --format explicitly.",
          file=sys.stderr)
    sys.exit(1)


def _sep_for(path: str, explicit_sep: str = None) -> str:
    if explicit_sep is not None:
        return explicit_sep
    return "\t" if os.path.splitext(path)[1].lower() == ".tsv" else ","


# ── Inspect mode ─────────────────────────────────────────────────────────────

def _print_candidate_previews(get_col, columns) -> None:
    for label, candidates in CANDIDATE_COLS.items():
        found = [c for c in candidates if c in columns]
        if not found:
            print(f"\n[{label}] no candidate column found among {candidates}")
            print(f"          -> inspect the full column list above and pass the real name explicitly.")
            continue
        for col in found:
            series = get_col(col)
            print(f"\n[{label}] '{col}'  dtype={series.dtype}")
            if series.dtype == object or str(series.dtype) == "category":
                print(series.value_counts().head(30).to_string())
            else:
                print(series.describe().to_string())


def inspect_tabular(path: str, sep: str) -> None:
    print(f"[inspect] Reading {path} ...")
    df = pd.read_csv(path, sep=sep, nrows=200_000)
    print(f"[inspect] {len(df):,} rows read (may be truncated for inspection) x {len(df.columns)} columns\n")
    print("-" * 70)
    print("Columns")
    print("-" * 70)
    print(list(df.columns))
    print("\n" + "-" * 70)
    print("Candidate column matches + value previews")
    print("-" * 70)
    _print_candidate_previews(lambda c: df[c], df.columns)
    print("\n[inspect] Done. Nothing written. Re-run without --inspect, passing "
          "--cell-type-col / --x-col / --y-col / --experiment-group-col / "
          "--sample-col with CONFIRMED names.")


def inspect_h5ad(path: str) -> None:
    try:
        import anndata as ad
    except ImportError:
        print("ERROR: anndata not installed. `pip install anndata --break-system-packages`",
              file=sys.stderr)
        sys.exit(1)

    print(f"[inspect] Reading {path} ...")
    adata = ad.read_h5ad(path)
    print(f"[inspect] n_obs={adata.n_obs:,}  n_vars={adata.n_vars}\n")

    print("-" * 70)
    print("var_names (marker panel) —", adata.n_vars, "markers")
    print("-" * 70)
    print(list(adata.var_names))

    print("\n" + "-" * 70)
    print("obs.columns")
    print("-" * 70)
    print(list(adata.obs.columns))

    print("\n" + "-" * 70)
    print("obsm keys (coordinates sometimes live here instead of obs)")
    print("-" * 70)
    for k in adata.obsm.keys():
        try:
            print(f"  {k}: shape={adata.obsm[k].shape}")
        except Exception:
            print(f"  {k}")

    print("\n" + "-" * 70)
    print("Candidate column matches + value previews")
    print("-" * 70)
    _print_candidate_previews(lambda c: adata.obs[c], adata.obs.columns)

    print("\n[inspect] Done. Nothing written. Re-run without --inspect, passing "
          "--cell-type-col / --x-col / --y-col / --experiment-group-col / "
          "--sample-col with CONFIRMED names.")


# ── Shared relabeling + write logic ──────────────────────────────────────────

def _apply_experiment_group_map(df: pd.DataFrame, group_map: dict) -> None:
    if not group_map:
        return
    df["experiment_group"] = df["experiment_group"].astype(str).map(
        lambda v: group_map.get(v, v)
    )


def _apply_cell_type_merge_map(df: pd.DataFrame, merge_map: dict) -> None:
    """
    merge_map: {regex_pattern: replacement}, applied in order via a full
    regex substitution over the cell_type column. E.g.
        {"^CD4\\+ T cells.*": "CD4+ T cells"}
    merges every "CD4+ T cells (...)" variant into one label. Original
    values are NOT preserved separately here — if you need the raw label
    kept, don't rename the source column before calling this, or add your
    own '_raw_cell_type' copy before running this script.
    """
    if not merge_map:
        return
    for pattern, replacement in merge_map.items():
        before = df["cell_type"].astype(str)
        after = before.str.replace(pattern, replacement, regex=True)
        n_changed = (before != after).sum()
        df["cell_type"] = after
        if n_changed:
            print(f"[prep] Cell-type merge '{pattern}' -> '{replacement}': {n_changed:,} rows changed")


def _validate_required_columns(df: pd.DataFrame) -> None:
    missing = [c for c in REQUIRED_OUTPUT_COLS if c not in df.columns]
    if missing:
        print(f"ERROR: after renaming, these required output columns are still "
              f"missing: {missing}. Check --cell-type-col/--x-col/--y-col/"
              f"--experiment-group-col map to real columns (see --inspect output).",
              file=sys.stderr)
        sys.exit(1)
    n_null = df[REQUIRED_OUTPUT_COLS].isnull().any(axis=1).sum()
    if n_null:
        print(f"[prep] ⚠️  {n_null:,} row(s) have a null value in one of "
              f"{REQUIRED_OUTPUT_COLS} — these rows will still be written "
              f"(triads.py may drop or mishandle them; worth checking the "
              f"source data if this count is unexpectedly high).")


def _write_per_sample_csvs(df: pd.DataFrame, output_dir: str, sample_col: str,
                            force: set) -> int:
    os.makedirs(output_dir, exist_ok=True)

    if "cell_id" not in df.columns:
        df = df.copy()
        df.insert(0, "cell_id", df.index.astype(str))

    n_written, n_skipped = 0, 0
    # Group by (experiment_group, sample) — not sample alone. Sample/region
    # identifiers are sometimes reused across experiment_groups (e.g. two
    # different physical cores both named "reg005"); grouping on sample
    # alone would silently merge unrelated tissue into one file.
    group_cols = ["experiment_group", sample_col] if sample_col != "experiment_group" else ["experiment_group"]
    for group_key, group in df.groupby(group_cols):
        if isinstance(group_key, tuple):
            experiment_group, sample_id = group_key
        else:
            experiment_group, sample_id = group_key, group_key
        if len(group) == 0:
            continue

        safe_sample = str(sample_id).replace("/", "-").replace(" ", "_")
        out_name = f"{experiment_group}_{safe_sample}_matched_with_boundaries.csv"
        out_path = os.path.join(output_dir, out_name)

        already_done = (
            os.path.exists(out_path) and os.path.getsize(out_path) > 0
            and out_name not in force and str(sample_id) not in force
        )
        if already_done:
            n_skipped += 1
            continue

        group.reset_index(drop=True).to_csv(out_path, index=False)
        n_written += 1

    print(f"\n[prep] Wrote {n_written} file(s), skipped {n_skipped} already-present "
          f"file(s) -> {output_dir}")
    sample_files = sorted(os.listdir(output_dir))[:5]
    print(f"[prep] Sample files: {sample_files}")
    print("\n[prep] cell_type values written:")
    print(df["cell_type"].value_counts().to_string())
    print("\n[prep] experiment_group values written:")
    print(df["experiment_group"].value_counts().to_string())
    return n_written


# ── Convert: tabular (CSV/TSV) ───────────────────────────────────────────────

def convert_tabular(
    path: str, sep: str, output_dir: str,
    cell_type_col: str, x_col: str, y_col: str,
    experiment_group_col: str, sample_col: str,
    experiment_group_map: dict, cell_type_merge_map: dict, force: set,
) -> None:
    print(f"[prep] Reading {path} ...")
    df = pd.read_csv(path, sep=sep)
    print(f"[prep] {len(df):,} rows read")

    required = {
        "cell_type_col": cell_type_col, "x_col": x_col, "y_col": y_col,
        "experiment_group_col": experiment_group_col, "sample_col": sample_col,
    }
    missing = {k: v for k, v in required.items() if v not in df.columns}
    if missing:
        print(f"ERROR: these columns were not found in the input: {missing}", file=sys.stderr)
        print("Run with --inspect first to see the real column names.", file=sys.stderr)
        sys.exit(1)

    df = df.rename(columns={
        x_col: "centroid_x", y_col: "centroid_y",
        cell_type_col: "cell_type", experiment_group_col: "experiment_group",
    })

    _apply_experiment_group_map(df, experiment_group_map)
    _apply_cell_type_merge_map(df, cell_type_merge_map)
    _validate_required_columns(df)
    _write_per_sample_csvs(df, output_dir, sample_col, force)


# ── Convert: AnnData h5ad ────────────────────────────────────────────────────

def convert_h5ad(
    path: str, output_dir: str,
    cell_type_col: str, x_col: str, y_col: str,
    experiment_group_col: str, sample_col: str,
    experiment_group_map: dict, cell_type_merge_map: dict, force: set,
) -> None:
    try:
        import anndata as ad
    except ImportError:
        print("ERROR: anndata not installed. `pip install anndata --break-system-packages`",
              file=sys.stderr)
        sys.exit(1)

    print(f"[prep] Reading {path} ...")
    adata = ad.read_h5ad(path)
    print(f"[prep] {adata.n_obs:,} cells x {adata.n_vars} markers")

    # x/y are allowed to come from adata.obs OR fall through to obsm lookup
    # by the caller before this function runs (kept simple here: obs only,
    # matching --inspect's obsm listing so the user can pick the right one
    # and pass it in as a plain obs column if needed, or extend this
    # function directly if coordinates genuinely only live in obsm).
    required = {
        "cell_type_col": cell_type_col, "x_col": x_col, "y_col": y_col,
        "experiment_group_col": experiment_group_col, "sample_col": sample_col,
    }
    missing = {k: v for k, v in required.items() if v not in adata.obs.columns}
    if missing:
        print(f"ERROR: these obs columns were not found: {missing}", file=sys.stderr)
        print("Run with --inspect first to see the real column names (and check "
              "the obsm listing if coordinates live there instead).", file=sys.stderr)
        sys.exit(1)

    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    marker_df = pd.DataFrame(np.asarray(X), columns=list(adata.var_names), index=adata.obs.index)
    df = pd.concat([adata.obs.reset_index(drop=True), marker_df.reset_index(drop=True)], axis=1)

    df = df.rename(columns={
        x_col: "centroid_x", y_col: "centroid_y",
        cell_type_col: "cell_type", experiment_group_col: "experiment_group",
    })

    _apply_experiment_group_map(df, experiment_group_map)
    _apply_cell_type_merge_map(df, cell_type_merge_map)
    _validate_required_columns(df)
    _write_per_sample_csvs(df, output_dir, sample_col, force)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generic prep: raw per-cell data -> *_matched_with_boundaries.csv for triads.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", "-i", required=True, help="Path to .csv/.tsv or .h5ad")
    parser.add_argument("--format", choices=["tabular", "h5ad"], default=None,
                         help="Override auto-detection from file extension")
    parser.add_argument("--sep", default=None, help="Delimiter for tabular input (default: ',' unless .tsv)")
    parser.add_argument("--inspect", action="store_true",
                         help="Print column structure and candidate matches; write nothing.")
    parser.add_argument("--output", "-o", default="data/matched_cells",
                         help="Output directory for per-image CSVs (Convert mode only)")
    parser.add_argument("--cell-type-col", default=None)
    parser.add_argument("--x-col", default=None)
    parser.add_argument("--y-col", default=None)
    parser.add_argument("--experiment-group-col", default=None)
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--experiment-group-map", default=None,
                         help='Optional JSON, e.g. \'{"1": "CLR", "2": "DII"}\'')
    parser.add_argument("--cell-type-merge-map", default=None,
                         help='Optional JSON of {regex: replacement}, e.g. '
                              '\'{"^CD4\\\\+ T cells.*": "CD4+ T cells"}\'')
    parser.add_argument("--force", default="",
                         help="Comma-separated sample values (or full output filenames) "
                              "to force-rewrite even if already present on disk")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    fmt = args.format or _detect_format(args.input)

    if args.inspect:
        if fmt == "tabular":
            inspect_tabular(args.input, _sep_for(args.input, args.sep))
        else:
            inspect_h5ad(args.input)
        return

    required = [args.cell_type_col, args.x_col, args.y_col, args.experiment_group_col, args.sample_col]
    if any(r is None for r in required):
        print("ERROR: --cell-type-col/--x-col/--y-col/--experiment-group-col/--sample-col "
              "are all required for conversion. Run with --inspect first to find the "
              "real column names — do not guess.", file=sys.stderr)
        sys.exit(1)

    experiment_group_map = json.loads(args.experiment_group_map) if args.experiment_group_map else {}
    cell_type_merge_map = json.loads(args.cell_type_merge_map) if args.cell_type_merge_map else {}
    force = {s.strip() for s in args.force.split(",") if s.strip()}

    if fmt == "tabular":
        convert_tabular(
            args.input, _sep_for(args.input, args.sep), args.output,
            cell_type_col=args.cell_type_col, x_col=args.x_col, y_col=args.y_col,
            experiment_group_col=args.experiment_group_col, sample_col=args.sample_col,
            experiment_group_map=experiment_group_map, cell_type_merge_map=cell_type_merge_map,
            force=force,
        )
    else:
        convert_h5ad(
            args.input, args.output,
            cell_type_col=args.cell_type_col, x_col=args.x_col, y_col=args.y_col,
            experiment_group_col=args.experiment_group_col, sample_col=args.sample_col,
            experiment_group_map=experiment_group_map, cell_type_merge_map=cell_type_merge_map,
            force=force,
        )


if __name__ == "__main__":
    main()
