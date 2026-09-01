#!/usr/bin/env python3
"""
pool_tissues_for_celltyping.py
===============================
Concatenates preprocessing.py's per-tissue combined h5ads
({tissue_id}_combined_all_experiment_groups.h5ad) into one cohort-level
h5ad, ready to be pointed at as cell_typing.py's `input_file`.

Why this exists
----------------
preprocessing.py pools *within* a tissue (across experiment_groups/crops that
share one tissue_id -- see extract_tissue_identifier(), which strips the
experiment_group label and coordinate block off image_id) -- not *across*
tissues. If your cohort spans multiple distinct tissue_ids (different
timepoints, different patients, different TMA blocks, etc.), preprocessing.py
writes one combined h5ad per tissue, and nothing else in the pipeline merges
them further.

cell_typing.py's `input_file` is a single hardcoded path -- it does not glob
or loop over multiple files. Point it at one tissue's h5ad and GMM thresholds
(and, in semi_automatic mode, Leiden clusters) are fit only on that tissue's
cells -- meaning cell-type calls aren't necessarily comparable across tissues
typed separately (any batch drift in staining/imaging between tissues gets
baked into each tissue's own thresholds). Run this script first if you want
thresholds/clusters fit jointly across every tissue in your cohort, then
point cell_typing.input_file at this script's --output-file instead.

Dataset-agnostic: pools however many tissue files are in --input-dir, however
many there are, named however extract_tissue_identifier() produced them.
Not specific to any one study's day/timepoint/patient naming scheme.

Usage
-----
    python pool_tissues_for_celltyping.py \\
        --input-dir /path/to/combined_processed_data/individual_processed_data \\
        --output-file /path/to/combined_processed_data/cohort_pooled.h5ad

Marker panel handling
----------------------
Tissues are joined on marker (var) name via anndata.concat(). Default
--join inner keeps only markers present in EVERY tissue file, so GMM
threshold fitting downstream never sees NaN-padded columns -- if tissue
marker panels differ, this prints exactly which markers get dropped and
from which tissues. Pass --join outer to keep the union of all markers
instead (NaN-padded for tissues missing a given marker) -- only do this if
you plan to handle the NaNs yourself before cell_typing.

A `tissue_id` column is added to adata.obs (derived from each input
filename) so provenance survives pooling -- useful for later checking
whether cell-type composition differs by tissue.
"""

import argparse
import glob
import os
import sys

import anndata as ad


def _tissue_id_from_filename(path: str) -> str:
    base = os.path.basename(path)
    return base.replace("_combined_all_experiment_groups.h5ad", "")


def pool_tissues(input_dir: str, output_file: str, glob_pattern: str, join: str) -> ad.AnnData:
    files = sorted(glob.glob(os.path.join(input_dir, glob_pattern)))
    if not files:
        raise FileNotFoundError(f"No files matching '{glob_pattern}' found in {input_dir}")

    print(f"[pool_tissues] Found {len(files)} tissue file(s):")
    adatas, tissue_ids, marker_sets = [], [], {}
    for f in files:
        tissue_id = _tissue_id_from_filename(f)
        try:
            a = ad.read_h5ad(f)
        except Exception as e:
            print(f"[pool_tissues]   ⚠️  {os.path.basename(f)}: could not read ({e}) — skipping.")
            continue
        print(f"[pool_tissues]   {tissue_id}: {a.n_obs:,} cells x {a.n_vars} markers")
        adatas.append(a)
        tissue_ids.append(tissue_id)
        marker_sets[tissue_id] = set(a.var_names)

    if not adatas:
        raise RuntimeError("No tissue files could be loaded — nothing to pool.")

    if len(adatas) == 1:
        print(f"[pool_tissues]   ⚠️  Only 1 tissue file found — pooling is a no-op here, "
              f"but proceeding (output will just be a copy with tissue_id added).")

    all_markers = set.union(*marker_sets.values())
    common_markers = set.intersection(*marker_sets.values())
    if all_markers != common_markers:
        print(f"\n[pool_tissues] ⚠️  Marker panels differ across tissues:")
        for m in sorted(all_markers - common_markers):
            missing_from = [t for t in tissue_ids if m not in marker_sets[t]]
            print(f"[pool_tissues]      '{m}' missing from: {missing_from}")
        if join == "inner":
            print(f"[pool_tissues]   join=inner -> these markers will be DROPPED from the pooled "
                  f"file ({len(common_markers)}/{len(all_markers)} markers kept, common to all tissues).")
        else:
            print(f"[pool_tissues]   join=outer -> these markers will be KEPT but NaN-padded for "
                  f"tissues missing them — make sure downstream GMM fitting handles NaNs before "
                  f"using this file.")

    pooled = ad.concat(
        adatas, join=join, label="tissue_id", keys=tissue_ids, index_unique="-",
    )
    print(f"\n[pool_tissues] Pooled: {pooled.n_obs:,} cells x {pooled.n_vars} markers "
          f"from {len(adatas)} tissue(s)")

    out_dir = os.path.dirname(os.path.abspath(output_file))
    os.makedirs(out_dir, exist_ok=True)
    pooled.write_h5ad(output_file)
    print(f"[pool_tissues] ✓ Saved -> {output_file}")
    return pooled


def main():
    p = argparse.ArgumentParser(
        description="Pool multiple tissue-level combined h5ads (from preprocessing.py) into one "
                     "cohort-level h5ad for cell_typing.py's input_file."
    )
    p.add_argument("--input-dir", required=True,
                    help="Directory containing {tissue_id}_combined_all_experiment_groups.h5ad "
                         "files (preprocessing.py's combined_processed_data/individual_processed_data/).")
    p.add_argument("--output-file", required=True,
                    help="Path to write the pooled cohort-level h5ad.")
    p.add_argument("--glob-pattern", default="*_combined_all_experiment_groups.h5ad",
                    help="Filename pattern to match tissue files. "
                         "Default: '*_combined_all_experiment_groups.h5ad' (preprocessing.py's naming).")
    p.add_argument("--join", choices=["inner", "outer"], default="inner",
                    help="How to handle differing marker panels across tissues. 'inner' (default) "
                         "keeps only markers common to every tissue — safe for GMM fitting. 'outer' "
                         "keeps the union, NaN-padded — only use if you'll handle NaNs yourself.")
    args = p.parse_args()

    if not os.path.isdir(args.input_dir):
        print(f"[pool_tissues] ERROR: --input-dir does not exist: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    try:
        pool_tissues(args.input_dir, args.output_file, args.glob_pattern, args.join)
    except (FileNotFoundError, RuntimeError) as e:
        print(f"[pool_tissues] ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
