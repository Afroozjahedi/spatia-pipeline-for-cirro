#!/usr/bin/env python3
"""
build_spatialdata.py
====================
Packages one TMA core's SPATIA outputs as a SpatialData store, for
spatial-data-studio's cli.py to turn into a viewer checkpoint.

Inputs, all written by run_pipeline.py:
  --seg-pickle   segmentation_results/<slide>/<core>_seg_output.pickle
                 (Mesmer label mask + every image channel)
  --cells        the core's typed cell table: <group>_<core>_cells_with_triad_flags.csv
                 when the core has triads, else matched_cells/<group>_<core>_matched_with_boundaries.csv
                 (triads.py writes no flags file for a core with zero triads)
  --triad-pairs  <group>_<core>_triad_pairs.csv, when the core has triads

The store holds:
  images/image            every antibody channel plus one nuclear stain, pyramided
  labels/cell_labels      the Mesmer mask
  shapes/cell_boundaries  one polygon per cell, traced from the mask
  shapes/annotations      one triangle per triad (anchor, partner1, partner2)
  tables/table            one row per typed cell, annotating cell_boundaries by mask label

shapes/annotations follows spatial-data-studio's shape-annotation layout
(backend/app/sessions/shape_annotations.py there), which its viewer draws over
the image and cells.

Usage:
    python bin/build_spatialdata.py \
        --seg-pickle results/segmentation_results/TMA_B/reg001_X01_Y01_Z08_seg_output.pickle \
        --cells results/CLR_reg001_X01_Y01_Z08_cells_with_triad_flags.csv \
        --triad-pairs results/CLR_reg001_X01_Y01_Z08_triad_pairs.csv \
        --output reg001_X01_Y01_Z08.zarr
"""
import argparse
import json
import pickle
import re

import anndata as ad
import geopandas as gpd
import numpy as np
import pandas as pd
import spatialdata as sd
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
from skimage.measure import regionprops_table
from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel, TableModel

SHAPES_KEY = "cell_boundaries"
LABEL_KEY = "cell_label"
# The element name spatial-data-studio reads shape annotations from.
ANNOTATIONS_KEY = "annotations"
# One outline colour per triad combination (anchor/partner1/partner2 types),
# assigned in sorted order of the combination label.
TRIAD_COLORS = ["#ffd700", "#00e5ff", "#ff4fd8", "#7cff4f", "#ff8c00"]
TRIAD_FILL_ALPHA = 0.15
SCALE_FACTORS = [2, 2, 2]
# Acquisition slots with no antibody in the CODEX panel, plus the composite
# segmentation.py adds to image_dict (it duplicates the pickle's `img`).
SKIPPED_CHANNEL = re.compile(r"^(blank|empty)|^segmentation_channel$", re.IGNORECASE)
# CODEX re-images the nuclear stain every cycle (HOECHST1..23, one spelled
# HOCHST13) and adds DRAQ5 at the end; only the first one is kept.
NUCLEAR_CHANNEL = re.compile(r"^HOE?CHST\d+ |^DRAQ5 ", re.IGNORECASE)
TRIAD_FLAGS = ["is_triad_anchor", "is_triad_partner1", "is_triad_partner2", "in_any_triad"]
# A typed cell's centroid is the Mesmer centroid after a CSV round trip, so it
# matches the mask's own centroid to float precision, not exactly.
CENTROID_TOLERANCE_PX = 1e-3
# Per-marker GMM outputs from cell_typing.py; the *_pos calls are kept as obs,
# these are not.
DROPPED_COLUMN = re.compile(r"_(posterior|intensity)$|^Unnamed: 0$")


class _NumpyOnlyUnpickler(pickle.Unpickler):
    """The seg pickle is a dict of numpy arrays; refuse anything else it names."""

    def find_class(self, module, name):
        # Pickled under numpy 1.x, read here under numpy 2.x.
        if module == "numpy.core.multiarray":
            module = "numpy._core.multiarray"
        if module not in ("numpy", "numpy._core.multiarray"):
            raise pickle.UnpicklingError(f"refusing to load {module}.{name} from a seg pickle")
        return super().find_class(module, name)


def load_segmentation(path):
    with open(path, "rb") as fh:
        seg = _NumpyOnlyUnpickler(fh).load()
    masks = np.squeeze(seg["masks"])
    names = [name for name in seg["image_dict"] if not SKIPPED_CHANNEL.search(name)]
    nuclear = [name for name in names if NUCLEAR_CHANNEL.search(name)]
    channels = {name: seg["image_dict"][name] for name in names if name not in nuclear[1:]}
    return masks, channels


def match_cells_to_labels(cells, masks):
    """Mask label of each typed cell, matched on centroid."""
    props = pd.DataFrame(regionprops_table(masks, properties=("label", "centroid")))
    dist, idx = cKDTree(props[["centroid-1", "centroid-0"]].to_numpy()).query(
        cells[["centroid_x", "centroid_y"]].to_numpy()
    )
    unmatched = int((dist > CENTROID_TOLERANCE_PX).sum())
    if unmatched:
        raise ValueError(
            f"{unmatched} of {len(cells)} cells have no mask label within {CENTROID_TOLERANCE_PX} px "
            "of their centroid -- the cell table and seg pickle are not from the same segmentation"
        )
    labels = props["label"].to_numpy()[idx]
    if len(np.unique(labels)) != len(labels):
        raise ValueError("two cells matched the same mask label")
    return labels


def build_table(cells, labels, channel_names):
    marker_cols = [c for c in cells.columns if c in channel_names]
    obs_cols = [
        c for c in cells.columns
        if c not in marker_cols and c not in ("centroid_x", "centroid_y")
        and not DROPPED_COLUMN.search(c) and not NUCLEAR_CHANNEL.search(c)
    ]
    obs = cells[obs_cols].copy()
    for flag in TRIAD_FLAGS:
        if flag not in obs:
            obs[flag] = False
    obs["cell_type"] = obs["cell_type"].astype("category")
    obs["experiment_group"] = obs["experiment_group"].astype("category")
    obs["region"] = pd.Categorical([SHAPES_KEY] * len(obs))
    obs[LABEL_KEY] = labels
    obs.index = labels.astype(str)

    table = ad.AnnData(
        X=cells[marker_cols].to_numpy(dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=marker_cols),
    )
    table.obsm["spatial"] = cells[["centroid_x", "centroid_y"]].to_numpy()
    return TableModel.parse(table, region=SHAPES_KEY, region_key="region", instance_key=LABEL_KEY)


def build_triad_annotations(pairs):
    """One filled triangle per triad, as spatial-data-studio annotation rows."""
    combos = sorted(pairs["triad_combo_label"].unique())
    color = {combo: TRIAD_COLORS[i % len(TRIAD_COLORS)] for i, combo in enumerate(combos)}
    rows, triangles = [], []
    for row in pairs.itertuples():
        vertices = [
            [row.anchor_x, row.anchor_y],
            [row.partner1_x, row.partner1_y],
            [row.partner2_x, row.partner2_y],
        ]
        triangles.append(Polygon(vertices))
        rows.append({
            "kind": "polygon",
            "params": json.dumps({"vertices": vertices}),
            "label": f"{row.anchor_type} / {row.partner1_type} / {row.partner2_type}",
            "stroke_color": color[row.triad_combo_label],
            "stroke_width": 2.0,
            "stroke_dash": "solid",
            "stroke_arrow_start": False,
            "stroke_arrow_end": False,
            "stroke_arrow_size": 10.0,
            "stroke_z": 0,
            "fill_enabled": True,
            "fill_color": color[row.triad_combo_label],
            "fill_alpha": TRIAD_FILL_ALPHA,
            "fill_z": 0,
        })
    # Numbered like triads.py's QuPath export (triad_00001, ...).
    index = [f"triad_{i:05d}" for i in range(1, len(rows) + 1)]
    return ShapesModel.parse(gpd.GeoDataFrame(rows, geometry=triangles, index=index))


def build_spatialdata(seg_pickle, cells_csv, triad_pairs_csv=None):
    masks, channels = load_segmentation(seg_pickle)
    cells = pd.read_csv(cells_csv)
    labels = match_cells_to_labels(cells, masks)

    image = Image2DModel.parse(
        np.stack(list(channels.values())),
        dims=("c", "y", "x"),
        c_coords=list(channels),
        scale_factors=SCALE_FACTORS,
    )
    cell_labels = Labels2DModel.parse(masks, dims=("y", "x"), scale_factors=SCALE_FACTORS)
    # Only cells that survived QC have a row in the table; outline just those.
    polygons = sd.to_polygons(Labels2DModel.parse(np.where(np.isin(masks, labels), masks, 0), dims=("y", "x")))
    polygons.index.name = None
    shapes = {SHAPES_KEY: ShapesModel.parse(polygons)}
    if triad_pairs_csv:
        shapes[ANNOTATIONS_KEY] = build_triad_annotations(pd.read_csv(triad_pairs_csv))

    return sd.SpatialData(
        images={"image": image},
        labels={"cell_labels": cell_labels},
        shapes=shapes,
        tables={"table": build_table(cells, labels, set(channels))},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seg-pickle", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--triad-pairs", help="the core's *_triad_pairs.csv; omit for a core with no triads")
    parser.add_argument("--output", required=True, help="path of the .zarr store to write")
    args = parser.parse_args()

    sdata = build_spatialdata(args.seg_pickle, args.cells, args.triad_pairs)
    sdata.write(args.output)
    print(sdata)


if __name__ == "__main__":
    main()
