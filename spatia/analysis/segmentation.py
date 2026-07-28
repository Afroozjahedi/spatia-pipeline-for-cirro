"""
spatia/analysis/segmentation.py
==================================
Config-driven cell segmentation (Mesmer/Cellpose via spacec) on masked,
multichannel ROI TIFFs, with CSV export and QC overlay generation.

Refactored from: 03_segmentation.ipynb

Pipeline position: tif_conversion (step 0) -> roi_masking (step 1) ->
segmentation (step 2, this module) -> preprocessing (existing step 3).

IMPORTANT: this module's CSV output (*_mesmer_result.csv, one per slide
subfolder) is written to paths.segmentation_results_dir -- the exact
directory run_preprocessing() already expects as its input. This closes
the loop: the pipeline no longer starts "after preprocessing" -- it starts
here.

Known environment dependency
-----------------------------
Requires `spacec` (Mesmer or Cellpose backend). This is a heavy dependency
(pulls torch/tensorflow) -- see SPATIA_PIPELINE_LOG.md Day 2 for the
conda-forge-before-pip build note, and this project's later entries for
confirmation that a plain `pip install spacec` exhausts disk space in a
constrained sandbox. Import is module-level (matching preprocessing.py's
existing convention) so failures surface immediately and clearly rather
than partway through a run.

Entry point
-----------
run_segmentation(cfg) -- accepts the parsed YAML config dict.
"""

from __future__ import annotations

import os
import sys
import pickle
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage.measure import regionprops_table

import spacec as sp


# ─────────────────────────────────────────────────────────────────────────────
# PART 1 -- SEGMENTATION
# ─────────────────────────────────────────────────────────────────────────────

def _find_masked_tifs(processed_rois_dir: str) -> List[str]:
    """Recursively find masked ROI TIFFs, excluding preview images."""
    tif_files = []
    for root, _dirs, files in os.walk(processed_rois_dir):
        for f in files:
            if f.endswith(".tif") and "preview" not in f:
                tif_files.append(os.path.join(root, f))
    return tif_files


def run_cell_segmentation(
    processed_rois_dir: str,
    channel_file_path: str,
    output_dir: str,
    seg_method: str = "mesmer",
    nuclei_channel: str = "DAPI",
    membrane_channel_list: List[str] = ["CD45"],
    compartment: str = "whole-cell",
    input_format: str = "Multichannel",
    resize_factor: int = 1,
    size_cutoff: int = 0,
) -> Dict[str, dict]:
    """
    Run spacec cell segmentation over every masked ROI TIFF found under
    processed_rois_dir. Idempotent: existing *_seg_output.pickle files are
    skipped. Ported from 03_segmentation.ipynb Part 1.
    """
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(channel_file_path):
        raise FileNotFoundError(
            f"Channel file not found: {channel_file_path}. "
            "Segmentation cannot proceed without a channel-name mapping."
        )

    tif_files = _find_masked_tifs(processed_rois_dir)
    if not tif_files:
        raise FileNotFoundError(
            f"No masked ROI .tif files found in {processed_rois_dir} "
            "(or its subdirectories). Run roi_masking first."
        )
    print(f"Found {len(tif_files)} .tif files to process.")

    segmentation_outputs: Dict[str, dict] = {}
    print("Starting cell segmentation...")

    for input_file in tif_files:
        filename = os.path.basename(input_file)
        output_fname = os.path.splitext(filename)[0]
        slide_id = os.path.basename(os.path.dirname(input_file))
        slide_output_dir = os.path.join(output_dir, slide_id)
        os.makedirs(slide_output_dir, exist_ok=True)

        pickle_file = os.path.join(slide_output_dir, f"{output_fname}_seg_output.pickle")
        if os.path.exists(pickle_file):
            print(f"Segmentation output already exists, skipping: {pickle_file}")
            continue

        print(f"Segmenting: {input_file}")
        try:
            seg_output = sp.tl.cell_segmentation(
                file_name=input_file,
                channel_file=channel_file_path,
                output_dir=slide_output_dir,
                seg_method=seg_method,
                nuclei_channel=nuclei_channel,
                output_fname=output_fname,
                membrane_channel_list=membrane_channel_list,
                compartment=compartment,
                input_format=input_format,
                resize_factor=resize_factor,
                size_cutoff=size_cutoff,
            )
            segmentation_outputs[output_fname] = seg_output
            with open(pickle_file, "wb") as f:
                pickle.dump(seg_output, f)
            print(f"Segmentation completed for: {filename}")
        except Exception as e:
            print(f"Error during segmentation for {filename}: {e}")
            print(f"Stack trace: {sys.exc_info()}")

    print("All segmentation processes completed!")
    return segmentation_outputs


# ─────────────────────────────────────────────────────────────────────────────
# PART 2 -- EXPORT SEGMENTATION RESULTS TO *_mesmer_result.csv
# (this is the exact filename/shape spatia.analysis.preprocessing expects)
# ─────────────────────────────────────────────────────────────────────────────

def export_segmentation_to_csv(output_dir: str) -> List[str]:
    """
    Convert every *_seg_output.pickle under output_dir into a per-image
    *_mesmer_result.csv (label, centroid_x/y, area, eccentricity, perimeter,
    axis lengths, plus per-channel mean intensity) -- the exact input
    contract spatia.analysis.preprocessing.run_preprocessing() scans for.

    Ported from 03_segmentation.ipynb's CSV-export cell.
    """
    written_csvs: List[str] = []
    print("\nExporting segmentation results to *_mesmer_result.csv...")

    for slide_folder in os.listdir(output_dir):
        slide_dir = os.path.join(output_dir, slide_folder)
        if not os.path.isdir(slide_dir) or slide_folder == "csv_exports":
            continue

        pickle_files = [f for f in os.listdir(slide_dir) if f.endswith("_seg_output.pickle")]
        for pickle_file in pickle_files:
            pickle_path = os.path.join(slide_dir, pickle_file)
            base_name = pickle_file.replace("_seg_output.pickle", "")
            csv_path = os.path.join(slide_dir, f"{base_name}_mesmer_result.csv")
            if os.path.exists(csv_path):
                print(f"  ⏭️  Already exported: {csv_path}")
                written_csvs.append(csv_path)
                continue

            try:
                with open(pickle_path, "rb") as f:
                    seg_data = pickle.load(f)

                mask = seg_data["masks"]
                image_dict = seg_data.get("image_dict", {})

                props = regionprops_table(
                    mask,
                    properties=[
                        "label", "centroid", "area", "eccentricity",
                        "perimeter", "major_axis_length", "minor_axis_length",
                    ],
                )
                df = pd.DataFrame(props).rename(columns={
                    "centroid-0": "y", "centroid-1": "x",
                    "major_axis_length": "axis_major_length",
                    "minor_axis_length": "axis_minor_length",
                })

                if image_dict:
                    for ch_name, ch_data in image_dict.items():
                        if isinstance(ch_data, np.ndarray) and ch_data.shape == mask.shape:
                            try:
                                ch_props = regionprops_table(
                                    mask, intensity_image=ch_data.astype(float),
                                    properties=["label", "intensity_mean"],
                                )
                                df[ch_name] = ch_props["intensity_mean"]
                            except Exception as e:
                                print(f"    ⚠️  Failed to extract {ch_name}: {e}")

                df.insert(0, "image_ID", base_name)
                df.insert(1, "slide_folder", slide_folder)

                df.to_csv(csv_path, index=False)
                written_csvs.append(csv_path)
                print(f"  ✓ Exported {len(df)} cells → {os.path.basename(csv_path)}")
            except Exception as e:
                print(f"  ❌ Failed to process {pickle_file}: {e}")
                continue

    print(f"✓ {len(written_csvs)} *_mesmer_result.csv file(s) written under {output_dir}")
    return written_csvs


# ─────────────────────────────────────────────────────────────────────────────
# PART 3 -- QC OVERLAYS
# ─────────────────────────────────────────────────────────────────────────────

def generate_overlays(
    output_dir: str,
    nucleus_channel: str = "DAPI",
    additional_channels: List[str] = ["CD45"],
    show_subsample: bool = True,
    n: int = 2,
    tilesize: int = 400,
    rand_seed: int = 4,
) -> List[str]:
    """Generate <roi>_overlay.png QC images from *_seg_output.pickle files."""
    saved: List[str] = []
    slide_folders = [f for f in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, f))]
    if not slide_folders:
        print(f"No slide folders found in {output_dir}")
        return saved

    for slide_folder in slide_folders:
        slide_dir = os.path.join(output_dir, slide_folder)
        pickle_files = [f for f in os.listdir(slide_dir) if f.endswith("_seg_output.pickle")]
        for pickle_file in pickle_files:
            pickle_path = os.path.join(slide_dir, pickle_file)
            base_name = pickle_file.replace("_seg_output.pickle", "")
            save_path = os.path.join(slide_dir, f"{base_name}_overlay.png")
            if os.path.exists(save_path):
                saved.append(save_path)
                continue

            try:
                with open(pickle_path, "rb") as f:
                    seg_output = pickle.load(f)
            except Exception as e:
                print(f"  Error loading pickle file: {e}")
                continue

            try:
                plt.figure(figsize=(16, 12))
                sp.pl.show_masks(
                    seg_output=seg_output,
                    nucleus_channel=nucleus_channel,
                    additional_channels=additional_channels,
                    show_subsample=show_subsample,
                    n=n,
                    tilesize=tilesize,
                    rand_seed=rand_seed,
                )
                plt.savefig(save_path, dpi=150, bbox_inches="tight")
                plt.close()
                saved.append(save_path)
                print(f"  Overlay saved to: {save_path}")
            except Exception as e:
                print(f"  Error generating overlay: {e}")
                plt.close()
                continue

    return saved


# ─────────────────────────────────────────────────────────────────────────────
# PIPELINE ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def run_segmentation(cfg: dict) -> dict:
    """
    Config-driven wrapper matching the run_pipeline.py step contract.

    Config keys
    -----------
    paths.masked_roi_dir             -- input. Defaults to <output_dir>/masked_rois/
                                         (roi_masking step's output).
    paths.channel_file               -- required. Path to channelnames.txt.
    paths.segmentation_results_dir   -- output. This is the SAME directory
                                         preprocessing.run_preprocessing() reads
                                         from -- no separate wiring needed.
    segmentation.seg_method                default "mesmer"
    segmentation.nuclei_channel            default "DAPI"
    segmentation.membrane_channel_list     default ["CD45"]
    segmentation.compartment               default "whole-cell"
    segmentation.input_format              default "Multichannel"
    segmentation.resize_factor             default 1
    segmentation.size_cutoff               default 0
    segmentation.generate_overlays         default True

    Returns
    -------
    dict with keys: segmentation_outputs (dict), csv_files (list),
    overlay_files (list), output_dir (path)
    """
    paths = cfg.get("paths", {})
    processed_rois_dir = paths.get("masked_roi_dir") or os.path.join(
        paths.get("output_dir", "."), "masked_rois"
    )
    channel_file_path = paths.get("channel_file")
    if not channel_file_path:
        raise KeyError(
            "paths.channel_file is required for the segmentation step "
            "(path to channelnames.txt matching the marker panel)."
        )

    output_dir = paths.get("segmentation_results_dir")
    if not output_dir:
        raise KeyError(
            "paths.segmentation_results_dir is required -- this is the same "
            "directory spatia.analysis.preprocessing reads from."
        )

    seg_cfg = cfg.get("segmentation", {})
    seg_method = seg_cfg.get("seg_method", "mesmer")
    nuclei_channel = seg_cfg.get("nuclei_channel", "DAPI")
    membrane_channel_list = seg_cfg.get("membrane_channel_list", ["CD45"])
    compartment = seg_cfg.get("compartment", "whole-cell")
    input_format = seg_cfg.get("input_format", "Multichannel")
    resize_factor = seg_cfg.get("resize_factor", 1)
    size_cutoff = seg_cfg.get("size_cutoff", 0)
    do_overlays = seg_cfg.get("generate_overlays", True)

    print("=" * 72)
    print("SPATIA CELL SEGMENTATION")
    print("=" * 72)
    print(f"Input (masked ROIs) : {processed_rois_dir}")
    print(f"Channel file        : {channel_file_path}")
    print(f"Output dir          : {output_dir}")
    print(f"Method              : {seg_method}")

    if not os.path.isdir(processed_rois_dir):
        raise FileNotFoundError(f"masked_roi_dir not found: {processed_rois_dir}")

    segmentation_outputs = run_cell_segmentation(
        processed_rois_dir=processed_rois_dir,
        channel_file_path=channel_file_path,
        output_dir=output_dir,
        seg_method=seg_method,
        nuclei_channel=nuclei_channel,
        membrane_channel_list=membrane_channel_list,
        compartment=compartment,
        input_format=input_format,
        resize_factor=resize_factor,
        size_cutoff=size_cutoff,
    )

    csv_files = export_segmentation_to_csv(output_dir)

    overlay_files: List[str] = []
    if do_overlays:
        overlay_files = generate_overlays(output_dir, nucleus_channel=nuclei_channel,
                                           additional_channels=membrane_channel_list)

    return {
        "segmentation_outputs": segmentation_outputs,
        "csv_files": csv_files,
        "overlay_files": overlay_files,
        "output_dir": output_dir,
    }
