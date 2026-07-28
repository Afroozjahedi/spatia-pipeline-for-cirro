"""
spatia/analysis/roi_masking.py
=================================
Config-driven extraction and masking of ROI sub-images from whole-slide /
whole-core multiplexed TIFFs, using masks exported from QuPath (either the
existing manual ROI workflow, 00_ROI_extract_mask_project.groovy, or the
automatic TMA-dearray script, qupath_scripts/00b_auto_tma_dearray.groovy).

Refactored from: 02_masked_tif.ipynb

Pipeline position: tif_conversion (step 0) -> roi_masking (step 1) ->
segmentation (step 2) -> preprocessing (existing step 3, unchanged).

Directory contract
-------------------
Input:
    converted_tif_dir/
        *.tif / *.tiff                 <- from tif_conversion
        metadata/{basename}.json       <- optional channel-name metadata
    qupath_roi_dir/
        {slide_id}_binary_mask.png     <- whole-tissue binary mask
        {slide_id}_roi_labels.txt      <- ROI coordinate table (tab-separated)
        individual_masks/{slide_id}/{condition}_*.png

Output:
    masked_roi_dir/{slide_id}/
        {condition}_{slide_id}_x{x}_y{y}_w{w}_h{h}_masked.tif
        {condition}_{slide_id}_x{x}_y{y}_preview.png
        {condition}_{slide_id}_x{x}_y{y}_channel_info.txt
    masked_roi_dir/roi_processing_summary.csv

Entry point
-----------
run_roi_masking(cfg) -- accepts the parsed YAML config dict.

Note: the original notebook prompted interactively for condition prefixes.
That prompt is removed here -- conditions come from experiment.conditions
in the config, so this step can run headlessly (required for Cirro).
"""

from __future__ import annotations

import os
import glob
import re
import json
import gc
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import tifffile
import imageio.v2 as imageio
from skimage import transform, io
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from PIL import Image
Image.MAX_IMAGE_PIXELS = None  # allow large whole-slide/whole-core images


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS (ported from 02_masked_tif.ipynb)
# ─────────────────────────────────────────────────────────────────────────────

def find_tiff_files(directory: str) -> List[str]:
    tiff_files = []
    for ext in [".tif", ".tiff"]:
        tiff_files.extend(glob.glob(os.path.join(directory, f"*{ext}")))
    return tiff_files


def load_channel_names_from_json(tiff_path: str, metadata_dir: str) -> Optional[List[str]]:
    """Load channel names from a JSON metadata file, trying several common schemas."""
    tiff_basename = os.path.splitext(os.path.basename(tiff_path))[0]
    potential_json_files = [
        os.path.join(metadata_dir, f"{tiff_basename}_metadata.json"),
        os.path.join(metadata_dir, f"{tiff_basename}.json"),
        os.path.join(metadata_dir, f"metadata_{tiff_basename}.json"),
    ]
    json_path = next((p for p in potential_json_files if os.path.exists(p)), None)
    if not json_path:
        print(f"No JSON metadata file found for {tiff_basename} in {metadata_dir}")
        return None

    try:
        with open(json_path, "r") as f:
            metadata = json.load(f)

        channel_names = None
        if "channel_names" in metadata:
            channel_names = metadata["channel_names"]
        elif "channels" in metadata and isinstance(metadata["channels"], list):
            channel_names = [ch.get("name", ch.get("Name", f"Channel_{i}"))
                              for i, ch in enumerate(metadata["channels"])]
        elif "image" in metadata and "channels" in metadata.get("image", {}):
            channels = metadata["image"]["channels"]
            if isinstance(channels, list):
                channel_names = [ch.get("name", ch.get("Name", f"Channel_{i}"))
                                  for i, ch in enumerate(channels)]
        elif "ome" in metadata:
            ome = metadata["ome"]
            pixels = ome.get("Image", {}).get("Pixels", {})
            channels = pixels.get("Channel")
            if isinstance(channels, list):
                channel_names = [ch.get("Name", f"Channel_{i}") for i, ch in enumerate(channels)]
        elif "panel" in metadata:
            panel = metadata["panel"]
            if isinstance(panel, list):
                channel_names = panel
            elif isinstance(panel, dict) and "targets" in panel:
                channel_names = panel["targets"]
        elif "markers" in metadata and isinstance(metadata["markers"], list):
            channel_names = metadata["markers"]
        elif "ome_xml" in metadata:
            try:
                import xml.etree.ElementTree as ET
                xml_string = metadata["ome_xml"]
                if isinstance(xml_string, str):
                    xml_string = xml_string.encode("utf-8")
                root = ET.fromstring(xml_string)
                ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
                channels = root.findall(".//ome:Channel", ns) or root.findall(".//{*}Channel")
                channel_names = [ch.attrib.get("Name", f"Channel_{i}") for i, ch in enumerate(channels)]
            except Exception as e:
                print(f"Error parsing OME-XML from JSON: {e}")
                channel_names = None

        if channel_names:
            print(f"Loaded {len(channel_names)} channel names from JSON: {os.path.basename(json_path)}")
            return channel_names
        print(f"Could not extract channel names from JSON structure in {os.path.basename(json_path)}")
        return None
    except Exception as e:
        print(f"Error loading JSON metadata from {json_path}: {e}")
        return None


def get_slide_id(filename: str, conditions: List[str]) -> str:
    """Strip condition prefix + coordinate suffix from a TIFF filename to get a stable slide ID."""
    basename = os.path.splitext(os.path.basename(filename))[0]
    for condition in conditions:
        if basename.startswith(f"{condition}_"):
            parts = basename.split("_")
            prefix = parts[0] + "_"
            for part in parts[1:]:
                if part.startswith(("x", "y")) and len(part) > 1 and part[1:].isdigit():
                    return basename[len(prefix):basename.find(f"_{part}")]
    coord_match = re.search(r"_x\d+_y\d+", basename)
    if coord_match:
        return basename[:coord_match.start()]
    return basename


def find_related_files(slide_id: str, masks_dir: str) -> Dict[str, Optional[str]]:
    """Locate the binary mask, ROI label table, and individual masks dir for a slide."""
    binary_mask_path = os.path.join(masks_dir, f"{slide_id}_binary_mask.png")
    if not os.path.exists(binary_mask_path):
        candidates = glob.glob(os.path.join(masks_dir, f"*{slide_id}*binary*mask*.png"))
        binary_mask_path = candidates[0] if candidates else None

    roi_labels_path = os.path.join(masks_dir, f"{slide_id}_roi_labels.txt")
    if not os.path.exists(roi_labels_path):
        candidates = glob.glob(os.path.join(masks_dir, f"*{slide_id}*roi*labels*.txt"))
        roi_labels_path = candidates[0] if candidates else None

    individual_masks_dir = os.path.join(masks_dir, "individual_masks", slide_id)
    if not os.path.exists(individual_masks_dir):
        individual_masks_dir = None
        for root, _dirs, _files in os.walk(masks_dir):
            if slide_id in root and "individual" in root.lower():
                individual_masks_dir = root
                break

    return {
        "binary_mask": binary_mask_path,
        "roi_labels": roi_labels_path,
        "individual_masks_dir": individual_masks_dir,
    }


def is_tissue_processed(output_dir: str, slide_id: str) -> bool:
    slide_output_dir = os.path.join(output_dir, slide_id)
    return (
        len(glob.glob(os.path.join(slide_output_dir, "*_masked.tif"))) > 0
        and len(glob.glob(os.path.join(slide_output_dir, "*_preview.png"))) > 0
    )


def extract_roi_coordinates(roi_name: str, roi_df: Optional[pd.DataFrame] = None) -> Optional[Tuple[int, int, int, int]]:
    """Extract (x, y, w, h) from an ROI name, falling back to the ROI label table."""
    match = re.search(r"x(\d+)_y(\d+)_w(\d+)_h(\d+)", roi_name)
    if match:
        return tuple(int(g) for g in match.groups())  # type: ignore[return-value]

    if roi_df is not None:
        roi_info = roi_df[roi_df["Name"] == roi_name]
        if not roi_info.empty:
            cols = roi_df.columns
            x_col = next((c for c in cols if c.lower() == "x"), None)
            y_col = next((c for c in cols if c.lower() == "y"), None)
            w_col = next((c for c in cols if "width" in c.lower()), None)
            h_col = next((c for c in cols if "height" in c.lower()), None)
            if all([x_col, y_col, w_col, h_col]):
                return (
                    int(roi_info[x_col].values[0]), int(roi_info[y_col].values[0]),
                    int(roi_info[w_col].values[0]), int(roi_info[h_col].values[0]),
                )
    return None


def process_wsi(
    wsi_path: str,
    binary_mask_path: str,
    roi_labels_path: Optional[str],
    individual_masks_dir: str,
    output_dir: str,
    metadata_dir: str,
    conditions: List[str],
) -> List[dict]:
    """Crop + mask every ROI in one whole-slide/whole-core TIFF. Returns one dict per ROI."""
    slide_id = get_slide_id(wsi_path, conditions)
    print(f"\n===== Processing {slide_id} =====")

    if is_tissue_processed(output_dir, slide_id):
        print(f"Skipping {slide_id}: Processed ROI files already exist")
        return []
    if not binary_mask_path or not os.path.exists(binary_mask_path):
        print(f"ERROR: Binary mask not found for {slide_id}")
        return []
    if not individual_masks_dir or not os.path.exists(individual_masks_dir):
        print(f"ERROR: Individual masks directory not found for {slide_id}")
        return []

    slide_output_dir = os.path.join(output_dir, slide_id)
    os.makedirs(slide_output_dir, exist_ok=True)

    print(f"Loading binary mask: {os.path.basename(binary_mask_path)}...")
    try:
        binary_mask = imageio.imread(binary_mask_path)
        binary_mask = (binary_mask > 0).astype(np.uint8)
    except Exception as e:
        print(f"Failed to load binary mask: {e}")
        return []

    roi_df = None
    if roi_labels_path and os.path.exists(roi_labels_path):
        try:
            roi_df = pd.read_csv(roi_labels_path, sep="\t")
            print(f"Loaded ROI labels, found {len(roi_df)} ROIs")
        except Exception as e:
            print(f"Error loading ROI labels: {e}")
    else:
        print("No ROI labels file found, will use filename coordinates only")

    channel_names = load_channel_names_from_json(wsi_path, metadata_dir)

    try:
        with tifffile.TiffFile(wsi_path) as tif:
            image_shape = tif.series[0].shape
            num_channels = image_shape[0] if len(image_shape) > 2 else 1
            height, width = image_shape[-2], image_shape[-1]
            print(f"WSI dimensions: {image_shape}, {num_channels} channels")

            if not channel_names and num_channels > 1 and tif.ome_metadata is not None:
                try:
                    import xml.etree.ElementTree as ET
                    root = ET.fromstring(tif.ome_metadata)
                    ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
                    image_element = root.find(".//ome:Image", ns)
                    if image_element is not None:
                        channels = image_element.findall(".//ome:Channel", ns)
                        channel_names = [c.attrib.get("Name", f"Channel_{i}") for i, c in enumerate(channels)]
                    if len(channel_names or []) != num_channels:
                        channel_names = [f"Channel_{i}" for i in range(num_channels)]
                except Exception as e:
                    print(f"Warning: Could not parse OME-XML for channel names: {e}")

            if not channel_names or len(channel_names) != num_channels:
                channel_names = [f"Channel_{i}" for i in range(num_channels)]
                print(f"Using default channel names: {channel_names}")

            if binary_mask.shape != (height, width):
                print(f"Resizing binary mask from {binary_mask.shape} to ({height}, {width})")
                binary_mask = transform.resize(binary_mask, (height, width), preserve_range=True).astype(np.uint8)
    except Exception as e:
        print(f"Error analyzing TIFF file: {e}")
        return []

    roi_names = []
    for pattern in [f"{condition}_*.png" for condition in conditions]:
        roi_names.extend(glob.glob(os.path.join(individual_masks_dir, pattern)))
    roi_names = [os.path.basename(n) for n in roi_names]
    print(f"Found {len(roi_names)} ROI masks for conditions: {', '.join(conditions)}")

    processed_rois = []
    for i, roi_name in enumerate(roi_names):
        if i > 0 and i % 5 == 0:
            gc.collect()
            plt.close("all")

        roi_condition = next((c for c in conditions if roi_name.startswith(f"{c}_")), "Unknown")
        coords = extract_roi_coordinates(roi_name.replace(".png", ""), roi_df)
        if not coords:
            print(f"Could not extract coordinates for {roi_name}, skipping")
            continue

        x, y, w, h = coords
        print(f"Processing {roi_condition} ROI at x={x}, y={y}, w={w}, h={h}")
        x_end, y_end = min(x + w, width), min(y + h, height)
        x_width, y_height = x_end - x, y_end - y

        roi_mask_path = os.path.join(individual_masks_dir, roi_name)
        try:
            individual_roi_mask = io.imread(roi_mask_path)
            individual_roi_mask = (individual_roi_mask > 0).astype(np.uint8)
        except Exception as e:
            print(f"Error loading individual ROI mask: {e}")
            individual_roi_mask = None

        roi_output = np.zeros((num_channels, y_height, x_width), dtype=np.uint16)
        with tifffile.TiffFile(wsi_path) as tif:
            for c in range(num_channels):
                try:
                    if hasattr(tif.series[0], "levels") and len(tif.series[0].levels) > 0:
                        roi_data = tif.series[0].levels[0].asarray(key=c)[y:y_end, x:x_end]
                    else:
                        page_idx = c if num_channels > 1 else 0
                        roi_data = tif.pages[page_idx].asarray()[y:y_end, x:x_end]

                    mask_region = binary_mask[y:y_end, x:x_end]
                    if individual_roi_mask is not None:
                        if individual_roi_mask.shape != roi_data.shape:
                            ind_mask = transform.resize(
                                individual_roi_mask, roi_data.shape, preserve_range=True
                            ).astype(np.uint8)
                        else:
                            ind_mask = individual_roi_mask
                        masked_roi = roi_data * (mask_region > 0) * (ind_mask > 0)
                    else:
                        masked_roi = roi_data * (mask_region > 0)
                    roi_output[c, :roi_data.shape[0], :roi_data.shape[1]] = masked_roi
                except Exception as e:
                    print(f"Error processing channel {c}: {e}")

        output_filename = f"{roi_condition}_{slide_id}_x{x}_y{y}_w{w}_h{h}_masked.tif"
        output_path = os.path.join(slide_output_dir, output_filename)

        if len(channel_names) != num_channels:
            channel_names = [f"Channel_{i}" for i in range(num_channels)]

        ijmetadata = {
            "Info": "\n".join(f"Channel {i}: {name}" for i, name in enumerate(channel_names)),
            "Labels": channel_names,
            "channels": num_channels,
            "slices": 1,
            "frames": 1,
            "hyperstack": True,
            "mode": "composite",
        }
        tifffile.imwrite(
            output_path, roi_output, imagej=True, metadata=ijmetadata, photometric="minisblack",
        )

        channel_info_path = os.path.join(slide_output_dir, f"{roi_condition}_{slide_id}_x{x}_y{y}_channel_info.txt")
        with open(channel_info_path, "w") as f:
            f.write("Channel Names:\n")
            for idx, name in enumerate(channel_names):
                f.write(f"Channel {idx}: {name}\n")

        fig = plt.figure(figsize=(8, 8))
        plt.imshow(roi_output[0], cmap="viridis")
        plt.title(f"{roi_condition} ROI - {slide_id} - Channel 0\n{channel_names[0] if channel_names else ''}")
        plt.colorbar()
        plt.savefig(os.path.join(slide_output_dir, f"{roi_condition}_{slide_id}_x{x}_y{y}_preview.png"))
        plt.close(fig)

        print(f"Saved {num_channels}-channel TIFF: {output_filename}")
        processed_rois.append({
            "slide_id": slide_id, "roi_condition": roi_condition,
            "x": x, "y": y, "width": w, "height": h,
            "output_file": output_path,
            "channel_names": ",".join(channel_names) if channel_names else "",
        })

    return processed_rois


# ─────────────────────────────────────────────────────────────────────────────
# PIPELINE ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def run_roi_masking(cfg: dict) -> dict:
    """
    Config-driven wrapper matching the run_pipeline.py step contract.

    Config keys
    -----------
    experiment.conditions           -- reused from the existing config section
                                        (no separate/duplicated condition list)
    paths.converted_tif_dir         -- input. Defaults to <output_dir>/converted_tif/
                                        (tif_conversion step's output).
    paths.qupath_roi_dir            -- required. Output of either the manual
                                        (00_ROI_extract_mask_project.groovy) or
                                        automatic (00b_auto_tma_dearray.groovy)
                                        QuPath ROI/TMA export.
    paths.masked_roi_dir            -- output. Defaults to <output_dir>/masked_rois/.

    Returns
    -------
    dict with keys: processed_rois (list), condition_counts (dict),
    summary_csv (path), output_dir (path)
    """
    paths = cfg.get("paths", {})
    conditions = cfg["experiment"]["conditions"]

    tiff_dir = paths.get("converted_tif_dir") or os.path.join(
        paths.get("output_dir", "."), "converted_tif"
    )
    metadata_dir = os.path.join(tiff_dir, "metadata")
    masks_dir = paths.get("qupath_roi_dir")
    if not masks_dir:
        raise KeyError(
            "paths.qupath_roi_dir is required for the roi_masking step "
            "(output directory of the QuPath ROI/TMA export)."
        )
    if not os.path.isdir(masks_dir):
        raise FileNotFoundError(f"qupath_roi_dir not found: {masks_dir}")

    output_dir = paths.get("masked_roi_dir") or os.path.join(
        paths.get("output_dir", "."), "masked_rois"
    )
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 72)
    print("SPATIA ROI MASKING")
    print("=" * 72)
    print(f"TIFF dir      : {tiff_dir}")
    print(f"QuPath masks  : {masks_dir}")
    print(f"Output dir    : {output_dir}")
    print(f"Conditions    : {conditions}")

    tiff_files = find_tiff_files(tiff_dir)
    print(f"Found {len(tiff_files)} TIFF files")
    if not tiff_files:
        raise FileNotFoundError(
            f"No .tif/.tiff files found in {tiff_dir}. Run tif_conversion first, "
            "or check paths.converted_tif_dir."
        )

    all_processed_rois: List[dict] = []
    for tiff_file in tiff_files:
        slide_id = get_slide_id(tiff_file, conditions)
        print(f"\nProcessing WSI: {slide_id} ({os.path.basename(tiff_file)})")
        related = find_related_files(slide_id, masks_dir)
        if not related["binary_mask"]:
            print(f"ERROR: No binary mask found for {slide_id}, skipping")
            continue
        if not related["individual_masks_dir"]:
            print(f"ERROR: No individual masks directory found for {slide_id}, skipping")
            continue

        processed = process_wsi(
            tiff_file, related["binary_mask"], related["roi_labels"],
            related["individual_masks_dir"], output_dir, metadata_dir, conditions,
        )
        all_processed_rois.extend(processed)
        gc.collect()
        plt.close("all")

    condition_counts = {
        cond: sum(1 for roi in all_processed_rois if roi.get("roi_condition") == cond)
        for cond in conditions
    }

    print("\n===== ROI MASKING SUMMARY =====")
    print(f"Total WSI/core files scanned: {len(tiff_files)}")
    print(f"Total ROIs extracted: {len(all_processed_rois)}")
    for cond, count in condition_counts.items():
        print(f"  {cond}: {count}")

    summary_csv = None
    if all_processed_rois:
        summary_csv = os.path.join(output_dir, "roi_processing_summary.csv")
        pd.DataFrame(all_processed_rois).to_csv(summary_csv, index=False)
        print(f"Summary saved to: {summary_csv}")
    else:
        print(
            "⚠️  No ROIs were extracted. Check that qupath_roi_dir contains "
            "binary_mask.png / roi_labels.txt / individual_masks/ for these slide IDs."
        )

    return {
        "processed_rois": all_processed_rois,
        "condition_counts": condition_counts,
        "summary_csv": summary_csv,
        "output_dir": output_dir,
    }
