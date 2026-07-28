"""
spatia/analysis/tif_conversion.py
===================================
Config-driven conversion of specialty TIFF formats (QPTIFF, OME-TIFF) to
standard multi-page TIFF, with metadata preserved alongside as JSON.

Refactored from: 01_tif_conversion.ipynb

This is pipeline step 0 (upstream of preprocessing): raw scanner output
(QPTIFF / OME-TIFF) -> standard .tif + per-file metadata JSON, ready for
ROI masking (spatia.analysis.roi_masking) and then segmentation
(spatia.analysis.segmentation).

Entry point
-----------
run_tif_conversion(cfg) -- reads paths.raw_image_dir, writes to
paths.converted_tif_dir (defaults to <output_dir>/converted_tif/).

Accepts the parsed YAML config dict produced by:
    import yaml
    with open("config_example.yaml") as f:
        cfg = yaml.safe_load(f)
"""

from __future__ import annotations

import os
import json
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple

import tifffile


# ─────────────────────────────────────────────────────────────────────────────
# CORE CONVERSION FUNCTION
# (ported verbatim from 01_tif_conversion.ipynb -- no hardcoded paths)
# ─────────────────────────────────────────────────────────────────────────────

def convert_specialty_tiff_to_standard_tiff(
    input_dir: str,
    output_dir: str,
    image_format: str = ".qptiff",
    preserve_metadata: bool = True,
    overwrite: bool = False,
    verbose: bool = True,
) -> Tuple[List[str], Dict[str, str]]:
    """
    Convert specialty TIFF formats (QPTIFF, OME-TIFF, etc.) to standard TIFF
    while preserving metadata.

    Parameters
    ----------
    input_dir : str
        Directory containing input image files
    output_dir : str
        Directory where converted images will be saved
    image_format : str, optional
        File extension of input images, by default ".qptiff"
    preserve_metadata : bool, optional
        Whether to extract and save metadata separately, by default True
    overwrite : bool, optional
        Whether to overwrite existing files, by default False
    verbose : bool, optional
        Whether to print progress information, by default True

    Returns
    -------
    Tuple[List[str], Dict[str, str]]
        List of successfully converted files and dictionary of errors

    Notes
    -----
    Metadata is preserved in separate JSON files if preserve_metadata is True.
    A file is skipped if both the output .tif and its metadata JSON already
    exist (unless overwrite=True). Files starting with "._" (macOS resource
    forks) are excluded automatically.
    """
    os.makedirs(output_dir, exist_ok=True)
    metadata_dir = os.path.join(output_dir, "metadata")
    if preserve_metadata:
        os.makedirs(metadata_dir, exist_ok=True)

    target_files = [
        f for f in os.listdir(input_dir)
        if f.endswith(image_format) and not f.startswith("._")
    ]

    if not target_files:
        if verbose:
            print(f"No {image_format} files found in {input_dir}.")
        return [], {}

    if verbose:
        print(f"Found {len(target_files)} {image_format} files to convert")
        print("Starting conversion of files to .tif...")

    converted_files: List[str] = []
    errors: Dict[str, str] = {}

    for file_name in target_files:
        input_file = os.path.join(input_dir, file_name)
        base_name = file_name.replace(image_format, "")
        output_file = os.path.join(output_dir, base_name + ".tif")
        metadata_file = os.path.join(metadata_dir, base_name + "_metadata.json")

        if not overwrite and os.path.exists(output_file) and (
            not preserve_metadata or os.path.exists(metadata_file)
        ):
            if verbose:
                print(f"Skipping {file_name} as output already exists")
            continue

        if verbose:
            print(f"Converting: {file_name} -> {os.path.basename(output_file)}")

        try:
            with tifffile.TiffFile(input_file) as tif:
                image = tif.asarray()

                metadata: dict = {}
                if preserve_metadata:
                    if tif.is_ome:
                        ome_metadata = tif.ome_metadata
                        metadata["ome_xml"] = ome_metadata
                        try:
                            root = ET.fromstring(ome_metadata)
                            ns = "{http://www.openmicroscopy.org/Schemas/OME/2016-06}"
                            metadata["summary"] = {
                                "images": len(root.findall(f".//{ns}Image")),
                                "channels": len(root.findall(f".//{ns}Channel")),
                                "pixels_attributes": {
                                    k: v for k, v in
                                    root.find(f".//{ns}Pixels").attrib.items()
                                },
                            }
                        except Exception as e:
                            metadata["xml_parse_error"] = str(e)

                    metadata["tiff_tags"] = {}
                    for page_idx, page in enumerate(tif.pages):
                        page_metadata = {}
                        try:
                            for tag in page.tags.values():
                                if hasattr(tag, "name") and hasattr(tag, "value"):
                                    if not isinstance(tag.value, bytes):
                                        page_metadata[tag.name] = tag.value
                        except AttributeError:
                            try:
                                if hasattr(page, "keyframe") and hasattr(page.keyframe, "tags"):
                                    for tag in page.keyframe.tags.values():
                                        if not isinstance(tag.value, bytes):
                                            page_metadata[tag.name] = tag.value
                            except Exception as e:
                                page_metadata["error"] = f"Could not extract tags: {e}"
                        metadata["tiff_tags"][f"page_{page_idx}"] = page_metadata

            photometric = "rgb" if len(image.shape) > 2 and image.shape[-1] == 3 else "minisblack"

            if verbose:
                print(f"  Image shape: {image.shape}, dtype: {image.dtype}")
                print(f"  Size: {image.nbytes / (1024 * 1024):.2f} MB")

            tifffile.imwrite(
                output_file,
                image,
                photometric=photometric,
                metadata={"Simplified": f"Converted from {image_format}"},
                description=f"Converted from {file_name}",
                compression=None,
            )

            if preserve_metadata:
                with open(metadata_file, "w") as f:
                    f.write(json.dumps(metadata, default=str, indent=2))

            converted_files.append(output_file)

            if verbose:
                output_size = os.path.getsize(output_file) / (1024 * 1024)
                print(f"  Successfully converted: {os.path.basename(output_file)} ({output_size:.2f} MB)")
                if preserve_metadata:
                    print(f"  Metadata saved to: {os.path.basename(metadata_file)}")

        except Exception as e:
            errors[file_name] = str(e)
            if verbose:
                print(f"  Error converting {file_name}: {e}")

    if verbose:
        print(f"\nConversion complete. {len(converted_files)} files converted.")
        if errors:
            print(f"Encountered {len(errors)} errors during conversion.")

    return converted_files, errors


# ─────────────────────────────────────────────────────────────────────────────
# PIPELINE ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def run_tif_conversion(cfg: dict) -> dict:
    """
    Config-driven wrapper around convert_specialty_tiff_to_standard_tiff(),
    matching the run_pipeline.py step contract (accepts cfg, returns a
    result dict, raises on unrecoverable errors).

    Config keys
    -----------
    paths.raw_image_dir        -- required. Directory of raw QPTIFF/OME-TIFF files.
    paths.converted_tif_dir    -- optional. Defaults to <output_dir>/converted_tif/.
    tif_conversion.image_format       -- default ".qptiff"
    tif_conversion.preserve_metadata  -- default True
    tif_conversion.overwrite          -- default False
    tif_conversion.verbose            -- default True

    Returns
    -------
    dict with keys: converted_files, errors, output_dir
    """
    paths = cfg.get("paths", {})
    raw_dir = paths.get("raw_image_dir")
    if not raw_dir:
        raise KeyError(
            "paths.raw_image_dir is required for the tif_conversion step "
            "(directory containing raw .qptiff/.ome.tiff scanner output)."
        )
    if not os.path.isdir(raw_dir):
        raise FileNotFoundError(f"raw_image_dir not found: {raw_dir}")

    out_dir = paths.get("converted_tif_dir") or os.path.join(
        paths.get("output_dir", "."), "converted_tif"
    )

    tc_cfg = cfg.get("tif_conversion", {})
    image_format = tc_cfg.get("image_format", ".qptiff")
    preserve_metadata = tc_cfg.get("preserve_metadata", True)
    overwrite = tc_cfg.get("overwrite", False)
    verbose = tc_cfg.get("verbose", True)

    print("=" * 72)
    print("SPATIA TIF CONVERSION")
    print("=" * 72)
    print(f"Input dir  : {raw_dir}")
    print(f"Output dir : {out_dir}")
    print(f"Format     : {image_format}")

    converted_files, errors = convert_specialty_tiff_to_standard_tiff(
        input_dir=raw_dir,
        output_dir=out_dir,
        image_format=image_format,
        preserve_metadata=preserve_metadata,
        overwrite=overwrite,
        verbose=verbose,
    )

    if not converted_files and not errors:
        # No matching files at all -- likely a config problem, not silent success.
        raise FileNotFoundError(
            f"No '{image_format}' files found in {raw_dir}. "
            "Check paths.raw_image_dir and tif_conversion.image_format."
        )

    return {
        "converted_files": converted_files,
        "errors": errors,
        "output_dir": out_dir,
    }
