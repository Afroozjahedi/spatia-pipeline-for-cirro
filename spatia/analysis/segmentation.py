"""
spatia/analysis/segmentation.py
==================================
Config-driven cell segmentation (Mesmer/Cellpose via spacec) on masked,
multichannel ROI TIFFs, with CSV export and QC overlay generation.

CHANGELOG -- 2026-09-04 (Afrouz + Claude consulting session)
--------------------------------------------------------------
Context: the 2026-09-03 real-data run (pipeline_20260903_133016.log) against
all four CRC TMA folders (276 raw .tif files) failed completely -- 0 succeeded,
139 failed, 137 skipped for "channel mismatch". Root cause, confirmed against
the actual files on rsrch6 (not assumed): the raw cores are CODEX/PhenoCycler
ImageJ hyperstacks, axes "TCYX", shape (23, 4, 1440, 1920) -- 23 acquisition
cycles x 4 channels/cycle = 92 planes, NOT a flat 92-channel stack. Two bugs
followed from that one fact, both fixed below. Full reasoning and the pixel-
level verification (nuclear-periodicity correlation, r=0.87-0.92 across all 23
HOECHST cycles + DRAQ5) are in the new spatia/analysis/channels.py docstring
and in this session's chat log -- not repeated here in full.

1. `_find_masked_tifs()`: now skips filenames starting with "._" -- macOS
   AppleDouble resource-fork sidecars created when this data was copied through
   an HFS+/APFS volume. 134 of these were in the 276-file batch; each is ~4 KB
   of metadata, not an image, and crashed tifffile with
   "not a TIFF file b'\\x00\\x05\\x16\\x07'" the moment segmentation tried to
   open one. Pure noise removal -- does not change segmentation results for
   any real file.

2. `_check_channel_order()`: rewritten. The OLD check compared this image's
   embedded ImageJ 'Labels' against channel_file NAME-FOR-NAME. That is correct
   for ROIs written by roi_masking.py (which embeds real marker names) but can
   NEVER pass for raw CODEX acquisitions, whose embedded labels are scanner
   placeholders ("ch1".."ch4" repeating 23x) carrying no marker identity at
   all -- that's what rejected all 137 real cores as "mismatched" when nothing
   was actually wrong with them. The check now dispatches on what the image
   actually carries (via the new channels.StackInfo):
     - real names embedded            -> compare name-for-name, as before
     - generic/absent labels          -> compare PLANE COUNT instead (the
                                          comparison that actually catches a
                                          truncated acquisition or wrong panel)
     - OME-TIFF/QPTIFF whose names
       could not be parsed            -> FAIL LOUD, does not fall through to
                                          the generic path (see channels.py
                                          changelog -- this is the safety net
                                          against silently relabeling a file
                                          that already had real names)

3. `run_cell_segmentation()` / `run_segmentation()`: new normalization step.
   spacec's `sp.tl.cell_segmentation()` opens the file path itself and treats
   axis 0 as the channel axis -- on a (23, 4, Y, X) array that means it sees 23
   "channels", which is the direct cause of the 2026-09-03
   "IndexError: index 23 is out of bounds for axis 0 with size 23" on
   reg021_X01_Y01_Z08.tif (that file was NOT corrupt or short a cycle -- every
   one of the 137 real cores has this same shape and would have hit the same
   crash once past the name-mismatch check). Each hyperstack is now flattened
   cycle-major to a flat, honestly-labeled (92, Y, X) TIFF via
   channels.normalize_stack_to_file() before being handed to spacec. New
   params: `normalized_dir` (default "<output_dir>/_normalized"),
   `keep_normalized` (default False -- each ~508 MB normalized copy is deleted
   right after its image is segmented, so peak extra disk is one file, not a
   second copy of the ~72 GB four-folder dataset), `auto_panel_template`
   (writes a positional "cyc01_ch1"... placeholder panel if channel_file is
   missing, so the mechanics can be smoke-tested -- loudly flagged as NOT real
   marker identity).

4. `nuclei_channel` / `membrane_channel_list`: now resolved once per run via
   channels.resolve_channels() instead of passed to spacec as literal config
   strings. Lets the yaml say "auto" (-> channel 1 of cycle 1, the CODEX
   nuclear convention) and short names ("CD45" -> "CD45 - hematopoietic
   cells") instead of requiring an exact match to the panel file's full text.
   Ambiguous short names (e.g. "CD4", which prefix-matches CD44/CD45/CD45RA/
   CD45RO in this panel) raise rather than silently picking one.

5. One pixel-level QC check added per batch: channels.verify_nuclear_periodicity()
   confirms channel-1-of-every-cycle actually correlates across cycles (r>=0.5)
   after flattening -- the automated version of the manual correlation check
   run by hand during this session. If it fails, the run HALTS rather than
   proceeding on a possibly-wrong flatten order, since a wrong order would
   silently swap marker identities for every cell.

Net effect: this is a bug fix (CODEX hyperstacks were unsupported, and the
validation logic was checking something that could never be true for them),
not a behavior change for any input this pipeline had previously segmented
successfully. channel_check=True remains the default and is not weakened --
it now checks the thing that can actually be wrong (plane count / real name
match / periodicity) instead of a comparison that could never pass.

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
import re
import sys
import pickle
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage.measure import regionprops_table
import tifffile

import spacec as sp

from spatia.analysis import channels as ch


# ─────────────────────────────────────────────────────────────────────────────
# PART 1 -- SEGMENTATION
# ─────────────────────────────────────────────────────────────────────────────

def _find_masked_tifs(processed_rois_dir: str) -> List[str]:
    """
    Recursively find masked ROI TIFFs, excluding preview images and OME-TIFF
    exports (e.g. "*_combined_markers.ome.tif").

    Added 2026-08-28: ".ome.tif" is a compound extension, so a naive
    f.endswith(".tif") check also matches OME-TIFF files -- these are
    QuPath-viewer exports (structured OME-XML metadata, not necessarily the
    same channel order/count as the plain multichannel TIFFs this pipeline
    is built around) and were never meant to be segmented. Without this
    exclusion, an "*.ome.tif" sitting in the same folder as its plain-TIFF
    counterpart would get segmented too -- wasted spacec/Mesmer compute at
    best, a channel-order mismatch halt at worst.
    """
    tif_files = []
    for root, _dirs, files in os.walk(processed_rois_dir):
        for f in files:
            # "._name.tif" files are macOS AppleDouble resource forks, created
            # whenever this data is copied through an HFS+/APFS volume. They
            # are ~4 KB of metadata, not images, and tifffile raises
            # "not a TIFF file b'\\x00\\x05\\x16\\x07'" on them. Added
            # 2026-09-03 after 134 of them entered a 276-file batch.
            if f.startswith("._"):
                continue
            if f.endswith(".tif") and "preview" not in f and ".ome.tif" not in f:
                tif_files.append(os.path.join(root, f))
    return tif_files


# ─────────────────────────────────────────────────────────────────────────────
# CHANNEL-ORDER VALIDATION
# run_cell_segmentation applies ONE channel_file to every masked ROI TIFF in
# the batch. That's only correct if every image really does share the same
# channel count/order.
#
# REVISED 2026-09-03. The original check compared the TIFF's embedded ImageJ
# 'Labels' against the panel file name-for-name. That works for ROIs written by
# roi_masking.py (which embeds real marker names), but it can NEVER pass for
# raw CODEX/PhenoCycler acquisitions, whose labels are scanner placeholders --
# "ch1".."ch4" repeating once per cycle -- carrying no marker identity at all.
# On the CRC TMA that rejected all 137 real cores.
#
# The check now dispatches on what the image actually is:
#   * real marker names embedded  -> compare name-for-name, as before
#   * generic placeholders        -> compare PLANE COUNT against the panel, and
#                                    verify nuclear periodicity on the pixels
#
# Plane count is the comparison with teeth: it catches a truncated acquisition,
# a panel from a different experiment, or a core genuinely missing cycles --
# the failure modes the name check was reaching for and missing.
# ─────────────────────────────────────────────────────────────────────────────

def _read_channel_file(channel_file_path: str) -> List[str]:
    """Parse channelnames.txt: one channel name per line, in stack order."""
    with open(channel_file_path) as f:
        return [line.strip() for line in f if line.strip()]


def _read_tiff_channel_labels(tiff_path: str) -> Optional[List[str]]:
    """
    Per-image channel order as written by roi_masking.py. Prefers the TIFF's
    own embedded ImageJ 'Labels' metadata (travels with the file even if it's
    moved/renamed); falls back to the sidecar "..._channel_info.txt" written
    next to it. Returns None if neither is present (e.g. masked ROIs produced
    before this check existed, or by some other tool) -- caller treats that
    as "can't verify," not as a mismatch.
    """
    try:
        with tifffile.TiffFile(tiff_path) as tif:
            ij_meta = tif.imagej_metadata or {}
            labels = ij_meta.get("Labels")
            if labels:
                # tifffile round-trips a single-element Labels list as a bare
                # string, not a length-1 list -- e.g. ["DAPI"] comes back as
                # "DAPI". list("DAPI") would silently explode that into
                # ['D','A','P','I'], so a real single-channel image (nuclear-
                # only panel, DAPI-only ROI) must be special-cased here.
                if isinstance(labels, str):
                    return [labels]
                return list(labels)
    except Exception:
        pass

    # Fallback: "{roi_experiment_group}_{slide_id}_x{x}_y{y}_channel_info.txt" --
    # roi_masking.py writes this WITHOUT the "_w{w}_h{h}_masked" suffix the
    # .tif filename itself has, so strip that suffix rather than assuming
    # a simple extension swap.
    basename = os.path.basename(tiff_path)
    m = re.match(r"^(.*_x\d+_y\d+)_w\d+_h\d+_masked\.tif$", basename)
    if not m:
        return None
    info_path = os.path.join(os.path.dirname(tiff_path), m.group(1) + "_channel_info.txt")
    if not os.path.exists(info_path):
        return None
    try:
        labels = {}
        with open(info_path) as f:
            for line in f:
                cm = re.match(r"Channel (\d+):\s*(.*)", line.strip())
                if cm:
                    labels[int(cm.group(1))] = cm.group(2)
        return [labels[i] for i in sorted(labels)] if labels else None
    except Exception:
        return None


def _check_channel_order(
    tiff_path: str,
    expected_channels: List[str],
    info: Optional["ch.StackInfo"] = None,
) -> Tuple[Optional[bool], str]:
    """
    Validate this image against the global channel_file.

    Returns (True, msg) on match, (False, msg) on a real mismatch (the file
    gets skipped), (None, msg) if nothing could be checked (proceeds
    unverified with a warning, doesn't block).

    Two paths, chosen by what the image carries:

    1. Embedded labels are real marker names (roi_masking.py output) --
       compare name-for-name, exactly as before.

    2. Embedded labels are scanner placeholders, or absent (raw CODEX) --
       marker identity was never in the file, so names cannot be compared.
       Compare the PLANE COUNT instead, using the file's true structure:
       a (23, 4, Y, X) hyperstack is 92 planes, not 23.
    """
    info = info or ch.inspect_stack(tiff_path)

    # Fail loud rather than guessing. An OME-TIFF or QPTIFF normally carries
    # its own channel names; if they can't be read, "unknown" must not be
    # treated as "unlabeled" -- applying the panel positionally to a file that
    # already had names would silently relabel every marker.
    if info.labels_unverifiable:
        return False, (
            "file format normally carries channel names (OME-TIFF/QPTIFF) but "
            "they could not be read — refusing to apply channel_file "
            "positionally, which would silently relabel markers. Export this "
            "image with readable channel metadata, or convert it to a plain "
            "multichannel .tif whose channel order you have verified."
        )

    if info.n_planes != len(expected_channels):
        return False, (
            f"channel COUNT mismatch — channel_file has {len(expected_channels)} "
            f"names, this image has {info.n_planes} planes ({info.describe()})"
        )

    if info.labels_are_generic():
        return None, (
            f"image carries only generic channel labels ({info.describe()}) — "
            f"marker identity comes from channel_file; plane count "
            f"({info.n_planes}) verified"
        )

    actual = _read_tiff_channel_labels(tiff_path)
    if actual is None:
        return None, "no per-image channel metadata found (older ROI export?) — cannot verify"

    mismatches = [
        f"position {i}: channel_file says '{e}', image says '{a}'"
        for i, (e, a) in enumerate(zip(expected_channels, actual)) if e != a
    ]
    if mismatches:
        shown = "; ".join(mismatches[:5])
        if len(mismatches) > 5:
            shown += f"; ... and {len(mismatches) - 5} more"
        return False, "channel ORDER/NAME mismatch — " + shown

    return True, "channels match"


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
    channel_check: bool = True,
    normalized_dir: Optional[str] = None,
    keep_normalized: bool = False,
    auto_panel_template: bool = True,
) -> Dict[str, object]:
    """
    Run spacec cell segmentation over every masked ROI TIFF found under
    processed_rois_dir. Idempotent: existing *_seg_output.pickle files are
    skipped. Ported from 03_segmentation.ipynb Part 1.

    Returns
    -------
    dict with keys:
        "outputs" : Dict[str, dict]   -- successful seg_output per image
        "errors"  : List[dict]        -- [{"file", "error"}, ...] for images
                                          that raised during sp.tl.cell_segmentation.
                                          Not raised automatically (unlike a
                                          channel mismatch, a one-off spacec/Mesmer
                                          failure on a single image shouldn't
                                          halt the whole batch) -- but always
                                          summarized in the printed output, and
                                          checked by validation.validate_segmentation
                                          after the step runs.

    channel_check : bool, default True
        channel_file_path is a SINGLE file applied to every image in the
        batch. If images don't all share the same channel count/order (e.g.
        different scan batches or panel revisions), applying one global
        channel_file silently mislabels markers for whichever images don't
        match -- spacec maps names onto the channel stack positionally, with
        no error. When True, each image's own channel order (embedded by
        roi_masking.py in the TIFF's 'Labels' metadata, or its sidecar
        "..._channel_info.txt") is checked against channel_file_path before
        segmenting. A mismatched image is SKIPPED (not segmented with wrong
        labels) and the whole run raises at the end so the mismatch can't
        pass silently. Images with no per-image metadata to check against
        (e.g. pre-existing ROI exports) proceed with a warning, unverified --
        same as the behavior before this check existed. Set False to restore
        that old, unchecked behavior entirely.

    normalized_dir : str, optional
        Where to write flattened, correctly-labeled copies of any image whose
        on-disk layout spacec cannot read directly (raw CODEX "TCYX"
        hyperstacks). Defaults to "<output_dir>/_normalized".

    keep_normalized : bool, default False
        A normalized copy is the same size as its source (~508 MB per CRC core,
        ~72 GB across all four TMA folders), so by default each one is deleted
        as soon as its image has been segmented -- peak extra disk is one file,
        not a second copy of the dataset. Set True to keep them for inspection,
        or to let a re-run skip the conversion.

    auto_panel_template : bool, default True
        If channel_file_path does not exist, write a positional template
        ("cyc01_ch1", ...) sized from the first image instead of hard-failing,
        so the mechanics can be smoke-tested. Loudly flagged: those are
        placeholders, not marker identities.
    """
    os.makedirs(output_dir, exist_ok=True)
    normalized_dir = normalized_dir or os.path.join(output_dir, "_normalized")

    if not os.path.exists(channel_file_path) and not auto_panel_template:
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
    channel_mismatches: List[dict] = []
    segmentation_errors: List[dict] = []
    n_unverified = 0
    n_already_done = 0
    n_normalized = 0

    # Resolved once, from the first readable image, then re-validated per image
    # by plane count. expected_channels stays the single source of marker
    # identity; panel_source records where it came from for the run log.
    expected_channels: Optional[List[str]] = None
    resolved_nuclei = nuclei_channel
    resolved_membrane = list(membrane_channel_list)
    periodicity_checked = False

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
            n_already_done += 1
            continue

        # ── Structure detection ──────────────────────────────────────────────
        # Reads only the TIFF header, so this is cheap enough to do per image
        # and it is what tells us a (23, 4, Y, X) hyperstack is 92 planes.
        # A file that isn't a readable TIFF at all (macOS "._" resource forks,
        # truncated copies) fails here and is recorded, instead of crashing the
        # batch several steps later inside spacec.
        try:
            info = ch.inspect_stack(input_file)
        except Exception as e:
            print(f"  ❌ Unreadable TIFF, skipping {filename}: {e}")
            segmentation_errors.append({"file": input_file, "error": f"unreadable TIFF: {e}"})
            continue

        # ── Panel resolution (once) ──────────────────────────────────────────
        if expected_channels is None:
            expected_channels = ch.ensure_panel(
                channel_file_path, info, auto_template=auto_panel_template
            )
            resolved = ch.resolve_channels(
                expected_channels, info,
                nuclei_channel=nuclei_channel,
                membrane_channel_list=membrane_channel_list,
            )
            resolved_nuclei = resolved["nuclei_channel"]
            resolved_membrane = resolved["membrane_channel_list"]
            groups = ch.classify_panel(expected_channels, info)
            print(
                f"Panel               : {len(expected_channels)} channels "
                f"({len(groups['marker'])} markers, {len(groups['nuclear'])} nuclear, "
                f"{len(groups['blank'])} blank/empty)"
            )
            print(f"Layout              : {info.describe()}")
            print(f"Nuclei channel      : {resolved_nuclei!r}  (config: {nuclei_channel!r})")
            print(f"Membrane channel(s) : {resolved_membrane}  (config: {list(membrane_channel_list)})")
            if groups["blank"]:
                print(
                    f"  note: planes {groups['blank'][:6]}"
                    f"{' ...' if len(groups['blank']) > 6 else ''} are blank/empty "
                    f"slots — expect zero variance; drop them before any "
                    f"per-channel normalization downstream."
                )

        if channel_check:
            ok, msg = _check_channel_order(input_file, expected_channels, info)
            if ok is False:
                print(f"  ❌ CHANNEL MISMATCH for {filename}: {msg}")
                print(f"     Skipping this image — segmenting it against "
                      f"channel_file_path as-is would mislabel its markers. "
                      f"Fix channelnames.txt (or this image's channel order), "
                      f"then re-run.")
                channel_mismatches.append({"file": input_file, "reason": msg})
                continue
            elif ok is None:
                n_unverified += 1
                print(f"  ⚠️  {filename}: {msg}")

        # ── Normalization ────────────────────────────────────────────────────
        # spacec reads the file itself and treats axis 0 as the channel axis, so
        # a (23, 4, Y, X) hyperstack makes it see 23 channels and index out of
        # bounds. Hand it a flat (92, Y, X) TIFF whose embedded labels are the
        # real marker names instead. Nothing to do for already-flat images.
        seg_input = input_file
        normalized_path = None
        if info.needs_flattening or info.labels_are_generic():
            normalized_path = os.path.join(normalized_dir, slide_id, filename)
            try:
                ch.normalize_stack_to_file(
                    input_file, normalized_path, expected_channels, info
                )
                seg_input = normalized_path
                n_normalized += 1
            except Exception as e:
                print(f"  ❌ Could not normalize {filename}: {e}")
                segmentation_errors.append({"file": input_file, "error": f"normalize failed: {e}"})
                continue

            # One pixel-level check per run that the cycle-major flatten was
            # right: in CODEX, channel 1 of every cycle re-images the same
            # nuclei, so those planes must correlate. If they don't, the panel
            # is being mapped onto a stack it doesn't describe.
            if channel_check and not periodicity_checked:
                periodicity_checked = True
                try:
                    qc = ch.verify_nuclear_periodicity(
                        ch.load_flat_stack(input_file, info), info
                    )
                    if qc.get("checked"):
                        status = "✓" if qc["passed"] else "❌"
                        print(f"  {status} nuclear periodicity: mean r={qc['mean_corr']:.3f} — {qc['reason']}")
                        if not qc["passed"]:
                            raise RuntimeError(
                                f"Nuclear-periodicity check failed on {filename}: "
                                f"{qc['reason']}. Refusing to segment the batch — "
                                f"every marker would carry the wrong name."
                            )
                except RuntimeError:
                    raise
                except Exception as e:
                    print(f"  ⚠️  periodicity check could not run: {e}")

        print(f"Segmenting: {seg_input}")
        try:
            seg_output = sp.tl.cell_segmentation(
                file_name=seg_input,
                channel_file=channel_file_path,
                output_dir=slide_output_dir,
                seg_method=seg_method,
                nuclei_channel=resolved_nuclei,
                output_fname=output_fname,
                membrane_channel_list=resolved_membrane,
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
            segmentation_errors.append({"file": input_file, "error": str(e)})
        finally:
            # Peak extra disk stays at one normalized file, not a second copy
            # of the dataset. Kept only if explicitly asked for.
            if normalized_path and not keep_normalized and os.path.exists(normalized_path):
                try:
                    os.remove(normalized_path)
                except OSError:
                    pass

    print("All segmentation processes completed!")

    # Always-visible summary -- this used to require grepping the console log
    # to find out whether anything failed; now it's one line, every run.
    print(
        f"\nSegmentation summary: {len(segmentation_outputs)} succeeded, "
        f"{len(segmentation_errors)} failed, {n_already_done} already done "
        f"(skipped), {len(channel_mismatches)} skipped for channel mismatch "
        f"— {len(tif_files)} total files found."
    )
    if n_normalized:
        kept = "kept" if keep_normalized else "deleted after segmenting"
        print(f"   {n_normalized} image(s) flattened to a labeled CYX stack first ({kept}).")
    if segmentation_errors:
        print(f"⚠️  {len(segmentation_errors)} image(s) FAILED segmentation "
              f"(see errors above) and produced no output:")
        for err in segmentation_errors:
            print(f"   {os.path.basename(err['file'])}: {err['error']}")
        print("   These will show up as missing *_mesmer_result.csv files — "
              "run_segmentation() does not raise on this automatically, but "
              "validation.validate_segmentation() checks for it after the step runs.")

    if n_unverified:
        print(f"⚠️  {n_unverified} image(s) segmented without per-image channel "
              f"verification (no embedded/sidecar channel metadata found).")

    if channel_mismatches:
        print("\n" + "=" * 72)
        print(f"❌ {len(channel_mismatches)} image(s) SKIPPED due to channel mismatch:")
        for m in channel_mismatches:
            print(f"   {os.path.basename(m['file'])}: {m['reason']}")
        print("=" * 72)
        raise RuntimeError(
            f"{len(channel_mismatches)} masked ROI TIFF(s) do not match "
            f"channel_file_path ({channel_file_path}) and were skipped rather "
            f"than segmented with wrong channel labels. See the list above. "
            f"Either fix channelnames.txt to match, or split these images into "
            f"a separate segmentation run with the correct channel_file. "
            f"A count mismatch usually means a truncated acquisition or a panel "
            f"from a different experiment — check the reported plane count "
            f"against len(channelnames.txt) before assuming the panel is wrong. "
            f"Pass channel_check=False to run_segmentation()/run_cell_segmentation() "
            f"to bypass this check (not recommended: it also disables the "
            f"nuclear-periodicity verification, which is the only pixel-level "
            f"evidence that the panel maps onto this stack correctly)."
        )

    return {
        "outputs": segmentation_outputs,
        "errors": segmentation_errors,
        # Resolved names, not the raw config values -- overlays and any
        # downstream step need the real panel entry, not the literal "auto".
        "nuclei_channel": resolved_nuclei,
        "membrane_channel_list": resolved_membrane,
        "channel_names": expected_channels or [],
    }


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
    segmentation.channel_check             default True -- verify each image's
                                            own channel order (from roi_masking.py's
                                            metadata) against channel_file before
                                            segmenting it; see run_cell_segmentation's
                                            docstring for what this catches.

    Returns
    -------
    dict with keys: segmentation_outputs (dict), segmentation_errors (list of
    {"file", "error"} for any image that raised during spacec segmentation --
    not raised automatically, but always summarized in the printed output and
    checked by validation.validate_segmentation() after this step runs),
    csv_files (list), overlay_files (list), output_dir (path)
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
    # "auto" resolves structurally: channel 1 of cycle 1, the CODEX nuclear
    # slot. Short names ("CD45") also resolve against the full panel entry
    # ("CD45 - hematopoietic cells"), so the config doesn't have to mirror
    # the panel file's exact punctuation.
    nuclei_channel = seg_cfg.get("nuclei_channel", "auto")
    membrane_channel_list = seg_cfg.get("membrane_channel_list", ["CD45"])
    compartment = seg_cfg.get("compartment", "whole-cell")
    input_format = seg_cfg.get("input_format", "Multichannel")
    resize_factor = seg_cfg.get("resize_factor", 1)
    size_cutoff = seg_cfg.get("size_cutoff", 0)
    do_overlays = seg_cfg.get("generate_overlays", True)
    channel_check = seg_cfg.get("channel_check", True)
    normalized_dir = seg_cfg.get("normalized_dir") or paths.get("converted_tif_dir")
    keep_normalized = seg_cfg.get("keep_normalized", False)
    auto_panel_template = seg_cfg.get("auto_panel_template", True)

    print("=" * 72)
    print("SPATIA CELL SEGMENTATION")
    print("=" * 72)
    print(f"Input (masked ROIs) : {processed_rois_dir}")
    print(f"Channel file        : {channel_file_path}")
    print(f"Output dir          : {output_dir}")
    print(f"Method              : {seg_method}")

    if not os.path.isdir(processed_rois_dir):
        raise FileNotFoundError(f"masked_roi_dir not found: {processed_rois_dir}")

    seg_result = run_cell_segmentation(
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
        channel_check=channel_check,
        normalized_dir=normalized_dir,
        keep_normalized=keep_normalized,
        auto_panel_template=auto_panel_template,
    )
    segmentation_outputs = seg_result["outputs"]
    segmentation_errors = seg_result["errors"]
    resolved_nuclei = seg_result.get("nuclei_channel", nuclei_channel)
    resolved_membrane = seg_result.get("membrane_channel_list", membrane_channel_list)

    csv_files = export_segmentation_to_csv(output_dir)

    overlay_files: List[str] = []
    if do_overlays:
        overlay_files = generate_overlays(output_dir, nucleus_channel=resolved_nuclei,
                                           additional_channels=resolved_membrane)

    return {
        "segmentation_outputs": segmentation_outputs,
        "segmentation_errors": segmentation_errors,
        "csv_files": csv_files,
        "overlay_files": overlay_files,
        "output_dir": output_dir,
        "nuclei_channel": resolved_nuclei,
        "membrane_channel_list": resolved_membrane,
        "channel_names": seg_result.get("channel_names", []),
    }
