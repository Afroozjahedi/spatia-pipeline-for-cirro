"""
spatia/analysis/channels.py
===========================
Channel-structure detection, panel resolution, and stack normalization.

NEW FILE -- 2026-09-04 (Afrouz + Claude consulting session)
--------------------------------------------------------------
Created to fix the 2026-09-03 real-data segmentation failure (all 137 real CRC
TMA cores rejected as "channel mismatch"; see segmentation.py's changelog for
the full incident). Nothing in this file replaces prior logic -- it is new
capability that segmentation.py now calls into.

Two capabilities were added in the same sitting and are both below:

  (a) CODEX/PhenoCycler hyperstack support -- inspect_stack(), load_flat_stack(),
      panel_template(), ensure_panel(), classify_panel(), resolve_channel(s),
      normalize_stack_to_file(), verify_nuclear_periodicity(). This is the fix
      for the 2026-09-03 failure itself. Full rationale in the next section.

  (b) OME-TIFF / QPTIFF name-reading + a safety net (_ome_channel_names(),
      _qptiff_channel_names(), StackInfo.labels_unverifiable) -- added the same
      day, BEFORE any real .ome.tif or .qptiff file was run through this
      pipeline, in response to a direct question ("does this work for OME-TIFF
      or QPTIFF?"). Verified against synthetic OME-TIFFs (named, unnamed, and
      3-level pyramidal) in this session; NOT verified against a real .qptiff --
      none was available. The distinction that matters: an OME-TIFF or QPTIFF
      whose channel names could not be parsed is marked `labels_unverifiable`,
      which is NOT the same as `labels_are_generic()` (raw CODEX legitimately
      has no marker identity in the file). Collapsing those two cases would
      mean applying the panel positionally to a file that may already have had
      real, different, channel names -- silently mislabeling every marker.
      `labels_are_generic()` is explicitly False whenever `labels_unverifiable`
      is True, and segmentation.py's channel check fails loud on that case
      rather than falling through to the count-only comparison.

WHY THIS MODULE EXISTS
-----------------------------------------
The CRC TMA cores are raw CODEX/PhenoCycler output. Each per-core .tif is an
ImageJ *hyperstack* with axes "TCYX" -- shape (23, 4, 1440, 1920) -- where the
"T" (frames) axis is really the CYCLE axis and "C" is the 4 fluorescence
channels imaged per cycle. 23 cycles x 4 channels = the 92-plane panel.

Two consequences, both of which broke the 2026-09-03 run:

1. The embedded ImageJ 'Labels' metadata reads ["ch1","ch2","ch3","ch4", ...]
   repeating 23 times. It carries NO biological marker identity -- the scanner
   never knew it. So a check that compares those labels against a marker panel
   can never pass, for any correctly-acquired file. That check was rejecting
   all 137 real cores.

2. spacec's segmentation reads the file itself and treats axis 0 as the channel
   axis. On a (23, 4, Y, X) array that means it sees 23 "channels" and raises
   "index 23 is out of bounds for axis 0 with size 23" the moment anything
   walks past plane 22.

The fix is to normalize each stack to a flat, honestly-labeled (C, Y, X) TIFF
before handing it to spacec, and to derive as much of the configuration as
possible from the file's own structure instead of asking the user to hand-type
channel names.

WHAT CAN AND CANNOT BE AUTOMATED
---------------------------------
Automatable, because the file knows it:
  - the axis layout and plane count
  - the cycle count and channels-per-cycle
  - which planes are the nuclear stain (channel 1 of every cycle, in CODEX)
  - which plane is the last marker (the last one)
  - a positional template panel ("cyc01_ch1", ...) when no panel is supplied

NOT automatable, because the file genuinely does not contain it:
  - the biological identity of each marker

Marker identity is experiment metadata that lives in the acquisition sheet, not
in the pixels. A pipeline that guessed it would be inventing data. So the panel
file stays an explicit, one-time, per-experiment input -- but it is the ONLY
manual artifact, and this module validates its length against the real stack so
a wrong panel fails loudly at the first image instead of silently mislabeling
every cell in the run.

Panel-order convention
----------------------
Flattening is cycle-major (T outer, C inner), matching how the panel is written
in acquisition order:
    plane 0  = cycle 1 / channel 1   (nuclear)
    plane 1  = cycle 1 / channel 2
    ...
    plane 4  = cycle 2 / channel 1   (nuclear)
Verified on TMA_A/reg011_X01_Y01_Z09.tif: all 23 nuclear positions correlate
r = 0.90-0.92 with plane 0, DRAQ5 (plane 91, also a nuclear dye) at r = 0.87,
and the three Blank planes are identically zero. See verify_nuclear_periodicity().
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import tifffile


# Channel 1 of each CODEX cycle is the nuclear stain. Used for "auto" nuclei
# resolution and for the periodicity QC below.
NUCLEAR_CHANNEL_WITHIN_CYCLE = 0

# Panel entries matching these are non-informative acquisition slots, not
# markers. Kept in the panel (positions must stay aligned) but reported so
# downstream steps can drop them before any per-channel normalization -- they
# are identically zero and will produce divide-by-zero NaNs.
_BLANK_RE = re.compile(r"^(blank|empty)", re.IGNORECASE)

# Nuclear-stain names, including the "HOCHST13" typo present in the CRC panel.
_NUCLEAR_RE = re.compile(r"^(hoechst|hochst|dapi|draq5)", re.IGNORECASE)


def _ome_channel_names(tif: "tifffile.TiffFile") -> Optional[List[str]]:
    """
    Channel names from OME-XML. OME stores them as <Channel Name="..."> inside
    <Pixels>, NOT in ImageJ's 'Labels' tag -- so a file can be fully labeled
    and still look nameless to an ImageJ-only reader.

    Parsed with the stdlib XML parser rather than a regex over the whole
    document: <Image Name="..."> and <Plate Name="..."> also carry Name
    attributes, and a regex happily collects those too, silently shifting
    every channel by one.
    """
    xml = tif.ome_metadata
    if not xml:
        return None
    try:
        import xml.etree.ElementTree as ET

        root = ET.fromstring(xml)
        names = [
            el.get("Name")
            for el in root.iter()
            if el.tag.rsplit("}", 1)[-1] == "Channel"
        ]
    except Exception:
        return None
    names = [n for n in names if n]
    return names or None


def _qptiff_channel_names(tif: "tifffile.TiffFile") -> Optional[List[str]]:
    """
    Channel names from an Akoya/PhenoImager .qptiff. Each full-resolution page
    carries its own <PerkinElmer-QPI-ImageDescription> XML with a <Name> (and
    often <Biomarker>) element. NOT VERIFIED against a real .qptiff -- no such
    file was available when this was written -- so a parse failure here is
    reported as "unverifiable", never silently treated as "unlabeled".
    """
    try:
        import xml.etree.ElementTree as ET

        names: List[str] = []
        for page in tif.series[0].pages:
            desc = getattr(page, "description", None)
            if not desc or "<" not in desc:
                return None
            root = ET.fromstring(desc)
            found = None
            for tag in ("Biomarker", "Name"):
                el = root.find(f".//{tag}")
                if el is not None and el.text:
                    found = el.text.strip()
                    break
            if not found:
                return None
            names.append(found)
    except Exception:
        return None
    return names or None


@dataclass
class StackInfo:
    """What a multichannel TIFF's own structure tells us about itself."""

    path: str
    axes: str
    shape: Tuple[int, ...]
    n_planes: int
    n_cycles: int
    n_per_cycle: int
    height: int
    width: int
    embedded_labels: Optional[List[str]]
    needs_flattening: bool
    # Where embedded_labels came from: "imagej", "ome", "qptiff", or None.
    label_source: Optional[str] = None
    # True when the file is in a format that normally DOES carry channel names
    # (OME-TIFF, QPTIFF) but they could not be read. Distinct from "this file
    # has no names": one means unknown, the other means none exist. Applying a
    # panel positionally is safe in the second case and dangerous in the first.
    labels_unverifiable: bool = False
    n_pyramid_levels: int = 1

    def describe(self) -> str:
        layout = f"{self.axes} {self.shape}"
        if self.needs_flattening:
            layout += (
                f" -> flattened to ({self.n_planes}, {self.height}, {self.width}) "
                f"[{self.n_cycles} cycles x {self.n_per_cycle} channels]"
            )
        if self.n_pyramid_levels > 1:
            layout += f", {self.n_pyramid_levels} pyramid levels (using level 0)"
        if self.label_source:
            layout += f", names from {self.label_source}"
        return layout

    def labels_are_generic(self) -> bool:
        """
        True when this file demonstrably carries no marker identity -- either
        placeholder labels ("ch1".."ch4", "Channel 1", "C1") or none at all in
        a format that never stores them. That is the normal state for raw
        CODEX output and is NOT an error: it is precisely why an external panel
        file is required.

        Explicitly False when labels could not be verified (see
        labels_unverifiable). "I don't know what this file says" must never be
        collapsed into "this file says nothing" -- that is the difference
        between applying the panel to an unlabeled stack, which is correct, and
        overwriting real marker names positionally, which silently mislabels
        every cell.
        """
        if self.labels_unverifiable:
            return False
        if not self.embedded_labels:
            return True
        return all(
            re.fullmatch(r"(ch|channel|c)\s*\d+", lab.strip(), re.IGNORECASE)
            for lab in self.embedded_labels
        )


# ─────────────────────────────────────────────────────────────────────────────
# STRUCTURE DETECTION
# ─────────────────────────────────────────────────────────────────────────────

def inspect_stack(path: str) -> StackInfo:
    """
    Read a TIFF's structure without loading its pixels.

    Handles the layouts this pipeline actually encounters:
      "CYX"  -- already flat (nothing to do)
      "TCYX" -- CODEX hyperstack, T = cycles          (the CRC TMA case)
      "ZCYX" -- z-stack per channel, same flattening
      "YX"   -- single plane

    Any other axis order raises rather than guessing. That strictness is
    deliberate: a "CTYX" file (channel outer, cycle inner) has the SAME shape
    family as "TCYX" but the opposite memory order, and reshaping it as though
    it were cycle-major would silently permute every marker in the panel. A
    loud failure on an unrecognized layout is far cheaper than a run that
    completes with scrambled channel identities.
    """
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        axes = series.axes
        shape = tuple(int(s) for s in series.shape)
        ij_meta = tif.imagej_metadata or {}
        n_levels = len(getattr(series, "levels", [series]))
        is_ome = bool(getattr(tif, "is_ome", False))
        is_qptiff = bool(getattr(tif, "is_qpi", False)) or path.lower().endswith(".qptiff")

        # Channel names live in a different place in every format. Try each in
        # turn and record which one answered, so a later failure can say
        # "OME file whose names I couldn't read" rather than "no names".
        labels = ij_meta.get("Labels")
        if isinstance(labels, str):      # tifffile returns a bare str for len-1
            labels = [labels]
        elif labels is not None:
            labels = list(labels)
        label_source = "imagej" if labels else None

        if not labels and is_ome:
            labels = _ome_channel_names(tif)
            label_source = "ome" if labels else None
        if not labels and is_qptiff:
            labels = _qptiff_channel_names(tif)
            label_source = "qptiff" if labels else None

    # A format that normally carries names, from which none could be read.
    labels_unverifiable = bool((is_ome or is_qptiff) and not labels)

    if axes == "YX":
        h, w = shape
        n_planes, n_cycles, n_per_cycle, needs_flat = 1, 1, 1, False
    elif axes == "CYX":
        c, h, w = shape
        n_planes, n_cycles, n_per_cycle, needs_flat = c, 1, c, False
    elif axes in ("TCYX", "ZCYX"):
        outer, c, h, w = shape
        n_planes, n_cycles, n_per_cycle, needs_flat = outer * c, outer, c, True
    else:
        raise ValueError(
            f"Unhandled TIFF axis order {axes!r} (shape {shape}) in {path}. "
            "spatia.analysis.channels supports YX, CYX, TCYX and ZCYX. Add an "
            "explicit branch for this layout rather than assuming a reshape "
            "order -- guessing would silently permute the marker panel."
        )

    return StackInfo(
        path=path,
        axes=axes,
        shape=shape,
        n_planes=n_planes,
        n_cycles=n_cycles,
        n_per_cycle=n_per_cycle,
        height=int(h),
        width=int(w),
        embedded_labels=labels,
        needs_flattening=needs_flat,
        label_source=label_source,
        labels_unverifiable=labels_unverifiable,
        n_pyramid_levels=n_levels,
    )


def load_flat_stack(path: str, info: Optional[StackInfo] = None) -> np.ndarray:
    """Load a TIFF as a flat (C, Y, X) array in acquisition (cycle-major) order."""
    info = info or inspect_stack(path)
    with tifffile.TiffFile(path) as tif:
        arr = tif.series[0].asarray()

    if info.axes == "YX":
        arr = arr[np.newaxis, ...]
    elif info.needs_flattening:
        arr = arr.reshape(info.n_planes, info.height, info.width)
    return arr


# ─────────────────────────────────────────────────────────────────────────────
# PANEL FILE
# ─────────────────────────────────────────────────────────────────────────────

def read_panel(panel_path: str) -> List[str]:
    """One channel name per line, in stack order. Blank lines ignored."""
    with open(panel_path) as f:
        return [line.strip() for line in f if line.strip()]


def panel_template(info: StackInfo) -> List[str]:
    """
    Positional placeholder names derived purely from the file's structure --
    "cyc01_ch1", "cyc01_ch2", ... These are honest: they say where a plane sits
    without pretending to know what it stains. Written when no panel file
    exists so a run can be smoke-tested, and so the user has a correctly-sized
    file to fill in rather than counting planes by hand.
    """
    if info.n_cycles <= 1:
        return [f"ch{i + 1}" for i in range(info.n_planes)]
    return [
        f"cyc{cycle + 1:02d}_ch{ch + 1}"
        for cycle in range(info.n_cycles)
        for ch in range(info.n_per_cycle)
    ]


def ensure_panel(panel_path: str, info: StackInfo, auto_template: bool = True) -> List[str]:
    """
    Return the marker panel for this stack, creating a template if none exists.

    The length check is the validation that actually has teeth: it catches a
    truncated acquisition, a panel file from a different experiment, or a core
    that really is missing cycles -- the failure modes the old name-equality
    check was trying and failing to catch.
    """
    if os.path.exists(panel_path):
        names = read_panel(panel_path)
        if len(names) != info.n_planes:
            raise RuntimeError(
                f"Panel/stack size mismatch: {panel_path} lists {len(names)} "
                f"channel names but {os.path.basename(info.path)} has "
                f"{info.n_planes} planes ({info.describe()}). Either this image "
                f"is from a different panel/acquisition, or the panel file is "
                f"wrong. Not segmenting -- a positional mapping between "
                f"mismatched lengths would mislabel markers."
            )
        return names

    if not auto_template:
        raise FileNotFoundError(
            f"Channel panel file not found: {panel_path}. Segmentation cannot "
            "proceed without a channel-name mapping."
        )

    names = panel_template(info)
    os.makedirs(os.path.dirname(os.path.abspath(panel_path)) or ".", exist_ok=True)
    with open(panel_path, "w") as f:
        f.write("\n".join(names) + "\n")
    print(
        f"  ⚠️  No panel file at {panel_path} — wrote a {len(names)}-line "
        f"positional template derived from {os.path.basename(info.path)} "
        f"({info.describe()}).\n"
        f"      Marker identity cannot be read from the image; replace these "
        f"placeholders with the real acquisition panel before trusting any "
        f"marker-level result."
    )
    return names


def classify_panel(names: Sequence[str], info: StackInfo) -> Dict[str, List[int]]:
    """
    Split panel positions into nuclear / blank / marker groups.

    'blank' matters downstream: those planes are identically zero, so any
    per-channel standard-deviation normalization divides by zero and yields NaN
    columns. Better to know their indices up front than to debug NaNs later.
    """
    nuclear = [i for i, n in enumerate(names) if _NUCLEAR_RE.match(n)]
    blank = [i for i, n in enumerate(names) if _BLANK_RE.match(n)]
    structural_nuclear = list(
        range(NUCLEAR_CHANNEL_WITHIN_CYCLE, info.n_planes, info.n_per_cycle)
    ) if info.n_per_cycle > 1 else []
    marker = [i for i in range(len(names)) if i not in set(nuclear) | set(blank)]
    return {
        "nuclear": nuclear,
        "structural_nuclear": structural_nuclear,
        "blank": blank,
        "marker": marker,
    }


# ─────────────────────────────────────────────────────────────────────────────
# CONFIG RESOLUTION -- lets the YAML say "auto" or a short marker name
# ─────────────────────────────────────────────────────────────────────────────

def resolve_channel(spec, names: Sequence[str], info: StackInfo, role: str = "channel") -> str:
    """
    Turn a config value into an exact panel name.

    Accepted forms:
      "auto"          nuclei  -> channel 1 of cycle 1 (CODEX convention)
                      last    -> the final plane
      "CD45"          unique case-insensitive prefix/substring match, so the
                      config can stay short while the panel keeps its full
                      "CD45 - hematopoietic cells" text
      "HOECHST1 (C1)" exact match
      12              zero-based plane index

    Ambiguity raises. "CD4" legitimately prefix-matches CD44, CD45, CD45RA and
    CD45RO in this panel; picking the first would be a coin flip on a result
    nobody would notice was wrong.
    """
    if isinstance(spec, int):
        if not 0 <= spec < len(names):
            raise IndexError(f"{role} index {spec} outside panel of {len(names)} channels")
        return names[spec]

    text = str(spec).strip()

    if text.lower() == "auto":
        if role == "last_marker":
            return names[-1]
        if role == "nuclei_channel":
            idx = NUCLEAR_CHANNEL_WITHIN_CYCLE if info.n_per_cycle > 1 else 0
            return names[idx]
        raise ValueError(f"'auto' is not defined for role {role!r}")

    if text in names:
        return text

    lowered = [n.lower() for n in names]
    target = text.lower()

    exact_ci = [n for n, low in zip(names, lowered) if low == target]
    if len(exact_ci) == 1:
        return exact_ci[0]

    prefix = [n for n, low in zip(names, lowered) if low.startswith(target)]
    if len(prefix) == 1:
        return prefix[0]

    contains = [n for n, low in zip(names, lowered) if target in low]
    if len(contains) == 1:
        return contains[0]

    candidates = prefix or contains
    if candidates:
        raise ValueError(
            f"{role} {spec!r} is ambiguous — matches {len(candidates)} panel "
            f"entries: {candidates}. Use the full channel name or a plane index."
        )
    raise ValueError(
        f"{role} {spec!r} not found in the {len(names)}-channel panel. "
        f"First few entries: {list(names[:5])}"
    )


def resolve_channels(
    names: Sequence[str],
    info: StackInfo,
    nuclei_channel="auto",
    membrane_channel_list: Optional[Sequence] = None,
    last_marker="auto",
) -> Dict[str, object]:
    """Resolve every channel-name config value in one place, with one error style."""
    resolved_nuclei = resolve_channel(nuclei_channel, names, info, role="nuclei_channel")
    resolved_membrane = [
        resolve_channel(m, names, info, role="membrane_channel_list")
        for m in (membrane_channel_list or [])
    ]
    resolved_last = resolve_channel(last_marker, names, info, role="last_marker")
    return {
        "nuclei_channel": resolved_nuclei,
        "membrane_channel_list": resolved_membrane,
        "last_marker": resolved_last,
    }


# ─────────────────────────────────────────────────────────────────────────────
# QC -- turns the manual 2026-09-03 correlation check into an automatic one
# ─────────────────────────────────────────────────────────────────────────────

def verify_nuclear_periodicity(
    arr: np.ndarray,
    info: StackInfo,
    min_corr: float = 0.5,
    max_planes: int = 6,
) -> Dict[str, object]:
    """
    Confirm the flattening order actually produced a cycle-major stack.

    In CODEX, channel 1 of every cycle re-images the nuclear stain of the same
    tissue. So after a correct flatten, planes 0, n, 2n, ... must all correlate
    strongly with plane 0. If the reshape order were wrong, that periodicity
    disappears -- which makes this a direct, pixel-level test of the one
    assumption that cannot be checked from metadata.

    Cheap by design (a handful of planes, subsampled) so it can run per-slide
    rather than being a thing someone remembers to do by hand.
    """
    if info.n_per_cycle <= 1 or info.n_cycles <= 1:
        return {"checked": False, "reason": "not a cycled acquisition"}

    step = info.n_per_cycle
    nuclear_idx = list(range(0, info.n_planes, step))[:max_planes]
    ref = arr[nuclear_idx[0]].astype(np.float64).ravel()[::4]

    correlations = {}
    for i in nuclear_idx[1:]:
        other = arr[i].astype(np.float64).ravel()[::4]
        if ref.std() == 0 or other.std() == 0:
            correlations[i] = float("nan")
            continue
        correlations[i] = float(np.corrcoef(ref, other)[0, 1])

    finite = [r for r in correlations.values() if np.isfinite(r)]
    mean_corr = float(np.mean(finite)) if finite else float("nan")
    passed = bool(finite) and mean_corr >= min_corr
    return {
        "checked": True,
        "passed": passed,
        "mean_corr": mean_corr,
        "per_plane": correlations,
        "reason": (
            "nuclear channels repeat as expected"
            if passed
            else f"nuclear periodicity weak (mean r={mean_corr:.3f} < {min_corr}) — "
                 "the flatten order or the panel may not match this acquisition"
        ),
    }


# ─────────────────────────────────────────────────────────────────────────────
# NORMALIZATION -- what spacec actually gets handed
# ─────────────────────────────────────────────────────────────────────────────

def normalize_stack_to_file(
    src_path: str,
    dst_path: str,
    names: Sequence[str],
    info: Optional[StackInfo] = None,
    overwrite: bool = False,
) -> str:
    """
    Write a flat (C, Y, X) ImageJ TIFF whose embedded 'Labels' are the real
    marker names, so downstream tools -- spacec included -- read the stack the
    way the panel describes it.

    Written as a plain .tif, deliberately NOT .ome.tif: segmentation's own
    _find_masked_tifs() excludes "*.ome.tif", so an OME export here would be
    silently skipped and the run would report zero files.

    ImageJ hyperstack format is uncompressed and contiguous, so the output is
    roughly the size of the input (~508 MB per CRC core). Callers that don't
    want a permanent second copy should write into a scratch directory and
    delete after segmenting -- see segmentation.run_cell_segmentation.
    """
    info = info or inspect_stack(src_path)
    if len(names) != info.n_planes:
        raise RuntimeError(
            f"Cannot normalize {src_path}: {info.n_planes} planes but "
            f"{len(names)} channel names."
        )

    if os.path.exists(dst_path) and not overwrite:
        return dst_path

    os.makedirs(os.path.dirname(os.path.abspath(dst_path)) or ".", exist_ok=True)
    arr = load_flat_stack(src_path, info)

    tmp_path = dst_path + ".partial"
    tifffile.imwrite(
        tmp_path,
        arr,
        imagej=True,
        metadata={"axes": "CYX", "Labels": list(names), "mode": "composite"},
    )
    os.replace(tmp_path, dst_path)   # atomic: a killed job never leaves a half file
    return dst_path
