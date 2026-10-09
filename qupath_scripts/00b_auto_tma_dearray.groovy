/*
===============================================================================
Automatic TMA Core Extraction (Dearray)
===============================================================================

OVERVIEW
--------
NEW script (added alongside, not replacing, the existing manual ROI workflow
in 00_ROI_extract_mask_project.groovy). That script exports mask PNGs from
annotations that ALREADY exist in a QuPath project -- it does not detect
tissue/cores itself. For whole-slide images (WSI), core/ROI detection is
assumed to already exist elsewhere in the WSI workflow. For Tissue
MicroArrays (TMA), this project had no automatic equivalent -- cores were
being defined by hand. This script closes that gap using QuPath's own
built-in TMA Dearrayer, run headlessly (no dialogs) so it fits the same
Cirro-reproducibility bar as the rest of the pipeline.

WHAT THIS DOES NOT DO
----------------------
It does not reimplement core detection from scratch. QuPath's dearrayer
already does grid-fitting + per-core intensity thresholding reliably; this
script just drives it with parameters instead of the GUI "TMA Dearrayer"
dialog, and then hands off to the existing mask-export logic so the rest
of the pipeline (roi_masking.py) sees the same directory contract either
way -- manual WSI ROI or automatic TMA core.

-------------------------------------------------------------------------------
REQUIREMENTS
-------------------------------------------------------------------------------
- QuPath (tested against the TMA Dearrayer plugin shipped with QuPath's
  ImageJ extension: qupath.imagej.detect.dearray.TMADearrayerPluginIJ).
- Run via QuPath's command-line interface for headless/Cirro use. The
  executable is version-specific (e.g. "QuPath-0.7.0", or the console/.app
  variant on Windows/Mac -- see QuPath's own command-line docs), and the
  script path is the LAST argument after the "script" subcommand and its
  flags, not the first:
      QuPath-<version> script --image <path/to/tma_slide> 00b_auto_tma_dearray.groovy
  or from Automate -> Script Editor inside an open project (no dialogs
  are shown either way -- all parameters are set below, not prompted for).
- Image must be a single scanned TMA slide containing the full grid of
  cores (this is the whole-TMA-scan stage, upstream of the individual
  per-core files already used elsewhere in this pipeline, e.g. the
  TMA_A / TMA_B / DII_TMA_A / DII_TMA_B per-core exports).

-------------------------------------------------------------------------------
PARAMETERS (edit these, no interactive prompts)
-------------------------------------------------------------------------------
CORE_DIAMETER_MM   Physical core diameter in mm. Get this from the TMA's
                   construction spec (ask whoever built the array if unknown --
                   do NOT guess; a wrong diameter silently mis-detects cores).
GRID_HORIZONTAL    Column label range, e.g. "1-16" for 16 columns.
GRID_VERTICAL      Row label range, e.g. "A-J" for 10 rows.
LABEL_ORDER        "Row first" or "Column first" -- must match how the array
                   was physically constructed/labeled.
DENSITY_THRESHOLD  0-1. Lower = more permissive core detection (keeps faint/
                   sparse cores); higher = stricter (may drop damaged cores).
                   Default below (5) is the QuPath-recommended starting point,
                   NOT verified against this project's specific TMAs -- treat
                   as a starting point to tune, not a settled value (see
                   SPATIA_PIPELINE_LOG.md decision queue).

-------------------------------------------------------------------------------
OUTPUT
-------------------------------------------------------------------------------
- One TMACoreObject per detected core, added to the current image's
  hierarchy, each auto-labeled by grid position (e.g. "A-1", "A-2", ...
  matching LABEL_ORDER).
- Missing/failed cores are still added (marked "missing") so downstream
  code sees a complete, consistently-shaped grid rather than silent gaps.
- After dearray, this script calls the SAME mask-export routine used by
  00_ROI_extract_mask_project.groovy (color mask, binary mask, roi_labels.txt,
  individual per-core masks) so roi_masking.py's directory contract
  (binary_mask.png / roi_labels.txt / individual_masks/) is identical
  regardless of which script produced it.

-------------------------------------------------------------------------------
STATUS
-------------------------------------------------------------------------------
Written but NOT executed against a real slide -- this sandbox has no QuPath
installation to run it in (see SPATIA_PIPELINE_LOG.md for the environment
note). Needs a real test run against one of the raw TMA slide scans
(wherever the pre-dearray, full-grid scan lives -- the per-core files
already in this project's OneDrive folder appear to be POST-dearray, so
the true input for this script may need to be located or re-exported;
flagged as an open question, not assumed).

-------------------------------------------------------------------------------
EXAMPLE USAGE
-------------------------------------------------------------------------------
1. Command line (headless, Cirro-style) -- script path comes LAST:
     QuPath-<version> script --image /data/TMA_A_full_scan.tif 00b_auto_tma_dearray.groovy
2. Or inside QuPath: Automate -> Script Editor -> paste -> Run
===============================================================================
*/

import qupath.lib.images.servers.ImageServer
import qupath.lib.objects.hierarchy.TMAGrid
import groovy.json.JsonOutput

// ── Parameters (edit per TMA; no dialogs) ───────────────────────────────────
def CORE_DIAMETER_MM  = 1.0        // TODO: confirm against this TMA's construction spec
def GRID_HORIZONTAL   = "1-16"     // column labels
def GRID_VERTICAL     = "A-J"      // row labels
def LABEL_ORDER       = "Row first"
def DENSITY_THRESHOLD = 5          // QuPath default starting point -- not yet tuned

// Output directory for the mask-export handoff (same layout 00_ROI_extract
// produces, so roi_masking.py doesn't care which script ran).
def OUTPUT_DIR = buildFilePath(PROJECT_BASE_DIR, "auto_tma_extraction")

// ── Run the dearrayer headlessly ────────────────────────────────────────────
def imageData = getCurrentImageData()
if (imageData == null) {
    print "ERROR: No current image. Run this via QuPath's CLI with --image, " +
          "or open the image in a project first."
    return
}

def dearrayArgs = JsonOutput.toJson([
    coreDiameterMM : CORE_DIAMETER_MM,
    labelsHorizontal: GRID_HORIZONTAL,
    labelsVertical  : GRID_VERTICAL,
    labelOrder      : LABEL_ORDER,
    densityThreshold: DENSITY_THRESHOLD,
    fixNonStraightGrids: false,
])

print "Running TMA Dearrayer with: ${dearrayArgs}"
runPlugin('qupath.imagej.detect.dearray.TMADearrayerPluginIJ', imageData, dearrayArgs)

def hierarchy = imageData.getHierarchy()
def tmaGrid = hierarchy.getTMAGrid()
if (tmaGrid == null) {
    print "ERROR: Dearray produced no TMA grid. Check CORE_DIAMETER_MM and " +
          "DENSITY_THRESHOLD against this slide before re-running."
    return
}

def cores = tmaGrid.getTMACoreList()
def nMissing = cores.count { it.isMissing() }
print "Dearray complete: ${cores.size()} core slots (${nMissing} flagged missing)."

// Persist the detected grid back into the project data so the mask-export
// step (or a human spot-check in QuPath) sees exactly what was detected.
fireHierarchyUpdate()
print "TMA grid written to image hierarchy. Save the project (or let the " +
      "CLI runner auto-save) before running the mask-export step."

// ── Hand off to the existing mask-export logic ──────────────────────────────
// Deliberately NOT duplicated here: 00_ROI_extract_mask_project.groovy already
// implements color-mask / binary-mask / individual-mask export against
// whatever annotations exist in the hierarchy. TMACoreObjects behave as
// annotated regions for that purpose, so once the grid above is written,
// run 00_ROI_extract_mask_project.groovy (with "All annotations" selected)
// against this same image to get the binary_mask.png / roi_labels.txt /
// individual_masks/ that roi_masking.py expects. Kept as two scripts rather
// than merged into one so the existing, working export logic isn't touched
// or duplicated -- only the detection step was missing, and that's the only
// thing this file adds.
print "Next: run 00_ROI_extract_mask_project.groovy against this image to " +
      "export masks from the ${cores.size()} detected cores."
