/**
 * # Notebook 04-2: Preprocessing Visualization (QuPath)
 *
 * **File:** `notebook/04-2_preprocessing_visualization.groovy`
 *
 * Imports pre-processed cell segmentation data (TSV export) into QuPath,
 * applies color-coded classification overlays, and corrects local crop
 * coordinates to global whole-slide image (WSI) coordinates.
 *
 * ---
 *
 * ## Overview
 *
 *   Load TSV export (centroid_x/y, classification, area, DAPI, offsets)
 *          ↓
 *   Parse column indices from TSV header
 *          ↓
 *   Determine tile offsets
 *   (from TSV columns OR filename fallback _x{X}_y{Y}_w{W}_h{H})
 *          ↓
 *   Convert local crop coords → global WSI coords
 *   (globalX = centroid_x + offset_x)
 *          ↓
 *   Create Ellipse ROI per cell (radius = sqrt(area/π), min 3px)
 *          ↓
 *   Assign PathClass with color coding
 *          ↓
 *   Add detections to QuPath hierarchy & refresh
 *
 * ---
 *
 * ## Configuration
 *
 * Only ONE line needs to be changed per run:
 *
 *   def tsvPath = "/path/to/your_qupath_cells.tsv"
 *
 * ---
 *
 * ## Input Requirements
 *
 * The TSV file must contain the following columns:
 *
 *   Column           | Type    | Description
 *   -----------------|---------|------------------------------------------
 *   centroid_x       | double  | Local X coordinate within the crop tile
 *   centroid_y       | double  | Local Y coordinate within the crop tile
 *   classification   | string  | Cell class label (see Classification Labels)
 *   area             | double  | Cell area in pixels²
 *   DAPI             | double  | DAPI channel intensity value
 *   offset_x         | double  | (optional) Global X tile offset in pixels
 *   offset_y         | double  | (optional) Global Y tile offset in pixels
 *
 * > **Note:** If `offset_x` / `offset_y` columns are absent, offsets are
 * > parsed automatically from the filename pattern:
 * >   `_x{X}_y{Y}_w{W}_h{H}`
 *
 * ---
 *
 * ## Classification Labels & Color Map
 *
 *   Label                   | Color        | Meaning
 *   ------------------------|--------------|------------------------------
 *   Included                | Green        | Valid cell kept for analysis
 *   Excl_SmallArea          | Red          | Excluded — area too small
 *   Excl_LowDAPI            | Orange       | Excluded — low DAPI signal
 *   Excl_SmallArea_LowDAPI  | Purple       | Excluded — both criteria
 *   Excl_Noise              | Blue         | Excluded — noise artifact
 *
 * ---
 *
 * ## Coordinate Correction
 *
 * Cells are segmented on individual crop tiles. To display them correctly
 * on the full WSI, local coordinates are converted:
 *
 *   globalX = centroid_x + offset_x
 *   globalY = centroid_y + offset_y
 *
 * Offset source priority:
 *   1. TSV column  `offset_x` / `offset_y`  (preferred)
 *   2. Filename pattern  `_x3547_y3256_w9046_h9152`  (fallback)
 *
 * ---
 *
 * ## ROI Construction
 *
 * Each cell is represented as an ellipse ROI:
 *
 *   radius = max( sqrt(area / π),  3.0 )   // minimum radius of 3px
 *   ROI = EllipseROI(globalX - r, globalY - r, 2r, 2r)
 *
 * ---
 *
 * ## Output
 *
 *   Detection objects added to QuPath hierarchy with measurements:
 *
 *   Measurement   | Description
 *   --------------|-------------------------------
 *   Area          | Cell area (pixels²)
 *   DAPI          | DAPI channel intensity
 *   centroid_x    | Global WSI X coordinate
 *   centroid_y    | Global WSI Y coordinate
 *   offset_x      | Tile X offset applied
 *   offset_y      | Tile Y offset applied
 *
 * ---
 *
 * ## Dependencies
 *
 *   Package                              | Purpose
 *   -------------------------------------|--------------------------------
 *   qupath.lib.objects.PathObjects       | Create detection objects
 *   qupath.lib.roi.ROIs                  | Create ellipse ROIs
 *   qupath.lib.regions.ImagePlane        | Define image plane for ROIs
 *   qupath.lib.objects.classes.PathClass | Assign classification labels
 */
import qupath.lib.objects.PathObjects
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.objects.classes.PathClass

// ============================================================
// ⚠️  USER CONFIGURATION — CHANGE THIS PATH
// ============================================================
def tsvPath = "/path/to/your_qupath_cells.tsv"

// ============================================================
// COLOR MAP
// ============================================================
def classColors = [
    "Included"              : getColorRGB(0,   210, 0),
    "Excl_SmallArea"        : getColorRGB(255, 60,  60),
    "Excl_LowDAPI"          : getColorRGB(255, 165, 0),
    "Excl_SmallArea_LowDAPI": getColorRGB(180, 0,   180),
    "Excl_Noise"            : getColorRGB(30,  144, 255),
]

classColors.each { className, color ->
    PathClass.fromString(className, color)
}

// ============================================================
// READ TSV HEADERS
// ============================================================
def detections = []
def lines   = new File(tsvPath).readLines()
def headers = lines[0].split('\t').toList()

def cx_idx    = headers.indexOf('centroid_x')
def cy_idx    = headers.indexOf('centroid_y')
def cls_idx   = headers.indexOf('classification')
def area_idx  = headers.indexOf('area')
def dapi_idx  = headers.indexOf('DAPI')
def ox_idx    = headers.indexOf('offset_x')   // ← reads offset from TSV column
def oy_idx    = headers.indexOf('offset_y')   // ← reads offset from TSV column

println "Headers: ${headers}"
println "Reading ${lines.size() - 1} cells..."

// ============================================================
// FALLBACK: parse offsets from filename if columns missing
// ============================================================
double fallbackOffsetX = 0
double fallbackOffsetY = 0
def matcher = (tsvPath =~ /_x(\d+)_y(\d+)_w\d+_h\d+/)
if (matcher) {
    fallbackOffsetX = matcher[0][1].toDouble()
    fallbackOffsetY = matcher[0][2].toDouble()
}
println "Fallback offsets from filename: x=${fallbackOffsetX}, y=${fallbackOffsetY}"

def plane = ImagePlane.getDefaultPlane()

// ============================================================
// CREATE DETECTIONS WITH OFFSET CORRECTION
// ============================================================
lines[1..-1].each { line ->
    def cols = line.split('\t')
    try {
        double cx   = Double.parseDouble(cols[cx_idx])
        double cy   = Double.parseDouble(cols[cy_idx])
        String cls  = cols[cls_idx]
        double area = Double.parseDouble(cols[area_idx])
        double radius = Math.max(Math.sqrt(area / Math.PI), 3.0)

        // Use offset from TSV column if available, else use filename fallback
        double offsetX = (ox_idx >= 0) ? Double.parseDouble(cols[ox_idx]) : fallbackOffsetX
        double offsetY = (oy_idx >= 0) ? Double.parseDouble(cols[oy_idx]) : fallbackOffsetY

        // KEY FIX: Convert local crop coordinates → global WSI coordinates
        double globalX = cx + offsetX
        double globalY = cy + offsetY

        def roi = ROIs.createEllipseROI(
            globalX - radius,
            globalY - radius,
            radius * 2,
            radius * 2,
            plane
        )

        def pathClass = PathClass.fromString(cls)
        def detection = PathObjects.createDetectionObject(roi, pathClass)

        detection.getMeasurementList().put("Area",       area)
        detection.getMeasurementList().put("DAPI",       Double.parseDouble(cols[dapi_idx]))
        detection.getMeasurementList().put("centroid_x", globalX)
        detection.getMeasurementList().put("centroid_y", globalY)
        detection.getMeasurementList().put("offset_x",   offsetX)
        detection.getMeasurementList().put("offset_y",   offsetY)
        detection.getMeasurementList().close()

        detections.add(detection)

    } catch (Exception e) {
        println "  ⚠️ Skipping row: ${e.message}"
    }
}

// ============================================================
// ADD TO QUPATH & REFRESH
// ============================================================
def hierarchy = getCurrentImageData().getHierarchy()
hierarchy.addPathObjects(detections)
fireHierarchyUpdate()

println "✓ Imported ${detections.size()} cells into QuPath"
