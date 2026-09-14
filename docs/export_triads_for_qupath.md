# Exporting triads for QuPath visualization

## What this is for

Converts triad detection's output into a format QuPath can import, so you can visually check detected triads against the actual tissue image.

## Usage

```bash
python export_triads_for_qupath.py --input-dir /path/to/triad_analysis_output

python export_triads_for_qupath.py \
    --input-dir /path/to/triad_analysis_output \
    --output-dir /path/to/triad_qupath_exports \
    --anchor-cell-type "Dendritic cells" \
    --partner1-cell-type "CD4 T cells" \
    --partner2-cell-type "CD8 T cells"
```

Run this after triad detection, then use the triad-visualization QuPath script to render the result.

## Things to know

- If a file for one image is corrupted or unreadable, it's skipped and reported rather than stopping the whole export — check the console output for any warning lines before assuming every image was exported.
