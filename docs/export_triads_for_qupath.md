# Exporting triads for QuPath visualization

## What this is for

Converts triad detection's output into a format QuPath can import, so you can visually check detected triads against the actual tissue image.

## Automatic (recommended)

As of 2026-09-23 this runs automatically at the end of triad detection if you set `analysis.triad.qupath_export.enabled: true` in your config — no separate command needed. See [Triad detection](spatia_analysis_triads.md#configuring-it). Anchor/partner1/partner2 cell-type labels default to whatever you set for `analysis.triad.anchor_type` / `.partner_type_1` / `.partner_type_2`, so in most cases the config block is just:

```yaml
analysis:
  triad:
    qupath_export:
      enabled: true
```

## Manual (standalone script)

The script below still works unchanged — useful for re-exporting an older run's output without re-running triad detection, or for a one-off export with different cell-type labels than the run used.

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
