# QuPath: cell-boundary visualization

**Where this fits:** a manual, optional QuPath step for visually reviewing segmentation and cell typing together, after both have run. Not part of the automatic pipeline.

## What it does

Draws each cell in QuPath using its actual segmented shape (not just a dot), colored by cell type — useful for checking segmentation quality alongside cell-type assignment. If a cell's shape data is missing, it draws a small circle at the cell's location instead of skipping it outright.

## How to run it

Open the relevant QuPath project, edit the file path at the top of the script to point at the data file for the image you're viewing, then run it from QuPath's Script Editor.

## Things to know

- You need to edit that file path yourself before every run — the script clearly marks where to do it.
- This script doesn't apply cell-type colors the way the plain dot-visualization script does — cells show up in whatever default color QuPath assigns, not a deliberately chosen palette. Fine if you're just checking cell shapes, less useful if you want colors to match other visualizations.
- If a lot of cells get skipped, the console only shows the first few reasons why — worth checking your source data file directly if the skipped count looks high.
