# QuPath: cell-type visualization (dot plot)

**Where this fits:** a manual, optional QuPath step for visually reviewing cell types after cell typing has run, directly inside your QuPath project. This is the original QuPath script the pipeline's Python visualization tool is based on.

## What it does

Draws each cell as a colored, sized dot at its location in QuPath, colored by cell type, with marker levels and other details attached so you can inspect any cell by clicking on it.

## How to run it

Open the relevant QuPath project, edit the file path at the top of the script to point at the data file for the image you're viewing, then run it from QuPath's Script Editor.

## Things to know

- The cell-type colors built into this script were chosen for the older LILRB2 mouse study and don't match this CRC study's cell types — running it as-is on CRC data renders every cell type in a generic default color rather than the intended palette, with no warning that this happened. Use the Python visualization tool instead for the CRC dataset, or update this script with CRC's own cell-type names first.
- You need to edit the file path at the top of the script before every run.
