# Triad detection

**Pipeline position:** Step 5 — after cell typing, before functional/survival analysis. The pipeline's core scientific output.

## What this step does

Finds "triads" — an anchor cell type with both of two partner cell types simultaneously nearby within a chosen radius (e.g. a dendritic cell with both a CD4 and a CD8 T cell nearby, forming a priming triad) — across every image in your experiment. Produces per-image triad tables, QC overlay images so you can visually confirm the triads found look right, comparison plots across experiment groups, and a trajectory analysis showing how sensitive the triad count is to your chosen radius.

## Configuring it

- `experiment.name`, `.experiment_groups`, `.image_experiment_group_map`
- `imaging.microns_per_pixel`, `.experiment_group_areas_um2` — needed to convert triad counts into density (triads per area); one value per group works if your images/cores are close to the same size
- `imaging.roi_labels_dir` (optional) — if your images vary meaningfully in area (e.g. whole-slide images rather than uniform TMA cores), point this at your QuPath ROI export to use each image's own measured area instead of one group-wide constant
- `paths.input_dir`, `.output_dir`
- `analysis.triad.radius_um` — the search radius that defines "nearby"
- `analysis.triad.anchor_type`, `.partner_type_1`, `.partner_type_2` — the three cell types that define your triad. Leaving these unset makes the pipeline try every possible 3-way combination of cell types instead, which is thorough but can be slow with more than a handful of cell types
- `analysis.triad.report_radius_um` — a tighter radius used just for the QC plots, if you want to visualize a stricter subset than what you searched with

## What you'll get

- One triad table and one flagged-cells table per image
- A combined summary across all images (`triad_summary.csv`) and experiment-group comparison plots/tables (counts, density, distances)
- QC overlay images per image so you can visually confirm the triads found actually make sense
- A trajectory plot showing triad count vs. distance threshold, useful for justifying your chosen radius

## Things to know

- If you don't set all three of anchor/partner1/partner2, the pipeline searches every 3-way combination of cell types it sees in your data — this can be slow with more than a handful of cell types, and it prints a warning when this fallback kicks in. Worth setting the three types explicitly once you know what you're looking for.
- Already-processed images are automatically skipped on a re-run (adding one new image to a large cohort doesn't reprocess everything) — but if you edit an image's underlying cell data after it's already been processed, delete that image's two output files yourself to force it to be redone; there's no automatic way to detect the input changed.
- If per-image area comes from `roi_labels_dir` and some images don't have a matching ROI file, their area is imputed from the group average (or dropped from the total if there's no group average either) — the run log prints exactly which method was used for every image, so you can trace where your density numbers came from.
