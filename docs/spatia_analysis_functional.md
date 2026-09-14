# Functional marker analysis

**Pipeline position:** Step 6a (optional) — runs after triad detection.

## What this step does

For each cell type in your triads (anchor, partner 1, partner 2), compares marker expression two ways: cells that are part of a triad vs. cells that aren't, and — within triad cells only — one experiment group vs. another. Useful for asking, e.g., "do T cells in a triad look more activated than T cells that aren't?" Statistical significance is assessed with Mann-Whitney U tests, corrected for the number of comparisons run.

## Configuring it

- `experiment.name`, `.experiment_groups`, `.image_experiment_group_map`
- `paths.output_dir`
- `analysis.triad.anchor_type` / `.partner_type_1` / `.partner_type_2` — must be set (this step doesn't have triad detection's "try every combination" fallback — leaving these unset means it finds nothing)
- `analysis.functional.enabled` — must be `true` to run this step at all
- `analysis.functional.report_radius_um` (default `20.0`) — can be tighter than the radius used for triad detection itself
- `analysis.functional.markers` — which markers to compare, as `{display name: actual column name}`
- `analysis.functional.markers_anchor` / `.markers_partner1` / `.markers_partner2` (optional) — use a different marker set for a specific role, e.g. exhaustion markers on a T-cell partner but maturation markers on a dendritic-cell anchor. A role without its own override uses the shared `markers` list

## What you'll get

- A violin plot per cell type per comparison (in-triad vs. not; group A vs. group B)
- `functional_marker_summary.csv` with a Bonferroni-corrected p-value for every comparison

## Things to know

- Marker columns need to actually be present in your per-cell data — that depends on marker intensities having been kept during an earlier data-prep step. If a configured marker is missing, this step reports which one and skips just that comparison rather than failing the whole run — worth checking the summary CSV's row count matches what you expect if some markers were misconfigured.
- The multiple-testing correction is applied across every comparison in the whole run at once, not separately per marker or per cell type — adding more markers to your config makes it statistically harder for any single marker to come out significant. Worth keeping in mind when deciding how many markers to include at once.
- Unlike triad detection, this step doesn't skip re-processing on a re-run — it always re-pools and re-tests everything, since its output is a single statistic across the whole cohort rather than something computable image-by-image. That's intentional, and it's fast either way.
