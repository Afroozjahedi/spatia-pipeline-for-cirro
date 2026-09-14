# Cell typing

**Pipeline position:** Step 4 — after preprocessing, before triad/functional/survival analysis.

## What this step does

Assigns a cell type to every cell in your data, based on which markers it's positive or negative for. There are two ways to run it:

- **Automatic** — the pipeline decides positivity for every marker itself, then matches each cell against the rules in your `cell_type_definitions.yaml` file (e.g. "CD3+ and CD8+ and CD4- → CD8 T cell"). No manual review needed.
- **Semi-automatic** — same marker-positivity step, but instead of rule-matching, cells are clustered (so similar cells group together) and the run stops so you can look at the clusters and assign a cell-type label to each one by hand. Re-running after you've filled in the labels finishes the job.

## Deciding whether a cell is positive for a marker

For each marker, the pipeline looks at the distribution of intensities across all cells and tries to split it into a "negative" group and a "positive" group. There are two settings for how strict that split is (`gmm.threshold_mode` in your config):

- **`std_multiplier`** (the long-standing default) — draws one cutoff line above the negative population; anything past it counts as positive. Simple and fast, but it only looks at the negative population's own spread.
- **`posterior`** — for each cell, asks "how likely is this cell to belong to the positive population, given both populations' shapes?" and uses `confidence_level` (default 0.8, i.e. 80% likely) as the bar. This is the more statistically complete option, especially when the positive population is more spread out than the negative one.

**Something to flag before this feeds the paper:** the CRC pipeline's current validated accuracy number (24.5% overall / 26.8% by core category, from comparing automatic calls against expert-labeled cells) was measured using `std_multiplier` mode. `crc_tma_full_pipeline.yaml` now uses `posterior` mode for the real production run — that combination hasn't been re-checked against the expert labels yet. If the manuscript states "posterior thresholding" as the method, it's worth re-running the accuracy comparison under posterior mode first so the reported number actually matches what was run.

## Markers with a lot of cells stuck at one value ("zero-inflated" markers)

Some markers have a large share of cells sitting at exactly the same low value — a real detection floor, not noise — rather than a smooth curve. Left alone, this confuses the negative/positive split (the fitted "negative" population can collapse onto that floor and the resulting cutoff lands too low, calling almost everything positive).

The fix is to apply a transform to reshape the marker's values before the positive/negative split is fit (`gmm.transform`, and `gmm.per_marker_transform_overrides` for applying it only to specific markers rather than the whole panel). Two transform options are available:

- **arcsinh** — a simple, fixed transformation. Same result every time you run it, which makes it easy to describe and reproduce in a methods section. **This is the current recommendation** for markers that need a transform.
- **Yeo-Johnson** — adapts itself to each marker's own data. This can fit a marker's shape better, but because it refits itself from the data each run, its exact behavior can shift slightly if you rerun it on a different subset of cells — worth keeping in mind for reproducibility if you switch to it.

The current list of which 14 markers get a transform (out of ~55) was built by hand, by looking at each marker's distribution — it isn't detected automatically, and applying a transform to a marker that doesn't need one can make its threshold worse, not better, so don't apply it blanket-style across the whole panel.

If you want an automatic version of this (the pipeline detects zero-inflated markers itself and picks a transform), that's a feature worth discussing before building — it would need a rule for how much "stuck at one value" is enough to trigger a transform, and could reduce reproducibility if the pick changes run to run.

## Configuring it

- `cell_typing.mode` — `"automatic"` or `"semi_automatic"`
- `cell_typing.markers.panel` — the list of markers to type cells on
- `cell_typing.markers.gating_only` — markers used only to gate cells in/out beforehand (CD45, by default), not typed themselves
- `cell_typing.markers.column_map` — only needed if your data's column names don't match the short marker names used elsewhere in the config (e.g. your CODEX export names columns like `"CD45 - hematopoietic cells (C14)"` instead of just `"CD45"`) — maps the short name to the real column name
- `cell_typing.gmm.threshold_mode` / `.confidence_level` — see above
- `cell_typing.gmm.transform` / `.per_marker_transform_overrides` — see above
- `cell_typing.gmm.cd45_std_multiplier` / `.default_std_multiplier` / `.per_marker_overrides` — how strict the cutoff is per marker, when using `std_multiplier` mode
- `cell_typing.gmm.skip_cd45_gate` — set `true` to run on all cells instead of gating to immune (CD45+) cells first; needed for panels covering both immune and non-immune cell types (like this CRC panel)
- `cell_typing.clustering.n_neighbors` / `.leiden_resolution` — only used in semi-automatic mode, control how clusters are formed
- `cell_typing.cell_type_definitions_file` — required for automatic mode; the rules mapping marker combinations to cell types
- `cell_typing.cluster_labels_file` — semi-automatic mode; leave as `null` for the first run, fill in and re-run to finish

**Worth double-checking:** automatic mode and semi-automatic mode have historically used different strictness settings for the CD45 gate (semi-automatic tighter, automatic more permissive). Both are explicit config values now — just confirm each of your configs sets what you actually intend, since the two modes aren't directly comparable unless this matches.

## What you'll get

- A cell-typed h5ad file with a `cell_type` column (and per-marker positive/negative + intensity-level columns) added to every cell
- `marker_thresholds.csv` — the cutoff used for every marker, which mode and transform were used, so you can see exactly what was applied
- `cell_type_counts.csv` and a bar chart of cell type counts
- A breakdown of cell type proportions by experiment group, both as counts and as a chart
- A full diagnostic plot set: one distribution plot per marker showing where its cutoff landed, a positivity heatmap and an expression heatmap by cell type, UMAPs colored by cell type/group/tissue, a dot plot, per-marker UMAPs, and a report bundling everything as one HTML page and one PDF

## Things to know

- **Ties between cell types go to whichever is listed first** in `cell_type_definitions.yaml`, if two types score equally on the same cell. If you want two types to be mutually exclusive, make sure their rules don't allow an equal-score tie in the first place, or order the file deliberately.
- **Semi-automatic mode's "stop" isn't an error.** After clustering, it deliberately pauses so you can review the clusters and fill in labels — re-run with `cluster_labels_file` pointing at your filled-in file to finish.
- If a marker in your panel doesn't show up in your data, it's most often a name-mismatch (see `column_map` above) rather than the marker being genuinely missing — worth checking that before assuming a marker isn't in the data at all.
