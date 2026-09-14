# Preparing a new dataset for triad detection

## What this is for

Triad detection needs your per-cell data in one specific format (one CSV per experiment-group/sample, with specific column names). This script converts your raw per-cell data (a CSV/TSV, or an `.h5ad` file) into that format for any new dataset — you don't need a dataset-specific script, just your real column names.

## Usage

Step 1 — always look first:

```bash
python prepare_matched_cells.py --inspect --input /path/to/data.csv
```

This prints your file's actual columns (and some likely guesses for which ones are your cell type, coordinate, group, and sample columns) without writing anything — use it to confirm the real column names before converting.

Step 2 — convert, using the column names confirmed in step 1:

```bash
python prepare_matched_cells.py \
    --input  /path/to/data.csv \
    --output data/my_dataset/matched_cells \
    --cell-type-col cell_type_name \
    --x-col X_centroid --y-col Y_centroid \
    --experiment-group-col tissue_type \
    --sample-col unique_region \
    --experiment-group-map '{"1": "CLR", "2": "DII"}' \
    --cell-type-merge-map '{"^CD4\\+ T cells.*": "CD4+ T cells"}'
```

Works the same way for `.h5ad` files.

## Useful options

- `--experiment-group-map` — relabel raw group codes (e.g. `1`/`2`) to readable names (e.g. `CLR`/`DII`); optional
- `--cell-type-merge-map` — merge similar cell-type labels into one (e.g. collapse several CD4+ T-cell subtypes into a single label); optional, applied as find-and-replace patterns
- `--force` — redo a sample even if its output file already exists (by default, already-converted samples are skipped so you can re-run safely)

## What you'll get

One CSV per experiment-group/sample combination, ready to point triad detection's input directory at directly.

## Things to know

- If any required column ends up with missing values after conversion, you'll get a warning rather than a hard failure — a large number of missing values there is usually a sign the wrong column was picked, worth checking before trusting your triad counts.
- If you use `--cell-type-merge-map`, the original (pre-merge) label isn't kept anywhere separately — if you'll want to audit which original labels got merged together later, keep a copy of your source data with the original column intact before running this.
- This hasn't yet been used on a dataset where cell coordinates live somewhere other than the standard table format some `.h5ad` files store coordinates separately from the main table. `--inspect` will show you if that's the case for your data; flag it if so, since it needs a small adjustment to handle correctly.
