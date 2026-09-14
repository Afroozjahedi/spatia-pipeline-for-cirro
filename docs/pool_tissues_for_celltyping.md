# Pooling multiple tissues before cell typing

## What this is for

If your cohort has more than one tissue (different timepoints, patients, or TMA blocks), preprocessing keeps each tissue's cells separate — it doesn't produce one combined file across your whole cohort. Cell typing, though, works on a single file. This script combines every tissue's processed data into one file, so cell-type thresholds are calibrated jointly across your whole cohort instead of separately per tissue.

**Whether joint typing is the right choice for your cohort is a real decision, not just a technical step.** Pooling assumes your tissues are comparable enough in staining/imaging that one shared threshold per marker makes sense. If some tissues have real batch effects, pooling could blur that rather than fix it — per-tissue typing, or a batch-correction step, might be more appropriate instead. This script just removes the "there's no way to do it" blocker; it doesn't make that call for you.

## Usage

```bash
python pool_tissues_for_celltyping.py \
    --input-dir /path/to/combined_processed_data/individual_processed_data \
    --output-file /path/to/combined_processed_data/cohort_pooled.h5ad
```

Then point cell typing at the result:

```yaml
cell_typing:
  input_file: /path/to/combined_processed_data/cohort_pooled.h5ad
```

## Options

- `--join inner` (default) — keeps only markers present in every tissue; safest choice, recommended
- `--join outer` — keeps every marker from every tissue, leaving gaps for tissues that don't have a given one; only use this if you're prepared to handle those gaps yourself before cell typing, since cell typing hasn't been tested against data with gaps in it

## What you'll get

One combined file covering your whole cohort, with a `tissue_id` column added so you can still check afterward whether cell-type composition or QC differs by tissue even after joint typing.

## Things to know

- If your tissues have different marker panels, the script tells you exactly which marker(s) are missing from which tissue(s) before combining — worth reading that message, since it affects which markers you'll be able to type on.
- Running this on just one tissue's data is fine — it just makes a copy with the `tissue_id` column added, so you can point cell typing at the same kind of file regardless of how many tissues you have.
