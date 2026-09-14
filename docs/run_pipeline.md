# Running the pipeline

The single command that runs the whole SPATIA pipeline end to end, in order: format conversion → ROI masking → segmentation → preprocessing → cell typing → triads → functional → survival. The first three steps are off by default, so an existing config that already starts from segmented data keeps working without changes.

## Usage

```bash
python run_pipeline.py --config experiments/your_experiment.yaml
python run_pipeline.py --config experiments/your_experiment.yaml --steps preprocessing cell_typing
python run_pipeline.py --config experiments/your_experiment.yaml --skip preprocessing
python run_pipeline.py --config experiments/your_experiment.yaml --dry-run
```

`--steps` runs only the named steps; `--skip` runs everything except the named steps; `--dry-run` shows what would run without actually running it.

Exit codes: `0` = everything requested succeeded, `1` = a step failed, `2` = the config file is missing or invalid.

## What happens on a failure

Every step is checked automatically right after it runs, and the whole pipeline stops at the first step that fails or fails its check — it does not skip ahead and try later steps. Anything already produced by earlier successful steps stays on disk; only the overall run is marked failed. Worth knowing this if, say, a problem in `survival` (a typo in your patient annotation file) stops it from running even though `functional` already succeeded in the same invocation — the two are independent of each other, but the pipeline currently treats a later failure as stopping the whole run.

## Things to know

- All 8 steps run one after another, never in parallel, even where two steps don't actually depend on each other (`functional` and `survival`, for instance, both just read `triads`' output). Not a problem at typical dataset sizes, but if a run ever feels slower than expected, this is why, not a sign something's broken.
- Config files can spell a step's settings either singular or plural (`analysis.triad.enabled` vs. `analysis.triads.enabled`) — both work, but it's worth picking one convention for your own configs so they're easy to compare against each other later.
