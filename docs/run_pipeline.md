# `run_pipeline.py`

**Pipeline position:** The orchestrator itself, not a step. Runs all 8 steps in order: `tif_conversion` → `roi_masking` → `segmentation` → `preprocessing` → `cell_typing` → `triads` → `functional` → `survival`.

## Purpose

The deterministic SPATIA pipeline orchestrator — the CLI entry point that chains all 8 pipeline steps in order, with per-step logging, output validation, and clean failure reporting. This is the script deployed to Cirro/HPC for the Nature Protocols submission, and is deliberately kept free of any LLM dependency (see `run_agentic.py`, which wraps this module rather than duplicating it).

Steps: `tif_conversion` → `roi_masking` → `segmentation` → `preprocessing` → `cell_typing` → `triads` → `functional` → `survival`. The first three are newer upstream steps, off by default (`enabled: false`) so existing configs pointing straight at a pre-populated `segmentation_results_dir` keep working unchanged.

## Usage

```bash
python run_pipeline.py --config experiments/PirB_D14.yaml
python run_pipeline.py --config experiments/PirB_D14.yaml --steps preprocessing cell_typing
python run_pipeline.py --config experiments/PirB_D14.yaml --skip preprocessing
python run_pipeline.py --config experiments/PirB_D14.yaml --dry-run
```

Exit codes: `0` all requested steps succeeded, `1` a step failed, `2` config file missing/invalid.

## Key functions

| Function | Role |
|---|---|
| `_Tee` | Small stream multiplexer so every `print()` goes to both the console and a log file. |
| `_start_logging(output_dir, run_ts)` | Redirects `sys.stdout`/`sys.stderr` through `_Tee` to `{output_dir}/logs/pipeline_{run_ts}.log`. |
| `ALL_STEPS` | The canonical ordered list of 8 step names — the source of truth `--steps`/`--skip` validate against. |
| `_import_steps(steps_to_run)` | Lazily imports only the `spatia.analysis.*` modules needed for the requested steps, so e.g. a missing `spacec` install doesn't block a run that only needs `triads`/`functional`/`survival`. |
| `_step_enabled_in_config(step, cfg)` | Determines whether a step should run: `preprocessing`/`cell_typing` are always on; `tif_conversion`/`roi_masking`/`segmentation` check `cfg[step]['enabled']` (top-level, off by default); everything else checks `cfg['analysis'][step]['enabled']`, also trying the singular form (`triad` for `triads`). |
| `main()` | Argparse CLI → loads config → starts logging → resolves the step list (`requested - skip`, filtered by `_step_enabled_in_config`) → imports step functions + `validate_step` → runs each step, calling `validate_step` after every successful step and halting on the first failure (either an exception or a failed validation) → prints a per-step summary table. |

## Config keys referenced

- `experiment.name`
- `paths.output_dir`
- `{step}.enabled` for `tif_conversion`/`roi_masking`/`segmentation`
- `analysis.{step}.enabled` (or singular form) for `triads`/`functional`/`survival`

## Notes / risks

- **A single step failure halts the entire remaining pipeline (confidence: high, by design — worth confirming this is the intended behavior, not just the current one).** Both an exception and a failed validation `break` out of the loop entirely rather than skipping just the failed step. For a long multi-image run this means, e.g., a `survival` step failing because the patient annotation file has a typo throws away the fact that `functional` already succeeded in the same invocation — though the outputs on disk from earlier steps are preserved, only the summary/exit-code treats the whole run as failed. Reasonable for a linear dependency chain (each step reads the last one's output), but worth being deliberate about since `functional` and `survival` are actually independent siblings, not a strict chain (both just read `triads` output) — they could in principle run even if the other fails.
- **Two normalization schemes for step names coexist (`triads` vs `triad`, confidence: medium).** `_step_enabled_in_config` tries both `analysis[step]` and `analysis[step.rstrip('s')]` to handle the plural/singular config-key mismatch. This is a reasonable defensive shim, but it does mean two different config authors could write `analysis.triad.enabled` and `analysis.triads.enabled` and get different-looking-but-equivalent configs — worth standardizing on one in new configs even though both are accepted.
- **No parallelism across independent steps (confidence: high, structural observation, not a defect).** Steps run strictly sequentially even where the dependency graph would allow concurrency (e.g. `functional` and `survival` could both start once `triads` finishes). Not a problem at current scale, but worth knowing before assuming the pipeline can't be sped up without restructuring `triads`' outputs.
- **This is the natural place a future "agentic step selection" would plug in (confidence: medium, forward-looking observation given the project's stated goal).** Since your project instructions mention wanting an *agentic pipeline*, note that `_step_enabled_in_config` + `ALL_STEPS` is currently a static, config-driven gate — there's no mechanism here (yet) for an agent to decide dynamically which steps to run based on prior step outputs (e.g. skip `survival` if too few patients have annotation data). `run_agentic.py` adds agentic *behavior within* one step (`cell_typing`), not agentic *orchestration* of the pipeline as a whole — worth deciding whether that gap matters for what "agentic pipeline" means for this project.
