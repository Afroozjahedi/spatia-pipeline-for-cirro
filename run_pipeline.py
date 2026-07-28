#!/usr/bin/env python3
"""
run_pipeline.py
===============
SPATIA pipeline orchestrator.

Chains all enabled steps in order:
    0. tif_conversion  → raw QPTIFF/OME-TIFF → standard .tif + metadata JSON
    1. roi_masking     → crop + mask ROIs/TMA cores using QuPath-exported masks
    2. segmentation    → Mesmer/Cellpose cell segmentation → *_mesmer_result.csv
    3. preprocessing   → filter, normalise, noise-remove, save h5ad
    4. cell_typing     → GMM thresholds, assign cell types (auto or semi-auto)
    5. triads          → detect DC–CD4–CD8 triads, compute density, plot
    6. functional      → in-triad vs out-of-triad marker expression (reads triads output)
    7. survival        → Kaplan-Meier OS/DFS by triad density (reads triads output)

Steps 0-2 are new upstream stages (ROI selection is a separate QuPath/Groovy
step, run before this script -- see qupath_scripts/). Before this, the
pipeline only had a runnable path starting at preprocessing; segmentation_results_dir
had to already exist. Steps 0-2 are OFF by default (see tif_conversion.enabled /
roi_masking.enabled / segmentation.enabled in the config) so existing configs
that already have a segmentation_results_dir keep working unchanged.

Usage
-----
    python run_pipeline.py --config experiments/PirB_D14.yaml

    # Run specific steps only
    python run_pipeline.py --config experiments/PirB_D14.yaml --steps preprocessing cell_typing

    # Skip a step (e.g. preprocessing already done)
    python run_pipeline.py --config experiments/PirB_D14.yaml --skip preprocessing

Exit codes
----------
    0  all requested steps completed successfully
    1  one or more steps failed
    2  config file not found or invalid
"""

import argparse
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import yaml


# ─────────────────────────────────────────────────────────────────────────────
# TEE LOGGING — mirrors stdout/stderr to a log file in output_dir
# ─────────────────────────────────────────────────────────────────────────────

class _Tee:
    """Write to multiple streams simultaneously."""
    def __init__(self, *streams):
        self._streams = streams

    def write(self, data):
        for s in self._streams:
            s.write(data)
        if data.endswith("\n"):
            self.flush()

    def flush(self):
        for s in self._streams:
            s.flush()

    def fileno(self):
        return self._streams[0].fileno()


def _start_logging(output_dir: str, run_ts: str) -> Path:
    """
    Open a log file at {output_dir}/logs/pipeline_{run_ts}.log and
    redirect stdout + stderr so every print() also lands in the file.
    Returns the log path.
    """
    log_dir = Path(output_dir) / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"pipeline_{run_ts}.log"

    log_file = open(log_path, "w", buffering=1)          # line-buffered
    sys.stdout = _Tee(sys.__stdout__, log_file)
    sys.stderr = _Tee(sys.__stderr__, log_file)
    return log_path


# ─────────────────────────────────────────────────────────────────────────────
# STEP REGISTRY
# Keys must match --steps / --skip argument values.
# ─────────────────────────────────────────────────────────────────────────────

ALL_STEPS = [
    "tif_conversion", "roi_masking", "segmentation",
    "preprocessing", "cell_typing", "triads", "functional", "survival",
]


def _import_steps(steps_to_run: list):
    """
    Import only the modules needed for the steps that will actually run.
    This avoids crashing on missing dependencies (e.g. spacec) when
    those steps are being skipped.
    """
    _MODULE_MAP = {
        "tif_conversion": ("spatia.analysis.tif_conversion", "run_tif_conversion"),
        "roi_masking":   ("spatia.analysis.roi_masking",    "run_roi_masking"),
        "segmentation":  ("spatia.analysis.segmentation",   "run_segmentation"),
        "preprocessing": ("spatia.analysis.preprocessing", "run_preprocessing"),
        "cell_typing":   ("spatia.analysis.cell_typing",   "run_cell_typing"),
        "triads":        ("spatia.analysis.triads",        "run_triad_analysis"),
        "functional":    ("spatia.analysis.functional",    "run_functional_analysis"),
        "survival":      ("spatia.analysis.survival",      "run_survival_analysis"),
    }
    step_fns = {}
    for step in steps_to_run:
        module_path, fn_name = _MODULE_MAP[step]
        import importlib
        mod = importlib.import_module(module_path)
        step_fns[step] = getattr(mod, fn_name)
    return step_fns


def _import_validator():
    from spatia.analysis.validation import validate_step
    return validate_step


# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _banner(text: str, char: str = "=", width: int = 72) -> str:
    border = char * width
    return f"\n{border}\n{text}\n{border}"


def _hms(seconds: float) -> str:
    h, rem = divmod(int(seconds), 3600)
    m, s   = divmod(rem, 60)
    if h:
        return f"{h}h {m}m {s}s"
    if m:
        return f"{m}m {s}s"
    return f"{s}s"


def _step_enabled_in_config(step: str, cfg: dict) -> bool:
    """
    Check whether a step is switched on inside the config.
    preprocessing and cell_typing are always on (no enable flag).
    tif_conversion / roi_masking / segmentation are OFF by default (they're
    new upstream steps -- existing configs that already point straight at a
    populated segmentation_results_dir must keep working unchanged) and are
    checked at the top level (cfg['tif_conversion']['enabled'], etc.), not
    under cfg['analysis'] like triads/functional/survival/tls.
    All other steps check cfg['analysis'][step]['enabled'].
    Also tries the singular form (e.g. "triad" for step "triads").
    """
    if step in ("preprocessing", "cell_typing"):
        return True
    if step in ("tif_conversion", "roi_masking", "segmentation"):
        return cfg.get(step, {}).get("enabled", False)
    analysis = cfg.get("analysis", {})
    return (
        analysis.get(step, {}).get("enabled", False)
        or analysis.get(step.rstrip("s"), {}).get("enabled", False)
    )


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="SPATIA pipeline orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--config", "-c",
        required=True,
        help="Path to experiment YAML config file",
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        choices=ALL_STEPS,
        default=None,
        metavar="STEP",
        help=(
            "Steps to run (default: all enabled steps). "
            f"Choices: {', '.join(ALL_STEPS)}"
        ),
    )
    parser.add_argument(
        "--skip",
        nargs="+",
        choices=ALL_STEPS,
        default=[],
        metavar="STEP",
        help="Steps to skip even if enabled in config.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the steps that would run without executing them.",
    )
    args = parser.parse_args()

    # ── Load config ───────────────────────────────────────────────────────
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"ERROR: config file not found: {config_path}", file=sys.stderr)
        sys.exit(2)

    try:
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
    except yaml.YAMLError as exc:
        print(f"ERROR: invalid YAML in {config_path}:\n{exc}", file=sys.stderr)
        sys.exit(2)

    exp_name = cfg.get("experiment", {}).get("name", config_path.stem)

    # ── Start logging ─────────────────────────────────────────────────────
    run_ts  = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = cfg.get("paths", {}).get("output_dir", ".")
    log_path = _start_logging(out_dir, run_ts)
    print(f"Log file    : {log_path}")

    # ── Determine which steps to run ──────────────────────────────────────
    requested = args.steps if args.steps else ALL_STEPS
    skip_set  = set(args.skip)

    steps_to_run = [
        s for s in requested
        if s not in skip_set and _step_enabled_in_config(s, cfg)
    ]

    steps_skipped = [
        s for s in requested
        if s in skip_set or not _step_enabled_in_config(s, cfg)
    ]

    # ── Print plan ────────────────────────────────────────────────────────
    print(_banner(f"SPATIA PIPELINE  |  {exp_name}"))
    print(f"Config  : {config_path.resolve()}")
    print(f"Started : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nSteps to run : {steps_to_run or '(none)'}")
    if steps_skipped:
        print(f"Steps skipped: {steps_skipped}")

    if args.dry_run:
        print("\n[DRY RUN] No steps executed.")
        sys.exit(0)

    if not steps_to_run:
        print("\nNothing to do — check your config or --steps argument.")
        sys.exit(0)

    # ── Import step functions ─────────────────────────────────────────────
    try:
        step_fns  = _import_steps(steps_to_run)
        validator = _import_validator()
    except ImportError as exc:
        print(f"\nERROR: could not import spatia modules.\n{exc}", file=sys.stderr)
        print("Make sure the spatia package is on your PYTHONPATH.", file=sys.stderr)
        sys.exit(1)

    # ── Execute steps ─────────────────────────────────────────────────────
    pipeline_start = time.time()
    results        = {}
    failed_steps   = []

    for step in steps_to_run:
        print(_banner(f"STEP: {step.upper()}", char="-"))
        t0 = time.time()

        try:
            result = step_fns[step](cfg)
            elapsed = time.time() - t0
            results[step] = result
            print(f"\n✓  {step} completed in {_hms(elapsed)}")

        except Exception as exc:
            elapsed = time.time() - t0
            failed_steps.append(step)
            print(f"\n✗  {step} FAILED after {_hms(elapsed)}")
            print(f"   {type(exc).__name__}: {exc}")
            traceback.print_exc()

            remaining = [s for s in steps_to_run if s not in failed_steps and s != step]
            skip_hint = f"--skip {' '.join(remaining)}" if remaining else ""
            print(f"\nPipeline halted. Fix the error above, then re-run with:\n"
                  f"  python run_pipeline.py --config {args.config} "
                  f"--skip {step} {skip_hint}".strip())
            break

        # ── Validate outputs ──────────────────────────────────────────────
        print(f"\n  Validating {step} outputs…")
        passed, val_errors = validator(step, cfg)
        if passed:
            print(f"  ✓  Validation passed")
        else:
            print(f"\n  ✗  Validation FAILED for step '{step}':")
            for err in val_errors:
                print(f"     • {err}")
            failed_steps.append(step)
            remaining = [s for s in steps_to_run if s not in results and s != step]
            skip_hint = " ".join(remaining)
            print(f"\nPipeline halted at validation. Fix the issues above, then re-run with:\n"
                  f"  python run_pipeline.py --config {args.config} "
                  f"--skip {step}{' --skip ' + skip_hint if skip_hint else ''}".strip())
            break

    # ── Final summary ─────────────────────────────────────────────────────
    total_elapsed = time.time() - pipeline_start
    print(_banner("PIPELINE SUMMARY"))
    print(f"Experiment : {exp_name}")
    print(f"Total time : {_hms(total_elapsed)}")
    print(f"Finished   : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    for step in steps_to_run:
        if step in failed_steps:
            status = "✗  FAILED"
        elif step in results:
            status = "✓  OK"
        else:
            status = "–  NOT REACHED"
        print(f"  {status:<12} {step}")

    if steps_skipped:
        for step in steps_skipped:
            print(f"  ⏭   SKIPPED    {step}")

    if failed_steps:
        sys.exit(1)

    print(f"\nAll steps completed successfully.")
    sys.exit(0)


if __name__ == "__main__":
    main()
