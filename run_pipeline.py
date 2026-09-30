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
import csv
import os
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
# SAMPLESHEET OVERRIDE — optional, added 2026-09-24 for Cirro packaging
# ─────────────────────────────────────────────────────────────────────────────

# Column-name aliases accepted in a samplesheet CSV. Cirro's own native
# samplesheet convention (confirmed 2026-09-30 against a real samplesheet
# Afrouz uploaded through Cirro's "Samples" feature) is sample,file,group --
# NOT the sample_id,image_path,group names this function originally required
# (those were this repo's own invented names, written before we'd seen what
# Cirro actually produces). Both are accepted below so nothing that already
# depends on the old names breaks; SAMPLE_COL/FILE_COL below is the
# preferred/canonical pair going forward -- see derive_samplesheet.py, which
# now writes sample,file,group directly.
_SAMPLE_COL_ALIASES = ("sample", "sample_id")
_FILE_COL_ALIASES = ("file", "image_path")


def _resolve_samplesheet_columns(fieldnames: list) -> tuple:
    """
    Picks whichever alias of the sample-id and file-path columns is present
    in this CSV's actual header, preferring Cirro's own names (sample, file)
    over this repo's older names (sample_id, image_path) when a file somehow
    has both. Returns (sample_col, file_col). Raises via sys.exit(2) with a
    clear message listing both accepted spellings if neither alias for a
    required column is present, or if "group" itself is missing.
    """
    fieldset = set(fieldnames or [])

    sample_col = next((c for c in _SAMPLE_COL_ALIASES if c in fieldset), None)
    file_col = next((c for c in _FILE_COL_ALIASES if c in fieldset), None)
    missing = []
    if sample_col is None:
        missing.append("/".join(_SAMPLE_COL_ALIASES))
    if file_col is None:
        missing.append("/".join(_FILE_COL_ALIASES))
    if "group" not in fieldset:
        missing.append("group")

    if missing:
        print(
            f"ERROR: samplesheet is missing required column(s): {missing} "
            f"(sample id column accepts either {_SAMPLE_COL_ALIASES}, "
            f"file/path column accepts either {_FILE_COL_ALIASES})",
            file=sys.stderr,
        )
        sys.exit(2)

    return sample_col, file_col


def _common_dir(paths: list) -> str:
    """
    Return the common ancestor directory of a list of file paths.

    Deliberately does NOT use os.path.commonpath or pathlib.Path: both
    treat consecutive slashes as redundant and collapse them, which turns
    's3://bucket/a/b.tif' into 's3:/bucket/a/b.tif' (single slash) -- not a
    valid S3 URI, and the exact bug that broke the 2026-09-30 run (masked_roi_dir
    resolved to 's3:/...', segmentation.py's S3FileSystem then couldn't find
    it). Plain string-splitting on '/' preserves the empty component between
    the two slashes after 's3:', so it round-trips correctly for both S3 URIs
    and ordinary POSIX paths.
    """
    def _dirname(path: str) -> str:
        idx = path.rfind("/")
        return path[:idx] if idx >= 0 else path

    split = [_dirname(p).split("/") for p in paths]
    common = []
    for parts in zip(*split):
        if len(set(parts)) == 1:
            common.append(parts[0])
        else:
            break
    if not common:
        raise ValueError(f"No common directory among: {paths}")
    return "/".join(common)


def _apply_samplesheet(cfg: dict, samplesheet_path: Path) -> dict:
    """
    Overrides cfg["experiment"]["image_experiment_group_map"] and
    cfg["paths"]["masked_roi_dir"] with sample metadata read from a
    samplesheet CSV, instead of the hand-written values in the YAML.
    Everything else in cfg (all pipeline parameters -- thresholds, which
    steps run, marker panel, etc.) is untouched.

    Column names: accepts either Cirro's own convention (sample, file,
    group) or this repo's earlier names (sample_id, image_path, group) --
    see _resolve_samplesheet_columns. Whichever pair is present in the
    actual CSV header is used; you do not need to pick one.

    This keeps the pipeline's existing single-batched-call design (Q10):
    it does not change how many times the pipeline runs, only where the
    sample list comes from. segmentation.py still walks ONE masked_roi_dir
    and picks up each sample's folder name -- so every file/image_path in
    the samplesheet must share the same parent directory. If your samples
    do not share one parent, this override does not apply cleanly -- run
    without --samplesheet and use the YAML's own image_experiment_group_map
    instead, or extend this function, rather than relying on it silently.
    """
    if not samplesheet_path.exists():
        print(f"ERROR: samplesheet not found: {samplesheet_path}", file=sys.stderr)
        sys.exit(2)

    rows = []
    with open(samplesheet_path, newline="") as fh:
        reader = csv.DictReader(fh)
        sample_col, file_col = _resolve_samplesheet_columns(reader.fieldnames)
        for row in reader:
            rows.append(row)

    if not rows:
        print(f"ERROR: {samplesheet_path} has no sample rows.", file=sys.stderr)
        sys.exit(2)

    # WIDENED 2026-09-30: this used to require every row's file to share the
    # EXACT SAME immediate parent directory, which broke the moment a
    # samplesheet covered more than one raw-data folder (e.g. TMA_A, TMA_B,
    # DII_TMA_A, DII_TMA_B -- Afrouz's actual CRC TMA layout). That
    # requirement was stricter than segmentation.py actually needs:
    # segmentation.py's process_images() walks masked_roi_dir with os.walk()
    # (recursive, any depth) and computes each file's own slide_id from its
    # OWN immediate parent folder name (os.path.basename(os.path.dirname(
    # input_file)), segmentation.py line 420) -- it never assumed every file
    # sits directly in masked_roi_dir itself. So masked_roi_dir only needs to
    # be a common ANCESTOR of every row's file, not their shared immediate
    # parent.
    #
    # FIXED 2026-09-30 (same day, real run failure): the first version of
    # this used os.path.commonpath(), which collapsed 's3://bucket/...'
    # down to 's3:/bucket/...' (single slash) -- an invalid S3 URI that
    # made segmentation.py 404 with "masked_roi_dir not found" on the very
    # next real Cirro run. Use the S3-safe _common_dir() helper above
    # instead (still equals the single shared parent in the common case of
    # one folder, so this remains a strict widening, not a behavior change,
    # for existing single-folder samplesheets).
    file_paths = [r[file_col] for r in rows]
    masked_roi_dir = _common_dir(file_paths)
    parents = {_common_dir([p]) for p in file_paths}
    if len(parents) != 1:
        print(
            f"Samplesheet rows span {len(parents)} folders {sorted(parents)} -- "
            f"using their common ancestor '{masked_roi_dir}' as masked_roi_dir "
            "(segmentation.py walks it recursively, so this is fine)."
        )
    group_map = {r[sample_col]: r["group"] for r in rows}

    cfg.setdefault("experiment", {})["image_experiment_group_map"] = group_map
    cfg.setdefault("paths", {})["masked_roi_dir"] = masked_roi_dir

    print(f"Samplesheet : {samplesheet_path} ({len(rows)} samples, columns "
          f"'{sample_col}'/'{file_col}'/'group') -- overrides "
          f"experiment.image_experiment_group_map and paths.masked_roi_dir from --config")
    for r in rows:
        print(f"  {r[sample_col]:14s} group={r['group']:4s} path={r[file_col]}")

    return cfg


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

IMAGE_GIT_SHA_FILE = Path("/app/.image_git_sha")


def _check_image_freshness(expected_commit: str) -> None:
    """
    Fail fast if this container image predates the git commit Nextflow is
    actually running -- added 2026-09-30 after two separate real runs
    silently executed a stale run_pipeline.py baked into an old image
    (COPY'd into the image at build time, unlike preprocess.py, which
    Cirro re-downloads fresh every run -- see preprocess.py's own
    docstring). Each one burned several minutes before failing deep inside
    segmentation with an error that had nothing to do with the real
    problem: the image was just out of date. This check turns that into
    an immediate, unambiguous failure instead.

    expected_commit comes from Nextflow's workflow.commitId (main.nf passes
    it as --expected-commit); IMAGE_GIT_SHA_FILE is written at image build
    time from the GIT_SHA build-arg (see Dockerfile.spacec-base and
    .github/workflows/build-spatia-image.yml). Skipped entirely for a
    local/manual run with no --expected-commit given.
    """
    if not expected_commit:
        return
    image_commit = (
        IMAGE_GIT_SHA_FILE.read_text().strip()
        if IMAGE_GIT_SHA_FILE.exists() else "unknown"
    )
    # Nextflow's commitId and github.sha are both full 40-char SHAs in
    # practice, but compare with startswith() in both directions so a
    # short SHA on either side still matches correctly.
    if image_commit == "unknown" or not (
        expected_commit.startswith(image_commit) or image_commit.startswith(expected_commit)
    ):
        print(
            "ERROR: stale container image -- this image was built from git "
            f"commit '{image_commit}', but Nextflow is running commit "
            f"'{expected_commit}'. Rebuild it: GitHub -> Actions -> "
            "'Build SPATIA Cirro image' -> Run workflow (on this branch), "
            "then re-run.",
            file=sys.stderr,
        )
        sys.exit(3)


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
        "--samplesheet",
        default=None,
        help=(
            "Optional samplesheet CSV (columns: sample or sample_id, file or "
            "image_path, group). When given, "
            "overrides experiment.image_experiment_group_map and "
            "paths.masked_roi_dir from --config with the samplesheet's contents "
            "-- everything else in --config (thresholds, steps, marker panel, "
            "etc.) is used unchanged. Omit this flag to run exactly as before "
            "(sample metadata read from --config as it always has been)."
        ),
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
    parser.add_argument(
        "--expected-commit",
        default=None,
        help=(
            "Git commit this run_pipeline.py is expected to match (Nextflow's "
            "workflow.commitId, passed by main.nf). If it doesn't match the "
            "commit this container image was actually built from, fail "
            "immediately instead of running against stale code. Omit for "
            "local/manual runs."
        ),
    )
    args = parser.parse_args()

    _check_image_freshness(args.expected_commit)

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

    if args.samplesheet:
        cfg = _apply_samplesheet(cfg, Path(args.samplesheet))

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
