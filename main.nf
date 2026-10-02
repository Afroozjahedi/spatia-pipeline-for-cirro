#!/usr/bin/env nextflow
/*
 * main.nf — SPATIA pipeline, single-process Nextflow wrapper (Q10, 2026-07-21)
 *
 * Packaging decision (Q10, decided by Afrouz 2026-07-21): ONE black-box
 * Nextflow process calls run_pipeline.py end-to-end inside the Docker image
 * built from ./Dockerfile — not five staged per-step processes. Lower
 * engineering cost, at the tradeoff of no per-step retry/resourcing/
 * Cirro-visualized intermediate outputs (that tradeoff was flagged to
 * Afrouz in SPATIA_PIPELINE_LOG.md Day 2 and she chose this option anyway).
 *
 * Scope: CRC TMA only (this project's confirmed dataset scope, Q4/Q6).
 * A second, separate-purpose experiment config was added 2026-07-29 (Q24):
 * experiments/gbm_c3_demo.yaml (GBM CyCIF, out of this project's CRC scope
 * per Q4/Q6, run at Afrouz's explicit request). No changes to this
 * workflow were needed for that — it's driven entirely by --config, e.g.:
 *
 *   nextflow run main.nf --config experiments/gbm_c3_demo.yaml \
 *       --container spatia-pipeline:latest
 *
 * (build spatia-pipeline:latest from ../Dockerfile, or from
 * ../Dockerfile.spacec-base if using the real spacec container as a base
 * — see that file's header for why both exist and what's unverified in
 * each. Pass --container spatia-pipeline-spacec-base:latest or similar if
 * you build the second one instead.)
 *
 * SAMPLESHEET — added 2026-09-24, in response to Cirro admin (Dima)
 * feedback: sample metadata (which image belongs to which sample/group)
 * now comes from a separate samplesheet CSV (sample_id,image_path,group)
 * instead of being hand-written into --config. This is STILL Q10's single
 * batched process — one call handles every sample in the samplesheet, the
 * same way run_pipeline.py already walked --config's masked_roi_dir
 * before. See run_pipeline.py's _apply_samplesheet() for exactly how the
 * override works, and its docstring for the one constraint it enforces
 * (all samplesheet rows must share one parent directory, since
 * segmentation.py still walks a single masked_roi_dir).
 * --config's own experiment.image_experiment_group_map / paths.masked_roi_dir
 * are IGNORED whenever --samplesheet is passed (samplesheet wins); if you
 * omit --samplesheet, behavior is unchanged from before this date.
 *
 * NOT YET RUN — no Nextflow/Docker runtime available in the sandbox this
 * was written in. Needs a real test (small CRC config, e.g.
 * experiments/crc_tma.yaml) before trusting the channel/publish wiring.
 */

nextflow.enable.dsl = 2

params.outdir       = params.outdir ?: null                         // declared explicitly to silence Nextflow's
                                                                     // "Access to undefined parameter `outdir`" warning
                                                                     // (harmless -- Cirro always supplies it via
                                                                     // .cirro/process-input.json's "outdir" mapping --
                                                                     // but noisy in every Cirro run's log otherwise)
params.config       = "experiments/crc_tma_full_pipeline.params.yaml"  // pipeline parameters only (no sample metadata)
params.samplesheet  = "experiments/samplesheet.csv"                    // sample_id,image_path,group -- optional, see run_pipeline.py --help
params.output_dir   = params.outdir ?: "results"     // FIXED (2026-10-02): Cirro injects the dataset's real
                                                      // S3 output path as params.outdir via .cirro/process-input.json's
                                                      // "outdir": "$.dataset.dataPath" mapping (confirmed against a working
                                                      // sibling Cirro pipeline at this org, btc-oncoanalyser, which uses the
                                                      // same pattern). Without this, params.output_dir stayed a LOCAL relative
                                                      // "results" path -- publishDir copied into the ephemeral head-job
                                                      // container's scratch space, not to Cirro's dataset, so every run that
                                                      // completed (even a clean exit 0, e.g. CRC_run 28) published nothing.
                                                      // "results" remains the default for manual/local runs where params.outdir
                                                      // is never set.
params.container    = params.container ?: "spatia-pipeline:latest"  // default for manual/local runs -- Cirro overrides this via process-compute.config's params.container (see 2026-09-30 fix: this used to be an unconditional assignment that silently stomped any container Cirro supplied)

process run_spatia_pipeline {
    tag "${params.config}"

    container params.container

    publishDir params.output_dir, mode: 'copy'

    // FIXED (2026-10-01): this process's own `results/**` glob is relative to
    // the task's work dir, but every write-target path in the yaml configs
    // (paths.output_dir, paths.segmentation_results_dir,
    // segmentation.normalized_dir, analysis.survival.image_patient_map,
    // paths.input_dir) used to be an ABSOLUTE /rsrch6/... path -- nothing
    // the pipeline wrote ever landed under this task's work dir, so this
    // glob matched zero files on every run (success or failure) and nothing
    // was ever published to Cirro. Those configs now use paths nested under
    // a relative "results/" tree so this glob actually finds them. Dropped
    // the separate logs/** emit channel below (run_pipeline.py's own logs
    // already nest under {output_dir}/logs/, i.e. results/logs/ -- covered
    // by results/** -- and a second top-level-only glob with nothing to
    // match it risked its own "missing output file(s)" failure).
    //
    // REVERTED same day: this process briefly also declared
    // `validExitStatus 0, 4`, meant to let a validation-halted run (exit 4
    // from run_pipeline.py) still publish its partial output. That directive
    // does not exist in real Nextflow -- confirmed the hard way, it broke
    // this file's parsing outright on a real Cirro run ("Unknown process
    // directive: validExitStatus"). Reverted; run_pipeline.py exits 1 for
    // any failure again. A clean exit-0 run publishes fine with the fix
    // above -- publishing a validation-halted run's partial output is a
    // real, still-open problem, not solved here.

    input:
    path config_file
    path samplesheet_file

    output:
    path "results/**", emit: results

    script:
    // Added 2026-09-30: pass the git commit Nextflow actually checked out
    // (workflow.commitId) into run_pipeline.py's --expected-commit, so it
    // can refuse to run if the container image is older than this commit
    // -- see run_pipeline.py's _check_image_freshness for why (two real
    // runs silently executed stale baked-in code otherwise). Empty string
    // for a run not launched from a git repo (workflow.commitId is null),
    // which skips the check entirely.
    """
    python /app/run_pipeline.py --config ${config_file} --samplesheet ${samplesheet_file} --expected-commit ${workflow.commitId ?: ''}
    """
}

workflow {
    config_ch      = Channel.fromPath(params.config)
    samplesheet_ch = Channel.fromPath(params.samplesheet)
    run_spatia_pipeline(config_ch, samplesheet_ch)
}
