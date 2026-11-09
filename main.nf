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

// steps_to_run / skip_steps -- added 2026-10-12 at Afrouz's request, so Cirro
// users can skip expensive upstream steps (segmentation, preprocessing) when
// iterating on cell_typing/triads/functional/survival parameters, instead of
// re-running the whole 30-45min pipeline every time. run_pipeline.py has
// supported --steps/--skip since before this date -- this just forwards
// Cirro's new "Steps to run" / "Steps to skip" form fields (comma-separated,
// .cirro/process-form.json) to those existing flags. Both default to empty,
// which reproduces the exact command this process ran before this change --
// no behavior change for anyone who leaves the new fields blank.
// NOTE (unverified, flag for Afrouz to confirm on first real test): skipping
// an upstream step does not fetch that step's input for you -- whoever runs
// with a restricted --steps/--skip is responsible for cfg's input path (e.g.
// cell_typing.input_file) already pointing at a completed prior run's output.
params.steps_to_run = params.steps_to_run ?: ''
params.skip_steps   = params.skip_steps   ?: ''
// previous_segmentation -- added 2026-11-05. S3 folder holding a previous
// run's results/segmentation_results (*_mesmer_result.csv + *_seg_output.pickle
// per slide folder). Set by .cirro/preprocess.py from the "Reuse segmentation
// from run" form field. When set, it is staged into this task and copied to
// results/segmentation_results, and preprocess.py has already added
// segmentation to skip_steps. Empty = segment from scratch (old behaviour).
params.previous_segmentation = params.previous_segmentation ?: ''
// previous_preprocessing -- added 2026-11-09. S3 folder with a previous run's
// results/combined_processed_data (per-tissue h5ads). Set by preprocess.py
// when "Reuse through step" = preprocessing; copied to
// results/combined_processed_data so cell typing can start directly.
params.previous_preprocessing = params.previous_preprocessing ?: ''

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
    path prev_seg, stageAs: 'prev_segmentation_input'
    path prev_pre, stageAs: 'prev_preprocessing_input'

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
    def steps_arg = params.steps_to_run?.trim() ? "--steps " + params.steps_to_run.trim().replaceAll(/\s*,\s*/, ' ') : ''
    def skip_arg   = params.skip_steps?.trim()   ? "--skip "  + params.skip_steps.trim().replaceAll(/\s*,\s*/, ' ')  : ''
    """
    if [ -d prev_segmentation_input ]; then
        echo "Reusing previous segmentation from ${params.previous_segmentation}"
        mkdir -p results/segmentation_results
        cp -rL prev_segmentation_input/. results/segmentation_results/
        echo "  copied \$(find results/segmentation_results -name '*_mesmer_result.csv' | wc -l) mesmer result files"
    fi
    if [ -d prev_preprocessing_input ]; then
        echo "Reusing previous preprocessing from ${params.previous_preprocessing}"
        mkdir -p results/combined_processed_data
        cp -rL prev_preprocessing_input/. results/combined_processed_data/
        echo "  copied \$(find results/combined_processed_data -name '*_combined_all_experiment_groups.h5ad' | wc -l) tissue h5ad files"
    fi
    python /app/run_pipeline.py --config ${config_file} --samplesheet ${samplesheet_file} --expected-commit ${workflow.commitId ?: ''} ${steps_arg} ${skip_arg}
    """
}

workflow {
    config_ch      = Channel.fromPath(params.config)
    samplesheet_ch = Channel.fromPath(params.samplesheet)
    // No previous results -> pass a tiny placeholder FILE; the process only
    // copies when the staged input is a DIRECTORY. 2026-11-09: placeholders
    // are now written to the work dir at run time -- the repo's .gitignore
    // ignores everything outside spatia/ etc., so assets/NO_FILE was never
    // pushed and Nextflow failed with "Can't stage file .../assets/NO_FILE".
    def noSeg = file("${workflow.workDir}/spatia_placeholders/NO_PREVIOUS_SEGMENTATION")
    def noPre = file("${workflow.workDir}/spatia_placeholders/NO_PREVIOUS_PREPROCESSING")
    if( !noSeg.exists() ) { noSeg.parent.mkdirs(); noSeg.text = 'none\n' }
    if( !noPre.exists() ) { noPre.parent.mkdirs(); noPre.text = 'none\n' }
    prev_seg_ch    = params.previous_segmentation
                       ? Channel.fromPath(params.previous_segmentation, type: 'dir', checkIfExists: true)
                       : Channel.value(noSeg)
    prev_pre_ch    = params.previous_preprocessing
                       ? Channel.fromPath(params.previous_preprocessing, type: 'dir', checkIfExists: true)
                       : Channel.value(noPre)
    run_spatia_pipeline(config_ch, samplesheet_ch, prev_seg_ch, prev_pre_ch)
}
