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
 * NOT YET RUN — no Nextflow/Docker runtime available in the sandbox this
 * was written in. Needs a real test (small CRC config, e.g.
 * experiments/crc_tma.yaml) before trusting the channel/publish wiring.
 */

nextflow.enable.dsl = 2

params.config       = "experiments/crc_tma.yaml"   // path to the SPATIA experiment YAML
params.output_dir   = "results"                     // Nextflow-local publish dir
params.container    = "spatia-pipeline:latest"      // image built from ./Dockerfile

process run_spatia_pipeline {
    tag "${params.config}"

    container params.container

    publishDir params.output_dir, mode: 'copy'

    input:
    path config_file

    output:
    path "results/**", emit: results
    path "logs/**",     emit: logs

    script:
    """
    python /app/run_pipeline.py --config ${config_file}
    """
}

workflow {
    config_ch = Channel.fromPath(params.config)
    run_spatia_pipeline(config_ch)
}
