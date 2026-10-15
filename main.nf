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

// Viewer: each core's image, cell outlines and typed cells as a
// spatial-data-studio checkpoint, published with the static viewer under
// ${output_dir}/viewer/. The checkpoint writer (backend/cli.py) and the
// viewer build come from the same spatial-data-studio release, so they agree
// on the checkpoint format; bump both together by changing sds_version.
params.sds_version          = "v1.2.0"
params.sds_source           = "https://github.com/CirroBio/spatial-data-studio/archive/refs/tags/${params.sds_version}.tar.gz"
params.viewer_dist          = "https://github.com/CirroBio/spatial-data-studio/releases/download/${params.sds_version}/viewer-dist.tar.gz"
params.viewer_title         = "SPATIA"
params.viewer_container     = "ghcr.io/astral-sh/uv:python3.11-bookworm"   // squidpy (a cli.py dependency) needs Python 3.11, not the SPATIA image's 3.9
params.viewer_os_packages   = "libgl1 libglib2.0-0 libgomp1"                 // image-tiling libraries missing from bookworm
params.lowres_max_image_mb  = 10                                             // image budget for each core's low-res copy
params.results_dir          = null                                           // viewer_only entry: an existing run's results/ folder

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
    def steps_arg = params.steps_to_run?.trim() ? "--steps " + params.steps_to_run.trim().replaceAll(/\s*,\s*/, ' ') : ''
    def skip_arg   = params.skip_steps?.trim()   ? "--skip "  + params.skip_steps.trim().replaceAll(/\s*,\s*/, ' ')  : ''
    """
    python /app/run_pipeline.py --config ${config_file} --samplesheet ${samplesheet_file} --expected-commit ${workflow.commitId ?: ''} ${steps_arg} ${skip_arg}
    """
}

process build_viewer_checkpoint {
    tag "${label}"

    container params.viewer_container
    cpus 2
    memory { 8.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${params.output_dir}/viewer", mode: 'copy'

    input:
    tuple val(label), path(seg_pickle), path(cell_table), path(triad_pairs)   // triad_pairs: [] for a core with no triads
    path build_script
    path sds_source

    output:
    path "checkpoints/${label}/*.sdata.zarr.zip", emit: checkpoints

    script:
    // cli.py with no recipes only loads the store and saves it as a
    // checkpoint: names it, sets how it opens (cells coloured by cell_type,
    // drawn as outlines once zoomed in) and adds the low-res copy.
    """
    apt-get update
    apt-get install -y --no-install-recommends ${params.viewer_os_packages}
    mkdir sds
    tar -xzf ${sds_source} -C sds --strip-components=1
    uv venv --python 3.11
    uv pip install -r sds/backend/requirements.txt

    .venv/bin/python ${build_script} \\
        --seg-pickle '${seg_pickle}' \\
        --cells '${cell_table}' \\
        ${triad_pairs ? "--triad-pairs '${triad_pairs}'" : ''} \\
        --output core.zarr

    mkdir -p .sds-work
    export SDS_WORK_DIR="\$PWD/.sds-work"
    export SDS_CONTAINER_MEM_MB=${task.memory.toMega()}
    export MPLBACKEND=Agg
    .venv/bin/python sds/backend/cli.py \\
        --parser spatialdata \\
        --input core.zarr \\
        --output 'checkpoints/${label}' \\
        --name results \\
        --lowres-max-image-mb ${params.lowres_max_image_mb} \\
        --lowres-name lowres \\
        --session-name '${label}' \\
        --display-color-by obs:cell_type \\
        --display-render-mode points+shapes
    """
}

process publish_viewer {
    container params.viewer_container

    // Published file by file: the checkpoints land in this same folder, and
    // publishing a directory would replace it rather than merge into it.
    publishDir "${params.output_dir}/viewer", mode: 'copy',
               saveAs: { f -> f.toString().tokenize('/').drop(1).join('/') }

    input:
    path viewer_dist
    val checkpoints

    output:
    path 'site/**'

    script:
    // index.json lists the checkpoints for the viewer's landing page; each
    // path resolves against index.json's own URL. Sorted so the listing does
    // not depend on which core finished first.
    def entries = checkpoints
        .collect { f -> [label: f.parent.name, file: f.name] }
        .sort { a, b -> a.label <=> b.label ?: a.file.startsWith('lowres-') <=> b.file.startsWith('lowres-') }
        .collect { e -> [
            path : "checkpoints/${e.label}/${e.file}".toString(),
            label: e.file.startsWith('lowres-') ? "${e.label} (low-res image)".toString() : e.label,
        ] }
    def index_json = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson([title: params.viewer_title, checkpoints: entries]))
    """
    mkdir site
    tar -xzf ${viewer_dist} -C site
    cat > site/index.json <<'INDEX_JSON'
${index_json}
INDEX_JSON
    """
}

def stripSuffix(String name, List suffixes) {
    def suffix = suffixes.find { s -> name.endsWith(s) }
    return name.substring(0, name.length() - suffix.length())
}

/*
 * One viewer checkpoint per core, from the files run_spatia_pipeline wrote.
 * A core is keyed by its seg pickle (<core>_seg_output.pickle); its cell table
 * is named <group>_<core>_..., and that <group>_<core> stem is the label the
 * viewer shows. Its <group>_<core>_triad_pairs.csv, written only when the core
 * has triads, becomes the triangles drawn over it.
 */
workflow SPATIA_VIEWER {
    take:
    results     // every file under results/

    main:
    def seg_pickle_suffix = "_seg_output.pickle"
    def triad_pairs_suffix = "_triad_pairs.csv"
    // Most informative first. triads.py writes no flags file for a core with
    // zero triads, so the matched-cells table stands in for it.
    def cell_table_suffixes = ["_cells_with_triad_flags.csv", "_matched_with_boundaries.csv"]

    def files = results.flatten()
    def pickles = files
        .filter { f -> f.name.endsWith(seg_pickle_suffix) }
        .map { f -> [stripSuffix(f.name, [seg_pickle_suffix]), f] }
    def cell_tables = files
        .filter { f -> cell_table_suffixes.any { s -> f.name.endsWith(s) } }
    def triad_pairs = files
        .filter { f -> f.name.endsWith(triad_pairs_suffix) }
        .map { f -> [stripSuffix(f.name, [triad_pairs_suffix]), f] }

    def cores = pickles
        .combine(cell_tables)
        .filter { core, _pickle, table -> stripSuffix(table.name, cell_table_suffixes).endsWith("_${core}") }
        .groupTuple(by: [0, 1])
        .map { _core, pickle, tables ->
            def table = tables.min { t -> cell_table_suffixes.findIndexOf { s -> t.name.endsWith(s) } }
            [stripSuffix(table.name, cell_table_suffixes), pickle, table]
        }
        .join(triad_pairs, remainder: true)
        // remainder pads with a single null, so an unmatched entry is
        // [label, null, pairs] (pairs whose core has no pickle: dropped) or
        // [label, pickle, table, null] (a core with no triads).
        .filter { entry -> entry[1] instanceof Path }
        .map { entry -> [entry[0], entry[1], entry[2], entry[3] ?: []] }

    build_viewer_checkpoint(
        cores,
        file("${projectDir}/bin/build_spatialdata.py"),
        file(params.sds_source),
    )
    publish_viewer(
        file(params.viewer_dist),
        build_viewer_checkpoint.out.checkpoints.flatten().toList(),
    )
}

workflow {
    config_ch      = Channel.fromPath(params.config)
    samplesheet_ch = Channel.fromPath(params.samplesheet)
    run_spatia_pipeline(config_ch, samplesheet_ch)
    SPATIA_VIEWER(run_spatia_pipeline.out.results)
}

// Builds only the viewer, from an existing run's results/ folder:
//   nextflow run main.nf -entry viewer_only --results_dir path/to/results
workflow viewer_only {
    if( !params.results_dir ) error "viewer_only needs --results_dir (an existing run's results/ folder)"
    SPATIA_VIEWER(channel.fromPath("${params.results_dir}/**"))
}
