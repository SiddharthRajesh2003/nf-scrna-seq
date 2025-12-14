#!usr/env/bin nextflow

nextflow.enable.dsl = 2

process QC {
    tag "Performing Quality Control on ${sample_id}"
    publishDir "${params.fastqqc_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(r2)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${r2} --threads ${task.cpus}
    """
}