#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process QC {
    tag "Performing Quality Control on ${sample_id}"
    publishDir "${params.fastqc_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_files)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:

    """
    fastqc ${fastq_files} \\
        --threads ${task.cpus} \\
        --outdir .
    """
}