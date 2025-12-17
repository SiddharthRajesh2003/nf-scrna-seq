#!usr/env/bin nextflow

nextflow.enable.dsl = 2

process QC {
    tag "Performing Quality Control on ${sample_id}"
    publishDir "${params.fastqqc_dir}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_dir)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${fastq_dir}/*.fastq.gz --threads ${task.cpus}
    """
}