#!usr/bin/env nextflow

nextflow.enable.dsl = 2

process CellRanger {
    tag "Running CellRanger count for ${sample_id}"
    publishDir "${params.cellranger_dir}", mode: "copy"

    input:
    tuple val(sample_id), path(fastq_dir)
    path transcriptome
    
    output:
    tuple val(sample_id),  path("${sample_id}_output/outs/filtered_feature_bc_matrix/")

    script:
    """
    cellranger count \\
        --id=${sample_id}_output \\
        --transcriptome=${transcriptome} \\
        --fastqs=${fastq_dir} \\
        --sample=${sample_id} \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        --create-bam=false
    """
}