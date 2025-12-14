#!usr/env/bin nextflow

nextflow.enable.dsl = 2

process Integration {
    tag "Integrating datasets using ${integration_method}"
    publishDir "${params.integration_dir}", mode: 'copy'

    input:
    tuple val(sample_names), path(input_dirs)
    val integration_method

    output:
    path "*.rds"
    path "qc_plots/"

    script:
    """
    Rscript ${params.integration_script} \\
        --input_dirs ${input_dirs.join(',')} \\
        --sample_names ${sample_names.join(',')} \\
        --method ${integration_method} \\
        --output_rds integrated_seurat.rds
    """
}