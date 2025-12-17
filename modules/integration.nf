#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process Integration {
    tag "Integrating ${sample_ids.size()} samples using ${params.integration_method}"
    publishDir "${params.integration_dir}", mode: 'copy'

    input:
    tuple val(sample_ids), path(matrix_dirs)
    path integration_script

    output:
    path "integrated_seurat.rds", emit: seurat_object
    path "qc_plots/*", emit: qc_plots
    path "integration_summary.txt", emit: summary

    script:
    // Join paths with commas for R script
    def input_paths = matrix_dirs.join(',')
    
    """
    Rscript ${integration_script} \\
        --input_dirs ${input_paths} \\
        --output integrated_seurat.rds \\
        --min_cells ${params.min_cells ?: 3} \\
        --min_features ${params.min_features ?: 200} \\
        --max_features ${params.max_features ?: 7000} \\
        --min_umi ${params.min_umi ?: 1000} \\
        --max_umi ${params.max_umi ?: 50000} \\
        --max_mt_percent ${params.max_mt_percent ?: 20} \\
        --dims ${params.dims ?: 30} \\
        --resolution ${params.resolution ?: 0.3} \\
        --integration_method ${params.integration_method ?: 'cca'}
    """
}