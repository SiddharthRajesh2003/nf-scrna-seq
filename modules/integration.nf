#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process Integration {
    tag "Integrating samples using ${params.integration_method}"
    publishDir "${params.integration_dir}", mode: 'copy'

    // Add environment variables
    beforeScript """
    export LANG=en_US.UTF-8
    export LC_ALL=en_US.UTF-8
    export R_LIBS=/N/project/Krolab/Siddharth/Pipelines/scrna-seq/conda_envs/env/renv/lib/R/library
    export RENV_CONFIG_SANDBOX_ENABLED=FALSE
    export RENV_PATHS_ROOT=/N/project/Krolab/Siddharth/Pipelines/scrna-seq/renv_root
    export RENV_PATHS_CACHE=/N/project/Krolab/Siddharth/Pipelines/scrna-seq/renv_cache
    export RENV_PATHS_SANDBOX=/N/project/Krolab/Siddharth/Pipelines/scrna-seq/renv_cache/sandbox
    """

    input:
    val sample_data  // List of [sample_id, matrix_path] tuples
    path integration_script

    output:
    path "integrated_seurat.rds", emit: seurat_object
    path "qc_plots/*", emit: qc_plots
    path "integration_summary.txt", emit: summary

    script:
    // Create comma-separated list of paths
    def input_paths = sample_data.collect { tuple -> 
        "${tuple[1]}"
    }.join(',')
    
    """
    # echo "R version:"
    R --version | head -n1
    
    echo "R library paths:"
    R --slave -e '.libPaths()'
    
    echo "Checking for Seurat:"
    R --slave -e 'if(file.exists("/N/project/Krolab/Siddharth/Pipelines/scrna-seq/conda_envs/env/renv/lib/R/library/Seurat")) cat("Seurat found in conda\\n")'
    
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