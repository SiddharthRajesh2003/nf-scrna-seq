#!usr/bin/env nextflow

nextflow.enable.dsl = 2

process Annotation {
    tag "Annotating clusters using ${params.annotation_method}"
    publishDir "${params.annotation_dir}", mode: 'copy'

    input:
    path seurat_rds
    path annotation_script

    output:
    path "annotated_seurat.rds", emit: annotated_object
    path "annotation_plots/*", emit: plots
    path "cell_type_markers.csv", emit: markers
    path "annotation_summary.txt", emit: summary
    
    script:
    """
    Rscript ${annotation_script} \\
        --seurat_object ${seurat_rds} \\
        --output annotated_seurat.rds \\
        --method ${params.annotation_method} \\
        --reference ${params.reference_dataset} \\
        --species ${params.species} \\
        --azimuth_reference ${params.azimuth_reference}
    """
}