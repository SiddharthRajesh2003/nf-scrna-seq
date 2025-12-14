#!usr/bin/env nextflow

nextflow.enable.dsl = 2

process Annotation {
    tag "Annotating clusters in ${input_rds}"
    publishDir "${params.annotation_dir}", mode: 'copy'

    input:
    path input_rds

    output:
    path "*.rds"
    path "annotation_plots/"
    path "annotation_reports.html", optional: true
    path "cell_type_markers.csv", optional: true
    
    script:
    """
    Rscript ${params.annotation_script} \\
        --seurat_object ${input_rds} \\
        --output annotated_seurat.rds \\
        --reference ${params.reference_dataset} \\
        --method ${params.annotation_methods}
    """
}