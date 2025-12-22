#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MultiQC {
    tag "Aggregating QC reports"
    publishDir "${params.multiqc_dir}", mode: 'copy'

    input:
    path('*')

    output:
    path('*.html')
    path('multiqc_report_data'), type: 'dir'

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --title "scRNA-seq QC Report"
    """

}