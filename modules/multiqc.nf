#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MultiQC {
    tag "Aggregating QC reports"
    publishDir "${params.multiqc_dir}", mode: 'copy'

    input:
    path fastqc_files, stageAs: 'fastqc_data/*'
    path cellranger_dirs, stageAs: 'cellranger_data/*'

    output:
    path('*.html')
    path('multiqc_report_data'), type: 'dir'

    script:
    """
    # Create organized directory structure
    mkdir -p input_data/fastqc
    mkdir -p input_data/cellranger

    # Copy FastQC files
    if [ -d "fastqc_data" ]; then
        cp fastqc_data/* input_data/fastqc/ 2>/dev/null || true
    fi

    # Copy Cell Ranger sample directories - each has unique name (e.g., OSA_001_output)
    for dir in cellranger_data/*; do
        if [ -d "\$dir" ]; then
            # Copy the entire sample directory which contains outs/web_summary.html
            cp -r "\$dir" input_data/cellranger/
        fi
    done

    # Run MultiQC
    multiqc input_data \\
        --filename multiqc_report.html \\
        --title "scRNA-seq QC Report" \\
        --force \\
        --config <(echo "
module_order:
  - fastqc
  - cellranger
")
    """
}