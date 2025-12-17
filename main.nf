#!usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { QC } from './modules/qc.nf'
include { CellRanger } from './modules/cellranger.nf'
include { Integration } from './modules/integration.nf'
include { Annotation } from './modules/annotation.nf'

// Help Message
def helpMessage() {
    log.info"""
    =====================================================================
                    scRNA-seq Processing Pipeline
    =====================================================================

    Nextflow pipeline for processing single-cell RNA sequencing data,
    including quality control, alignment, integration, and annotation.

    Usage:
        nextflow run main.nf [options]
    
    Required Options:
        --input_csv             Path to samplesheet CSV file (columns: sample_id, fastq_dir)
        --transcriptome         Path to Cell Ranger transcriptome reference
        --integration_script    Path to R script for data integration
        --annotation_script     Path to R script for cluster annotation
    
    Optional Parameters:
        --outdir                Output directory [default: ./results]
        --qc_dir                QC output directory [default: \${outdir}/qc]
        --cellranger_dir        Cell Ranger output directory [default: \${outdir}/cellranger]
        --integration_dir       Integration output directory [default: \${outdir}/integration]
        --annotation_dir        Annotation output directory [default: \${outdir}/annotation]
        
    Skip Options:
        --skip_qc               Skip FastQC step [default: false]
        --skip_cellranger       Skip Cell Ranger and use existing outputs [default: false]
        --fallback_to_cellranger Run Cell Ranger if existing outputs not found [default: true]
    
    Cell Ranger Options:
        --expect_cells          Expected number of cells [default: 3000]
        --chemistry             Chemistry version [default: auto]
    
    Example:
        nextflow run main.nf \\
            --input_csv samples.csv \\
            --transcriptome /path/to/refdata-gex-GRCh38-2020-A \\
            --integration_script scripts/integrate.R \\
            --annotation_script scripts/annotate.R
    """
}

// Validate required parameters
def validateParams() {
    def errors = []
    
    if (!params.input_csv) {
        errors << "Missing required parameter: --input_csv"
    }
    if (!params.transcriptome) {
        errors << "Missing required parameter: --transcriptome"
    }
    if (!params.integration_script) {
        errors << "Missing required parameter: --integration_script"
    }
    if (!params.annotation_script) {
        errors << "Missing required parameter: --annotation_script"
    }
    
    if (errors) {
        log.error "Validation errors:"
        errors.each { error -> log.error "  - ${error}" }
        log.error "\nRun 'nextflow run main.nf --help' for usage information"
        exit 1
    }
}

// Check if Cell Ranger should actually be skipped
def acutallyskipCellRanger() {
    if (!params.skip_cellranger) {
        return false
    }

    def cellrangerDir = file(params.cellranger_dir)
    if (!cellrangerDir.exists()) {
        if (params.fallback_to_cellranger) {
            log.warn "Cell Ranger directory ${params.cellranger_dir} does not exist! Will run cellranger!"
        } else {
            error "No Cell Ranger outputs found in directory ${params.cellranger_dir} and fallback_to_cellranger disabled!"
        }
    }
    return true
}

// Main workflow
workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    // Validate parameters
    validateParams()

    // Read sample sheet and create channel
    // Expected CSV format: sample_id, fastq_dir
    samples_ch = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id
            def fastq_dir = file(row.fastq_dir)

            if (!fastq_dir.exists()) {
                error "FASTQ directory does not exist for sample ${sample_id}: ${fastq_dir}"
            }

            tuple(sample_id, fastq_dir)
        }
    
    // Optional: Run FastQC
    if (!params.skip_qc) {
        log.info "Running FastQC on samples"
        QC(samples_ch)
    } else {
        log.info "Skipping QC as requested"
    }

    // Determine whether to run Cell Ranger or use existing outputs
    def skipCellRanger = acutallyskipCellRanger()
    
    if (skipCellRanger) {
        log.info "Using Existing Cell Ranger outputs from ${params.cellranger_dir}"

        // Create channel from existing Cell Ranger outputs
        cellranger_outputs_ch = channel
            .fromPath("${params.cellranger_dir}/*/outs/filtered_feature_bc_matrix", type: 'dir', checkIfExists: true)
            .map { matrix_dir ->
                // Extract sample_id from directory structure
                // Assumes structure: cellranger_dir/sample_id/outs/filtered_feature_bc_matrix
                def sample_id = matrix_dir.parent.parent.name
                tuple(sample_id, matrix_dir)
            }
    } else {
        log.info "Running Cell ranger on samples..."

        // Run Cell Ranger
        cellranger_outputs_ch = CellRanger(
            samples_ch,
            params.transcriptome
        )
    }

    // COllect all CellRanger outputs for integration
    // The collect() operator will wait for all samples to complete
    all_samples_ch = cellranger_outputs_ch.collect()

    // Run integration on all samples
    log.info "Running integration on all samples"
    integration_output = Integration(
        all_samples_ch,
        params.integration_script
    )

    // Run annotation on integrated data
    Annotation(
        integration_output,
        params.annotation_script
    )

    workflow.onComplete {
        log.info """
        =====================================================================
                            Pipeline Execution Summary
        =====================================================================
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        Work Dir    : ${workflow.workDir}
        Exit status : ${workflow.exitStatus}
        =====================================================================
        """
    }

    workflow.onError {
        log.error """
        =====================================================================
                                Pipeline Error
        =====================================================================
        Error message: ${workflow.errorMessage}
        Error report : ${workflow.errorReport}
        =====================================================================
        """
    }
}