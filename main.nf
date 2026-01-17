#!usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { QC } from './modules/qc.nf'
include { MultiQC } from './modules/multiqc.nf'
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
        --samplesheet             Path to samplesheet CSV file (columns: sample_id, fastq_dir)
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
        --skip_fastqc               Skip FastQC step [default: false]
        --skip_cellranger       Skip Cell Ranger and use existing outputs [default: false]
        --fallback_to_cellranger Run Cell Ranger if existing outputs not found [default: true]
        --skip_downstream_analysis   Skip integration and annotation steps [default: false]
    
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
    
    if (!params.samplesheet) {
        errors << "Missing required parameter: --samplesheet"
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
def actuallyskipCellRanger() {
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
    qc_samples_ch = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sample_id = row.sample_id
            def fastq_files = file("${row.fastq_dir}/*.fastq.gz")
            tuple(sample_id, fastq_files)
        }
    
    cellranger_samples_ch = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header:true)
        .map { row ->
            def sample_id = row.sample_id
            def fastq_dir = file(row.fastq_dir)
            tuple(sample_id, fastq_dir)
        }

    // Optional: Run FastQC
    if (!params.skip_fastqc) {
        log.info "Running FastQC on samples"
        QC(qc_samples_ch)
        fastqc_reports = QC.out.zip.collect()
    } else {
        log.info "Skipping FastQC as requested"
        fastqc_reports = channel.empty()
    }

    // Determine whether to run Cell Ranger or use existing outputs
    def skipCellRanger = actuallyskipCellRanger()

    if (skipCellRanger) {
        log.info "Using existing Cell Ranger outputs from ${params.cellranger_dir}"

        // Create channel from existing Cell Ranger outputs
        cellranger_outputs_ch = channel
            .fromPath("${params.cellranger_dir}/*/outs/filtered_feature_bc_matrix", type: 'dir', checkIfExists: true)
            .map { matrix_dir ->
                // Extract sample_id from directory structure
                // Assumes structure: cellranger_dir/sample_id_output/outs/filtered_feature_bc_matrix
                def sample_id = matrix_dir.parent.parent.name.replaceAll("_output", "")
                tuple(sample_id, matrix_dir)
            }

        // Create channel for existing sample output directories (preserves unique sample names)
        web_summary_ch = channel
            .fromPath("${params.cellranger_dir}/*_output", type: 'dir', checkIfExists: true)
            .toList()

    } else {
        log.info "Running Cell Ranger on samples..."

        // Run Cell Ranger
        CellRanger(
            cellranger_samples_ch,
            params.transcriptome
        )

        // Extract filtered_feature_bc_matrix for downstream analysis
        cellranger_outputs_ch = CellRanger.out.cellranger_out.map {
            sample_id, outs_dir ->
                tuple(sample_id, file("${outs_dir}/filtered_feature_bc_matrix"))
        }

        // Collect sample output directories for MultiQC (preserves unique sample names)
        web_summary_ch = CellRanger.out.cellranger_out
            .map { _sample_id, outs_dir -> outs_dir.parent }
            .toList()
    }

    // Run MultiQC if FastQC was run OR if we have Cell Ranger summaries
    if (!params.skip_fastqc || !skipCellRanger) {
        MultiQC(
            fastqc_reports.ifEmpty([]),
            web_summary_ch.ifEmpty([])
        )
    } else {
        log.info "Skipping MultiQC - no QC reports to aggregate"
    }

    // Run downstream analysis (integration and annotation)
    if (!params.skip_downstream_analysis) {
        log.info "Running downstream analysis (integration and annotation)"
        DownstreamAnalysis(
            cellranger_outputs_ch
        )
    } else {
        log.info "Skipping downstream analysis as requested"
    }
}

workflow DownstreamAnalysis {
    take:
        cellranger_outputs_ch

    main:
    // Collect all samples as list of tuples
    all_samples = cellranger_outputs_ch.toList()

    // Run integration on all samples
    log.info "Running integration on all samples"
    Integration(
        all_samples,
        params.integration_script
    )

    // Run annotation on integrated data
    Annotation(
        Integration.out.seurat_object,
        params.annotation_script
    )
}