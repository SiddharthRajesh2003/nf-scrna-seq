#!usr/bin/env nextflow

nextflow.enable.dsl = 2

include { QC } from './modules/qc.nf'
include { CellRanger } from './modules/cellranger.nf'
include { Integration } from './modules/integration.nf'
include { Annotation } from './modules/annotation.nf'

def helpMessage() {
    log.info"""
    =====================================================================
                    scRNA-seq Processing Pipeline
    =====================================================================

    Nextflow pipeline for processing single-cell RNA sequencing data,
    including quality control, alignment, integration, and annotation.

    Usage:
        nextflow run main.nf [options]
    
    Options:
        --samples               Path to samplesheet CSV file [required]
        --transcriptome         Path to cellranger transcriptome reference [required]
        --integration_script    Path to R script for data integration [required]
        --annotation_script     Path to R script for cluster annotation [required]
    """
}

def acutallyskipCellRanger() {
    if (params.skip_cellranger) {
        log.info "Skipping CellRanger step as per user request."
        return true
    }

}