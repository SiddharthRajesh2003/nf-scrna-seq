# nf-scrna-seq

A comprehensive Nextflow pipeline for automated single-cell RNA sequencing (scRNA-seq) data analysis, from raw FASTQ files to annotated cell types.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A525.10.2-brightgreen.svg)](https://www.nextflow.io/)
[![R](https://img.shields.io/badge/R-4.4.3-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.4.0-orange.svg)](https://satijalab.org/seurat/)

## Overview

This pipeline automates the complete analysis workflow for single-cell RNA-seq data, including:

- Quality control (FastQC, MultiQC)
- Alignment and quantification (Cell Ranger)
- Data integration and clustering (Seurat)
- Automated cell type annotation (Azimuth/SingleR)
- Interactive visualization (Shiny app)

The pipeline is designed for HPC environments with SLURM scheduling and supports flexible integration methods and annotation strategies.

## Table of Contents

- [Features](#features)
- [Pipeline Workflow](#pipeline-workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Specifications](#input-specifications)
- [Configuration](#configuration)
- [Usage](#usage)
- [Output Files](#output-files)
- [Advanced Options](#advanced-options)
- [Interactive Visualization](#interactive-visualization)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## Features

- **Automated workflow**: End-to-end analysis from FASTQ to annotated cell types
- **Modular design**: Skip steps as needed (QC, Cell Ranger)
- **Flexible integration**: Support for CCA, RPCA, Harmony, and FastMNN
- **Multiple annotation methods**: Azimuth and SingleR cell type annotation
- **Quality control**: Comprehensive filtering and QC metrics
- **Scalable**: Optimized for HPC clusters with SLURM
- **Reproducible**: Conda environment specification and Nextflow resume capability
- **Interactive visualization**: Built-in Shiny app for exploring results

## Pipeline Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│ 1. QUALITY CONTROL (Optional)                                   │
│    ├── FastQC: Raw read quality metrics                         │
│    └── MultiQC: Aggregate QC reports                            │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ 2. ALIGNMENT & QUANTIFICATION (Optional)                        │
│    └── Cell Ranger: Alignment, UMI counting, cell calling       │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ 3. INTEGRATION                                                  │
│    ├── Load Cell Ranger outputs                                 │
│    ├── QC filtering (features, UMI, mitochondrial %)            │
│    ├── Normalization and scaling                                │
│    ├── Integration (CCA/RPCA/Harmony/FastMNN)                   │
│    ├── Dimensionality reduction (PCA, UMAP, t-SNE)              │
│    └── Clustering                                               │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ 4. ANNOTATION                                                   │
│    ├── Cell type annotation (Azimuth/SingleR)                   │
│    ├── Find cluster markers                                     │
│    └── Generate visualization plots                             │
└─────────────────────────────────────────────────────────────────┘
```

## Requirements

### System Requirements

- Linux/Unix operating system
- HPC cluster with SLURM scheduler (recommended)
- Minimum 120 GB RAM for Cell Ranger and integration steps
- High-performance storage for large FASTQ and BAM files

### Software Dependencies

The pipeline uses a comprehensive conda environment with 713 dependencies. Key components include:

- Nextflow ≥25.10.2
- R ≥4.4.3
- Python ≥3.14.2
- MultiQC 1.33
- FastQC (latest)
- Cell Ranger (loaded via module system)

**R packages:**
- Seurat 5.4.0
- Azimuth 0.5.0
- SingleR 2.8.0
- CellDex 1.16.0
- tidyverse, ggplot2, patchwork, optparse

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/SiddharthRajesh2003/nf-scrna-seq
cd nf-scrna-seq
```

### 2. Set up the conda environment

The pipeline includes a comprehensive environment specification in [env.yaml](env.yaml):

```bash
# Create conda environment
conda create -p $base/conda_envs/env/renv \
    jq \
    bioconda::bioconductor-genomeinfodbdata \
    r-base \
    r-essentials \
    r-igraph \
    r-matrix \
    r-rcpp \
    r-rcppeigen \
    r-rcppprogress \
    r-seurat \
    bioconda::bioconductor-singler \
    r-azimuth \
    bioconda::bioconductor-celldex \
    r-tidyverse \
    r-optparse \
    r-ggplot2 \
    r-patchwork \
    r-devtools \
    nextflow \
    multiqc \
    nf-core

# Activate environment
conda activate /path/to/conda/env

# Fix Azimuth GSL library dependency
cd /path/to/conda/env/lib
# Azimuth needs this gsl library for its dependency but newest version is usually downloaded so we apply relation as gsl 2.7 is backwards compatible
ln -s libgsl.so.27 libgsl.so.25 
```

**Note:** The symbolic link for `libgsl.so.25` is required because Azimuth depends on an older GSL version, but the newest version (2.7) is backwards compatible.

### 3. Download Cell Ranger reference genome

```bash
# Download human reference (GRCh38-2024-A)
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
```

### 4. Configure paths

Edit [nextflow.config](nextflow.config) to set your paths:

```groovy
params {
    base = "/path/to/your/project"
    transcriptome = "${params.base}/refdata-gex-GRCh38-2024-A"
    outdir = "${params.base}/results"
}
```

## Quick Start

### Basic usage

```bash
# Using SLURM submission script
sbatch run_pipeline.sh

# Or run directly with Nextflow
nextflow run main.nf \
    --samplesheet metadata/samples.csv \
    -resume \
    -profile slurm
```

### Using existing Cell Ranger outputs

If you already have Cell Ranger outputs, skip the alignment step:

```bash
nextflow run main.nf \
    --samplesheet metadata/samples.csv \
    --skip_cellranger \
    -resume \
    -profile slurm
```

## Input Specifications

### Sample Sheet

Download an SRA Metadata file and use generate_samplesheet.ipynb to generate a script for downloading the FASTQ files, renaming them accordingly, and creating the samplesheet CSV file.

```csv
sample_id,fastq_dir
OSA_001,/path/to/fastqs/OSA_001
OSA_002,/path/to/fastqs/OSA_002
NoOSA_001,/path/to/fastqs/NoOSA_001
```

- **sample_id**: Unique identifier for each sample
- **fastq_dir**: Path to directory containing FASTQ files for that sample

See [metadata/samples.csv](metadata/samples.csv) for an example.

### FASTQ Files

FASTQ files should follow Cell Ranger naming conventions:
- `{sample_id}_S1_L001_R1_001.fastq.gz` (Read 1)
- `{sample_id}_S1_L001_R2_001.fastq.gz` (Read 2)

## Configuration

### Main Configuration Parameters

Edit [nextflow.config](nextflow.config) to customize:

```groovy
params {
    // Input/Output
    samplesheet = "metadata/samples.csv"
    outdir = "./results"
    transcriptome = "/path/to/refdata-gex-GRCh38-2024-A"

    // Workflow control
    skip_fastqc = false
    skip_cellranger = false
    fallback_to_cellranger = true

    // Integration parameters
    integration_method = 'CCA'  // Options: CCA, RPCA, Harmony, FastMNN
    min_cells = 3
    min_features = 200
    max_features = 7000
    min_umi = 1000
    max_umi = 50000
    max_mito = 20
    dims = 30
    resolution = 0.3

    // Annotation parameters
    annotation_method = 'Azimuth'  // Options: Azimuth, SingleR
    reference = 'pbmc'  // For Azimuth
    species = 'human'
    singler_ref = 'HumanPrimaryCellAtlasData'  // For SingleR
}
```

### Resource Allocation

Adjust resources in [nextflow.config](nextflow.config):

```groovy
process {
    withName: FASTQC {
        cpus = 4
        memory = '20 GB'
        time = '6h'
    }
    withName: CELLRANGER_COUNT {
        cpus = 16
        memory = '120 GB'
        time = '72h'
    }
    withName: INTEGRATE_SAMPLES {
        cpus = 12
        memory = '120 GB'
        time = '30h'
    }
}
```

## Usage

### Required Parameters

```bash
nextflow run main.nf \
    --samplesheet <path/to/samples.csv> \
    [options]
```

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--outdir` | Output directory | `./results` |
| `--skip_fastqc` | Skip FastQC step | `false` |
| `--skip_cellranger` | Skip Cell Ranger (use existing outputs) | `false` |
| `--fallback_to_cellranger` | Run Cell Ranger if outputs missing | `true` |
| `--integration_method` | Integration method (CCA/RPCA/Harmony/FastMNN) | `CCA` |
| `--annotation_method` | Annotation method (Azimuth/SingleR) | `Azimuth` |
| `--resolution` | Clustering resolution | `0.3` |

### Example Commands

#### Full pipeline with default settings
```bash
nextflow run main.nf --samplesheet metadata/samples.csv -resume
```

#### Skip QC and use existing Cell Ranger outputs
```bash
nextflow run main.nf \
    --samplesheet metadata/samples.csv \
    --skip_fastqc \
    --skip_cellranger \
    -resume
```

#### Use Harmony integration and SingleR annotation
```bash
nextflow run main.nf \
    --samplesheet metadata/samples.csv \
    --integration_method Harmony \
    --annotation_method SingleR \
    --singler_ref HumanPrimaryCellAtlasData \
    -resume
```

#### Custom clustering resolution
```bash
nextflow run main.nf \
    --samplesheet metadata/samples.csv \
    --resolution 0.5 \
    -resume
```

## Output Files

### Directory Structure

```
results/
├── fastqc/                     # FastQC reports
│   └── {sample_id}_fastqc.html
├── multiqc/                    # Aggregated QC report
│   └── multiqc_report.html
├── cellranger/                 # Cell Ranger outputs
│   └── {sample_id}_output/
│       └── outs/
│           ├── filtered_feature_bc_matrix/
│           ├── web_summary.html
│           └── metrics_summary.csv
├── integration/                # Integration results
│   ├── integrated_seurat_res{resolution}.rds
│   ├── qc_plots/
│   │   ├── umap.jpg
│   │   ├── tsne.jpg
│   │   ├── feature_plot.jpg
│   │   └── qc_metrics.jpg
│   └── integration_summary.txt
├── annotation/                 # Annotation results
│   ├── annotated_seurat.rds
│   ├── annotation_plots/
│   │   ├── umap_annotated.png
│   │   ├── cluster_markers.png
│   │   └── feature_plots.png
│   ├── cell_type_markers.csv
│   └── annotation_summary.txt
└── pipeline_info/              # Pipeline execution reports
    ├── execution_timeline_*.html
    ├── execution_report_*.html
    ├── execution_trace_*.txt
    └── pipeline_dag_*.html
```

### Key Output Files

| File | Description |
|------|-------------|
| [multiqc/multiqc_report.html](results/multiqc/multiqc_report.html) | Aggregated quality control report |
| [integration/integrated_seurat_res0.3.rds](results/integration/integrated_seurat_res0.3.rds) | Integrated Seurat object with clustering |
| [annotation/annotated_seurat.rds](results/annotation/annotated_seurat.rds) | Final annotated Seurat object with cell types |
| [annotation/cell_type_markers.csv](results/annotation/cell_type_markers.csv) | Marker genes for each cell type/cluster |
| [pipeline_info/execution_report_*.html](results/pipeline_info/) | Resource usage and execution statistics |

## Advanced Options

### Integration Methods

The pipeline supports four integration methods:

1. **CCA (Canonical Correlation Analysis)** - Default, recommended for most datasets
2. **RPCA (Reciprocal PCA)** - Better for large datasets with shared structure
3. **Harmony** - Fast, suitable for large datasets with batch effects
4. **FastMNN** - Mutual nearest neighbors, good for complex batch effects

```bash
# Example using Harmony
nextflow run main.nf \
    --samplesheet metadata/samples.csv \
    --integration_method Harmony \
    -resume
```

### Annotation Methods

#### Azimuth (Default)
- Reference-based annotation using pre-computed references
- Available references: `pbmc`, `pbmc_multimodal`, `bonemarrow`, `kidney`, etc.

```bash
nextflow run main.nf \
    --annotation_method Azimuth \
    --reference pbmc \
    -resume
```

#### SingleR
- Correlation-based annotation using reference datasets from CellDex
- Available references: `HumanPrimaryCellAtlasData`, `BlueprintEncodeData`, `MonacoImmuneData`, etc.

```bash
nextflow run main.nf \
    --annotation_method SingleR \
    --singler_ref HumanPrimaryCellAtlasData \
    -resume
```

### Quality Control Filters

Adjust QC thresholds in [nextflow.config](nextflow.config):

```groovy
params {
    min_features = 200      // Min genes per cell
    max_features = 7000     // Max genes per cell (filter doublets)
    min_umi = 1000          // Min UMI counts per cell
    max_umi = 50000         // Max UMI counts per cell
    max_mito = 20           // Max mitochondrial percentage
    min_cells = 3           // Min cells per gene
}
```

## Interactive Visualization

The pipeline includes a Shiny app for interactive exploration of results.

### Launch Shiny App

```bash
# Activate conda environment
conda activate /path/to/conda/env

# Run Shiny app
Rscript Rscripts/shiny_app.r
```

### Features

- Interactive UMAP and t-SNE plots
- Feature expression visualization
- Cluster and cell type exploration
- Differential expression analysis
- Marker gene identification
- Customizable plot parameters

## Troubleshooting

### Common Issues

#### 1. Cell Ranger not found
```
Error: cellranger command not found
```

**Solution:** Ensure Cell Ranger is loaded via module system:
```bash
module load cellranger/8.0.1
```

#### 2. GSL library error (Azimuth)
```
Error: cannot open shared object file: libgsl.so.25
```

**Solution:** Create symbolic link as shown in Installation step 2:
```bash
cd /path/to/conda/env/lib
ln -s libgsl.so.27 libgsl.so.25
```

#### 3. Out of memory during Cell Ranger
```
Error: SIGKILL (memory limit exceeded)
```

**Solution:** Increase memory allocation in [nextflow.config](nextflow.config):
```groovy
process {
    withName: CELLRANGER_COUNT {
        memory = '150 GB'  // Increase from 120 GB
    }
}
```

#### 4. Integration fails with low cell counts
```
Warning: Not enough cells for integration
```

**Solution:** Lower QC thresholds or adjust integration parameters:
```bash
nextflow run main.nf \
    --min_features 100 \
    --min_umi 500 \
    -resume
```

#### 5. Nextflow resume not working
```
Error: Cannot resume previous execution
```

**Solution:** Clear work directory and restart:
```bash
rm -rf work/
nextflow run main.nf --samplesheet metadata/samples.csv
```

### Getting Help

- Check Nextflow execution reports in `results/pipeline_info/`
- Review process logs in `work/` directories
- Examine R script outputs in integration/annotation directories
- Ensure all paths in [nextflow.config](nextflow.config) are correct

## Project Structure

```
nf-scrna-seq/
├── main.nf                     # Main Nextflow workflow
├── nextflow.config             # Pipeline configuration
├── run_pipeline.sh             # SLURM submission script
├── env.yaml                    # Conda environment specification
│
├── modules/                    # Nextflow process modules
│   ├── qc.nf                  # FastQC quality control
│   ├── multiqc.nf             # MultiQC aggregation
│   ├── cellranger.nf          # Cell Ranger alignment
│   ├── integration.nf         # Seurat integration
│   └── annotation.nf          # Cell type annotation
│
├── Rscripts/                   # R analysis scripts
│   ├── integration.r          # Integration & clustering
│   ├── annotation.r           # Cell type annotation
│   └── shiny_app.r           # Interactive visualization
│
├── metadata/                   # Sample metadata
│   ├── samples.csv            # Sample sheet
│   ├── download.sh            # Data download script
│   └── OSA_SraRunTable.csv    # SRA metadata
│
└── results/                    # Pipeline outputs
    ├── fastqc/
    ├── multiqc/
    ├── cellranger/
    ├── integration/
    ├── annotation/
    └── pipeline_info/
```

## Citation

If you use this pipeline in your research, please cite:

**Nextflow:**
- Di Tommaso P, et al. (2017) Nextflow enables reproducible computational workflows. Nat Biotechnol. 35(4):316-319.

**Cell Ranger:**
- 10x Genomics. Cell Ranger. https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

**Seurat:**
- Hao Y, et al. (2021) Integrated analysis of multimodal single-cell data. Cell. 184(13):3573-3587.e29.

**Azimuth:**
- Hao Y, et al. (2021) Integrated analysis of multimodal single-cell data. Cell. 184(13):3573-3587.e29.

**SingleR:**
- Aran D, et al. (2019) Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nat Immunol. 20(2):163-172.

## License

This pipeline is distributed under the MIT License. See LICENSE file for details.

## Contact

For questions or issues, please open an issue on the GitHub repository or contact the maintainer.

---

**Last updated:** January 2026
