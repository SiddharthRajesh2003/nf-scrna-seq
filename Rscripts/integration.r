#!/usr/bin/env Rscript

# === Environment Setup ===
Sys.setenv(
  LANG = "en_US.UTF-8",
  LC_ALL = "en_US.UTF-8"
)

# === Load Required Libraries ===
suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(optparse)
    library(ggplot2)
    library(patchwork)
})

cat("✓ All libraries loaded successfully\n\n")

option_list <- list(
  make_option(c("--input_dirs"), type = "character", default = NULL,
            help = "Comma separated list of Cellranger output directories"),
  make_option(c("--output"), type = "character", default = "integrated_seurat.rds",
            help = "Output RDS file path [default: %default]"),
  make_option(c("--min_cells"), type = "integer", default = 3,
            help = "Minimum cells per gene for Seurat::Read10X function [default: %default]"),
  make_option(c("--min_features"), type = "integer", default = 200,
            help = "Minimum features per cell for Seurat::CreateSeuratObject function [default: %default]"),
  make_option(c("--max_features"), type = "integer", default = 7000,
              help = "Maximum number of features per cell [default: %default]"),
  make_option(c("--max_mt_percent"), type = "numeric", default = 20,
              help = "Maximum mitochondrial percentage [default: %default]"),
  make_option(c("--dims"), type = "integer", default = 30,
              help = "Number of dimensions for integration [default: %default]"),
  make_option(c("--resolution"), type = "numeric", default = 0.3,
              help = "Clustering resolution [default: %default]"),
  make_option(c("--max_umi"), type = "integer", default = 50000,
              help = "Maximum Number of UMIs per cell [default: %default]"),
  make_option(c("--min_umi"), type = "integer", default = 1000,
              help = "Minimum number of UMIs per cell [default: %default]"),
  make_option(c("--integration_method"), type = "character", default = 'cca',
              help = "Integration method: cca[CCAIntegration], rpca[RPCAIntegration], harmony[HarmonyIntegration], mnn[FastMNNIntegration] [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_dirs)) {
  print_help(opt_parser)
  stop("--input_dirs must be provided", call. = FALSE)
}

input_dirs <- strsplit(opt$input_dirs, ",")[[1]]
sample_names <- basename(dirname(dirname(input_dirs)))
sample_names <- gsub("_output$", "", sample_names)

cat("Processing", length(input_dirs), "samples\n")
cat(paste(sample_names, collapse = "\n"), "\n\n")

# Create output directory for plots
if (!dir.exists("qc_plots")) dir.create("qc_plots")

# ============================================================================
# Load and create Seurat objects (Seurat v5 format)
# ============================================================================
cat("Loading data from Cell Ranger outputs...\n")

seurat_list <- list()

for (i in seq_along(input_dirs)) {
  cat("Loading sample:", sample_names[i], "\n")

  counts <- Read10X(data.dir = input_dirs[i])

  seurat_obj <- CreateSeuratObject(counts = counts,
                                    project = sample_names[i],
                                    min.cells = opt$min_cells,
                                    min.features = opt$min_features)
  
  seurat_obj$sample <- sample_names[i]

  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  
  seurat_list[[sample_names[i]]] <- seurat_obj
}

cat("Filtering cells based on QC metrics\n")
seurat_list <- lapply(seurat_list, function(x) {
    subset(x,
           subset = nFeature_RNA > opt$min_features &
                    nFeature_RNA < opt$max_features &
                    nCount_RNA > opt$min_umi &
                    nCount_RNA < opt$max_umi &
                    percent.mt < opt$max_mt_percent &
                    percent.rb < 40)
})

# Print QC stats
for (name in names(seurat_list)) {
  cat(sprintf("Sample %s: %d cells retained after QC\n", 
              name, ncol(seurat_list[[name]])))
}

cat('Merging Seurat Objects...\n')
merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[2:length(seurat_list)],
  add.cell.ids = names(seurat_list),
  project = 'Integrated_scrna'
) 

cat("Merged object contains", ncol(merged_seurat), "cells\n\n")

cat("Starting Seurat v5 integration workflow...\n")
cat("Method:", opt$integration_method, "\n\n")

if (inherits(merged_seurat[["RNA"]], "Assay5")) {
  merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])
}

merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'orig.ident', into = c("Condition", "Sample No."), sep = '_')

cat("Normalizing data...\n")
merged_seurat <- NormalizeData(merged_seurat)

cat("Finding variable features...\n")
merged_seurat <- FindVariableFeatures(merged_seurat, selection.method = "vst", nfeatures = 2000)

cat("Scaling data...\n")
merged_seurat <- ScaleData(merged_seurat)

cat("Running PCA...\n")
merged_seurat <- RunPCA(merged_seurat)

# Integrate layers based on chosen method
cat("\nIntegrating layers with", opt$integration_method, "...\n")

if (opt$integration_method == 'cca') {
  seurat_integrated <- IntegrateLayers(
    object = merged_seurat,
    method = CCAIntegration,
    orig.reduction = 'pca',
    new.reduction = 'integrated.cca',
    verbose = TRUE
  )
  reduction_name <- 'integrated.cca'

} else if (opt$integration_method == 'rpca') {
  seurat_integrated <- IntegrateLayers(
    object = merged_seurat,
    method = RPCAIntegration,
    orig.reduction = 'pca',
    new.reduction = 'integrated.rpca',
    verbose = TRUE
  )
  reduction_name <- 'integrated.rpca'

} else if (opt$integration_method == 'harmony') {
  seurat_integrated <- IntegrateLayers(
    object = merged_seurat,
    method = HarmonyIntegration,
    orig.reduction = 'pca',
    new.reduction = 'harmony',
    verbose = TRUE
  )
  reduction_name <- 'harmony'

} else if (opt$integration_method == 'mnn') {
  seurat_integrated <- IntegrateLayers(
    object = merged_seurat,
    method = FastMNNIntegration,
    orig.reduction = 'pca',
    new.reduction = 'integrated.mnn',
    verbose = TRUE
  )
  reduction_name <- 'integrated.mnn'
} else {
  stop("Unknown integration method: ", opt$integration_method) 
}

# Rejoin layers after integration
cat("\nRejoining layers...\n")
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])

# ============================================================================
# Dimensionality Reduction and Clustering
# ============================================================================

cat("\nPerforming dimensionality reduction and clustering...\n")

cat("Finding Neighbours...\n")
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction_name, dims = 1:opt$dims)

cat("Finding Clusters...\n")
seurat_integrated <- FindClusters(seurat_integrated, resolution = opt$resolution)

cat("Running UMAP...\n")
seurat_integrated <- RunUMAP(seurat_integrated, reduction = reduction_name, dims = 1:opt$dims)

cat("Running TSNE...\n")
seurat_integrated <- RunTSNE(seurat_integrated, reduction = reduction_name, dims = 1:opt$dims, perplexity = 30, check_duplicates = FALSE)

# ============================================================================
# Generate QC Plots
# ============================================================================

cat("\nGenerating QC plots...\n")

# UMAP by Sample
p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "sample") +
      ggtitle("UMAP by Sample") +
      theme_minimal() +
      theme(legend.position = 'right')
jpeg('qc_plots/umap_by_sample.jpg', width = 12, height = 8, units = 'in', res = 300)
print(p1)
dev.off()

# UMAP by Clusters
p2 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
      ggtitle("UMAP by Clusters") +
      theme_minimal() +
      theme(legend.position = 'right')
jpeg('qc_plots/umap_by_clusters.jpg', width = 12, height = 8, units = 'in', res = 300)
print(p2)
dev.off()

# Split by sample
p3 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters", split.by = 'sample', ncol = 2, label = TRUE, pt.size = 0.5) +
      ggtitle("UMAP by Clusters Split by Sample") +
      theme_minimal() +
      theme(legend.position = 'right')
jpeg('qc_plots/umap_by_clusters_split_by_sample.jpg', width = 16, height = 12, units = 'in', res = 300)
print(p3)
dev.off()

# TSNE by sample
p4 <- DimPlot(seurat_integrated, reduction = "tsne", group.by = "sample") +
      ggtitle("t-SNE by Sample") +
      theme_minimal() +
      theme(legend.position = 'right')
jpeg('qc_plots/tsne_by_sample.jpg', width = 12, height = 8, units = 'in', res = 300)
print(p4)
dev.off()

# TSNE by Clusters
p5 <- DimPlot(seurat_integrated, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
      ggtitle("t-SNE by Clusters") +
      theme_minimal() +
      theme(legend.position = 'right')
jpeg('qc_plots/tsne_by_clusters.jpg', width = 12, height = 8, units = 'in', res = 300)
print(p5)
dev.off()

# Split by sample
p6 <- DimPlot(seurat_integrated, reduction = "tsne", group.by = "seurat_clusters", split.by = 'sample', ncol = 2, label = TRUE, pt.size = 0.5) +
      ggtitle("t-SNE by Clusters Split by Sample") +
      theme_minimal() +
      theme(legend.position = 'right')
jpeg('qc_plots/tsne_by_clusters_split_by_sample.jpg', width = 16, height = 12, units = 'in', res = 300)
print(p6)
dev.off()

# QC metrics violin plot
p7 <- VlnPlot(seurat_integrated, 
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
              ncol = 3, pt.size = 0, group.by = "sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
jpeg('qc_plots/qc_metrics_by_sample.jpg', width = 16, height = 6, units = 'in', res = 300)
print(p7)
dev.off()

# Feature plots for key markers
common_markers <- c("CD3D", "CD8A", "CD4", "MS4A1", "CD14", "FCGR3A")
available_markers <- common_markers[common_markers %in% rownames(seurat_integrated)]

if (length(available_markers) > 0) {
  p8 <- FeaturePlot(seurat_integrated, features = available_markers, ncol = 3)
  jpeg("qc_plots/marker_genes.jpg", width = 15, height = 10, units = 'in', res = 300)
  print(p8)
  dev.off()
}

# PCA plot
p9 <- DimPlot(seurat_integrated, reduction = "pca", group.by = "sample") +
  ggtitle("PCA by Sample") +
  theme_minimal()
jpeg("qc_plots/pca_by_sample.jpg", width = 10, height = 8, units = 'in', res = 300)
print(p9)
dev.off()

# Elbow plot
p10 <- ElbowPlot(seurat_integrated, ndims = opt$dims) +
  ggtitle("PCA Elbow Plot") +
  theme_minimal()
jpeg("qc_plots/elbow_plot.jpg", width = 8, height = 6, units = 'in', res = 300)
print(p10)
dev.off()

# ============================================================================
# Save Results
# ============================================================================

cat("\nSaving Integrated Seurat Object...\n")
saveRDS(seurat_integrated, file = opt$output)

cat("\nIntegration complete!\n")
cat("Integrated object saved to:", opt$output, "\n")

# ============================================================================
# Print Summary to file
# ============================================================================

sink("integration_summary.txt")
cat(rep("=", 60), "\n", sep = "")
cat("INTEGRATION SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Integration method:", opt$integration_method, "\n")
cat("Reduction used:", reduction_name, "\n")
cat("\nData summary:\n")
cat("  Total cells:", ncol(seurat_integrated), "\n")
cat("  Total features:", nrow(seurat_integrated), "\n")
cat("  Number of samples:", length(unique(seurat_integrated$sample)), "\n")
cat("  Number of clusters:", length(unique(seurat_integrated$seurat_clusters)), "\n")
cat("\nCells per sample:\n")
print(table(seurat_integrated$sample))
cat("\nCells per cluster:\n")
print(table(seurat_integrated$seurat_clusters))
cat("\nQC metrics (median):\n")
cat("  Features per cell:", median(seurat_integrated$nFeature_RNA), "\n")
cat("  UMIs per cell:", median(seurat_integrated$nCount_RNA), "\n")
cat("  Mitochondrial %:", round(median(seurat_integrated$percent.mt), 2), "\n")
cat("\nOutput files:\n")
cat("  RDS object:", opt$output, "\n")
cat("  QC plots: qc_plots/\n")
cat(rep("=", 60), "\n", sep = "")
sink()

cat("\nPipeline complete!\n")