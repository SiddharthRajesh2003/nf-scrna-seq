#!/usr/bin/env Rscript

# Cell Type Annotation Script for scRNA-seq data (Seurat v5)
# Supports multiple annotation methods: SingleR, Azimuth, clustifyr

# Load renv
renv::load("/N/project/Krolab/Siddharth/Pipelines/scrna-seq/Renv")

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

option_list <- list(
  make_option(c("--seurat_object"), type = "character", default = NULL,
              help = "Path to integrated Seurat object RDS file"),
  make_option(c("--output"), type = "character", default = "annotated_seurat_object.rds",
              help = "Output RDS file path"),
  make_option(c("--method"), type = "character", default = "singler",
              help = "Annotation method: singler, azimuth, or clustifyr"),
  make_option(c("--reference"), type = "character", default = "celldex::HumanPrimaryCellAtlasData",
              help = "Reference dataset for annotation"),
  make_option(c("--species"), type = "character", default = "human",
              help = "Species: human or mouse"),
  make_option(c("--azimuth_reference"), type = "character", default = "pbmc",
              help = "Azimuth reference: pbmc, lung, motor_cortex, etc.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$seurat_object)) {
  print_help(opt_parser)
  stop("--seurat_object must be provided", call. = FALSE)
}

# Create output directory for plots
if (!dir.exists("annotation_plots")) dir.create("annotation_plots")

# Load Seurat Object
seurat_obj <- readRDS(opt$seurat_object)

cat("Loaded object with", ncol(seurat_obj), "cells\n")
cat("Number of clusters:", length(unique(seurat_obj$seurat_clusters)), "\n\n")

# ============================================================================
# Method 1: SingleR Annotation
# ============================================================================
annotate_with_singler <- function(seurat_obj, reference_string) {
  cat("Running SingleR annotation...\n")
  suppressPackageStartupMessages({
    library(SingleR)
    library(celldex)
  })

  # Parse reference dataset
  if (grepl("celldex::", reference_string)) {
    ref_name <- gsub("celldex::", "", reference_string)
    cat("Using celldex reference:", ref_name, "\n")

    # Get reference data
    ref_data <- switch(ref_name,
        "HumanPrimaryCellAtlasData" = celldex::HumanPrimaryCellAtlasData(),
        "BlueprintEncodeData" = celldex::BlueprintEncodeData(),
        "DatabaseImmuneCellExpressionData" = celldex::DatabaseImmuneCellExpressionData(),
        "MonacoImmuneData" = celldex::MonacoImmuneData(),
        "MouseRNAseqData" = celldex::MouseRNAseqData(),
        "ImmGenData" = celldex::ImmGenData(),
        stop("Unknown celldex reference: ", ref_name)
    )
  } else {
    cat("Loading custom reference from:", reference_string, "\n")
    ref_data <- readRDS(reference_string)
  }

  cat("Running SingleR prediction...\n")

  test_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

  predictions <- SingleR(
    test = test_data,
    ref = ref_data,
    labels = ref_data$label.main
  )

  # Add predictions to Seurat Object
  seurat_obj$singler_labels <- predictions$labels
  seurat_obj$singler_scores <- predictions$scores
  seurat_obj$singler_pruned_labels <- predictions$pruned.labels

  # Also get fine-grained labels if available
  if ("label.fine" %in% colnames(colData(ref_data))) {
    predictions_fine <- SingleR(
      test = test_data,
      ref = ref_data,
      labels = ref_data$label.fine
    )
    seurat_obj$singler_labels_fine <- predictions_fine$labels
  }

  # Summary
  cat("\nSingleR annotation summary:\n")
  print(table(seurat_obj$singler_labels))

  # Summary
  cat("\nSingleR annotation summary:\n")
  print(table(seurat_obj$singler_labels))
  
  # Generate plots
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "singler_labels", 
                label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle("SingleR Cell Type Annotation") +
    theme(legend.position = "bottom")
  png("annotation_plots/singler_umap.png", width = 12, height = 10)
  print(p1)
  dev.off()
  
  # Heatmap of scores
  p2 <- plotScoreHeatmap(predictions)
  ggsave("annotation_plots/singler_score_heatmap.png", width = 12, height = 10)
  print(p2)
  dev.off()
  
  # Cluster vs cell type
  cluster_celltype <- table(seurat_obj$seurat_clusters, seurat_obj$singler_labels)
  write.csv(cluster_celltype, "annotation_plots/cluster_celltype_table.csv")
  
  return(seurat_obj)
}

# ============================================================================
# Method 2: Azimuth Annotation
# ============================================================================
annotate_with_azimuth <- function(seurat_obj, reference_name) {
  cat("Running Azimuth annotation...\n")
  
  suppressPackageStartupMessages({
    library(Azimuth)
  })
  
  cat("Using Azimuth reference:", reference_name, "\n")
  seurat_obj <- RunAzimuth(
    seurat_obj,
    reference = reference_name,
    assay = 'RNA'
  )

  # Summary
  cat("\nAzimuth annotation summary:\n")
  if ("predicted.celltype.l1" %in% colnames(seurat_obj@meta.data)) {
    print(table(seurat_obj$predicted.celltype.l1))
    
    # Generate plots
    p1 <- DimPlot(seurat_obj, reduction = "umap", 
                  group.by = "predicted.celltype.l1", 
                  label = TRUE, repel = TRUE, label.size = 3) +
      ggtitle("Azimuth Cell Type Annotation (Level 1)") +
      theme(legend.position = "bottom")
    png("annotation_plots/azimuth_l1_umap.png", width = 12, height = 10)
    print(p1)
    dev.off()
  }
  
  if ("predicted.celltype.l2" %in% colnames(seurat_obj@meta.data)) {
    print(table(seurat_obj$predicted.celltype.l2))
    
    p2 <- DimPlot(seurat_obj, reduction = "umap", 
                  group.by = "predicted.celltype.l2", 
                  label = TRUE, repel = TRUE, label.size = 3) +
      ggtitle("Azimuth Cell Type Annotation (Level 2)") +
      theme(legend.position = "bottom")
    ggsave("annotation_plots/azimuth_l2_umap.png", width = 14, height = 12)
    print(p2)
  }
  
  return(seurat_obj)
}

# ============================================================================
# Find cluster markers
# ============================================================================
find_cluster_markers <- function(seurat_obj) {
  cat("\nFinding cluster markers...\n")
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Find all markers
  all_markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  # Save top markers
  top_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  write.csv(top_markers, "cell_type_markers.csv", row.names = FALSE)
  
  # Generate heatmap
  top5 <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
  
  p <- DoHeatmap(seurat_obj, features = top5$gene) +
    scale_fill_gradientn(colors = c("blue", "white", "red"))
  ggsave("annotation_plots/marker_heatmap.png", p, width = 14, height = 12)
  
  return(all_markers)
}

# ============================================================================
# Main execution
# ============================================================================

cat("Starting cell type annotation with method:", opt$method, "\n\n")

# Run annotation based on selected method
seurat_obj <- switch(opt$method,
  "singler" = annotate_with_singler(seurat_obj, opt$reference),
  "azimuth" = annotate_with_azimuth(seurat_obj, opt$azimuth_reference),
  stop("Unknown annotation method: ", opt$method)
)

# Find cluster markers
markers <- find_cluster_markers(seurat_obj)

# Generate comparison plots
cat("\nGenerating comparison plots...\n")

# Original clusters vs annotations
if (opt$method == "singler") {
  p_comp <- DimPlot(seurat_obj, reduction = "umap", 
                    group.by = c("seurat_clusters", "singler_labels"), 
                    ncol = 2)
  png("annotation_plots/clusters_vs_annotations.png", 
        width = 16, height = 6)
  print(p_comp)
  dev.off()
} else if (opt$method == "azimuth" && 
          "predicted.celltype.l1" %in% colnames(seurat_obj@meta.data)) {
  p_comp <- DimPlot(seurat_obj, reduction = "umap", 
                    group.by = c("seurat_clusters", "predicted.celltype.l1"), 
                    ncol = 2)
  ggsave("annotation_plots/clusters_vs_annotations.png",
          width = 16, height = 6)
  print(p_comp)
  dev.off()
}

# Split by sample
p_split <- DimPlot(seurat_obj, reduction = "umap", 
                    group.by = "sample", split.by = "sample", ncol = 2)
ggsave("annotation_plots/umap_by_sample_split.png",
       width = 14, height = ceiling(length(unique(seurat_obj$sample))/2) * 6)
print(p_split)
dev.off()

# Save annotated object
cat("\nSaving annotated Seurat object...\n")
saveRDS(seurat_obj, file = opt$output)

cat("\nCell type annotation complete!\n")
cat("Annotated object saved to:", opt$output, "\n")
cat("Plots saved to: annotation_plots/\n")
cat("Marker genes saved to: cell_type_markers.csv\n")

# Print final summary
cat("\n=== Final Summary ===\n")
cat("Total cells:", ncol(seurat_obj), "\n")
cat("Total features:", nrow(seurat_obj), "\n")
cat("Number of samples:", length(unique(seurat_obj$sample)), "\n")
cat("Number of clusters:", length(unique(seurat_obj$seurat_clusters)), "\n")

if (opt$method == "singler") {
  cat("Number of cell types (SingleR):", 
      length(unique(seurat_obj$singler_labels)), "\n")
} else if (opt$method == "azimuth") {
  if ("predicted.celltype.l1" %in% colnames(seurat_obj@meta.data)) {
    cat("Number of cell types (Azimuth L1):", 
        length(unique(seurat_obj$predicted.celltype.l1)), "\n")
  }
} else if (opt$method == "clustifyr") {
  cat("Number of cell types (clustifyr):", 
      length(unique(seurat_obj$clustifyr_type)), "\n")
}

cat("\nPipeline complete!\n")