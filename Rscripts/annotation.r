#!/usr/bin/env Rscript

# === Environment Setup ===
Sys.setenv(
  LANG = "en_US.UTF-8",
  LC_ALL = "en_US.UTF-8"
)

# Paths
base_dir <- "/N/project/Krolab/Siddharth/Pipelines/scrna-seq"
renv_project <- file.path(base_dir, "Renv")
conda_lib <- file.path(base_dir, "conda_envs/env/renv/lib/R/library")

# CRITICAL: Disable renv sandbox to prevent home directory access
Sys.setenv(
  RENV_CONFIG_SANDBOX_ENABLED = "FALSE",
  RENV_PATHS_ROOT = file.path(base_dir, "renv_root"),
  RENV_PATHS_CACHE = file.path(base_dir, "renv_cache"),
  RENV_PATHS_LIBRARY_ROOT = file.path(renv_project, "renv/library"),
  RENV_PATHS_LIBRARY = file.path(renv_project, "renv/library"),
  RENV_PATHS_SANDBOX = file.path(base_dir, "renv_cache/sandbox")
)

# Set library paths BEFORE loading renv
if (dir.exists(conda_lib)) {
  .libPaths(c(conda_lib, .libPaths()))
  cat("✓ Added conda library path\n")
}

# Try to load renv (but it's optional)
tryCatch({
  if (file.exists(file.path(renv_project, "renv/activate.R"))) {
    source(file.path(renv_project, "renv/activate.R"))
    cat("✓ renv loaded\n")
  }
}, error = function(e) {
  cat("Note: renv load had issues, continuing with conda libraries\n")
})

# Ensure both paths are available (renv first, then conda)
renv_lib <- file.path(renv_project, "renv/library/linux-rhel-8.10/R-4.5/x86_64-conda-linux-gnu")
if (dir.exists(renv_lib) && !renv_lib %in% .libPaths()) {
  .libPaths(c(renv_lib, .libPaths()))
}
if (dir.exists(conda_lib) && !conda_lib %in% .libPaths()) {
  .libPaths(c(.libPaths(), conda_lib))
}

# Verify library paths
cat("\nLibrary search paths:\n")
for (i in seq_along(.libPaths())) {
  cat(sprintf("  %d. %s\n", i, .libPaths()[i]))
}
cat("\n")

# Check where Seurat is located
if (file.exists(file.path(conda_lib, "Seurat"))) {
  cat("✓ Seurat found in conda library\n")
}
if (dir.exists(renv_lib) && file.exists(file.path(renv_lib, "Seurat"))) {
  cat("✓ Seurat found in renv library\n")
}
cat("\n")

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
  make_option(c("--output"), type = "character", default = "annotated_seurat.rds",
              help = "Output RDS file path [default: %default]"),
  make_option(c("--method"), type = "character", default = "singler",
              help = "Annotation method: singler or azimuth [default: %default]"),
  make_option(c("--reference"), type = "character", default = "celldex::HumanPrimaryCellAtlasData",
              help = "Reference dataset for annotation [default: %default]"),
  make_option(c("--species"), type = "character", default = "human",
              help = "Species: human or mouse [default: %default]"),
  make_option(c("--azimuth_reference"), type = "character", default = "pbmc",
              help = "Azimuth reference: pbmc, lung, motor_cortex, etc. [default: %default]")
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
cat("Loading Seurat object from:", opt$seurat_object, "\n")
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

  # Get normalized data
  test_data <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")

  # Run SingleR with main labels
  predictions <- SingleR(
    test = test_data,
    ref = ref_data,
    labels = ref_data$label.main
  )

  # Add predictions to Seurat Object
  seurat_obj$singler_labels <- predictions$labels
  seurat_obj$singler_pruned_labels <- predictions$pruned.labels

  # Also get fine-grained labels if available
  if ("label.fine" %in% colnames(colData(ref_data))) {
    cat("Running SingleR with fine labels...\n")
    predictions_fine <- SingleR(
      test = test_data,
      ref = ref_data,
      labels = ref_data$label.fine
    )
    seurat_obj$singler_labels_fine <- predictions_fine$labels
  }

  # Summary
  cat("\nSingleR annotation summary (main labels):\n")
  print(table(seurat_obj$singler_labels))
  
  # Generate plots
  cat("Generating SingleR plots...\n")
  
  # UMAP with cell types
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "singler_labels", 
                label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.5) +
    ggtitle("SingleR Cell Type Annotation") +
    theme_minimal() +
    theme(legend.position = "right")
  ggsave("annotation_plots/singler_umap.png", p1, width = 12, height = 10, dpi = 300)
  
  # Score heatmap
  png("annotation_plots/singler_score_heatmap.png", width = 12, height = 10, res = 300)
  plotScoreHeatmap(predictions)
  dev.off()
  
  # Cluster vs cell type table
  cluster_celltype <- table(Cluster = seurat_obj$seurat_clusters, 
                            CellType = seurat_obj$singler_labels)
  write.csv(cluster_celltype, "annotation_plots/cluster_celltype_table.csv")
  
  # Stacked bar plot of cell types per cluster
  celltype_df <- data.frame(
    cluster = seurat_obj$seurat_clusters,
    celltype = seurat_obj$singler_labels
  )
  
  p2 <- ggplot(celltype_df, aes(x = cluster, fill = celltype)) +
    geom_bar(position = "fill") +
    labs(title = "Cell Type Composition per Cluster",
         x = "Cluster", y = "Proportion", fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("annotation_plots/celltype_composition_per_cluster.png", p2, 
         width = 12, height = 8, dpi = 300)
  
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
  
  # Run Azimuth
  seurat_obj <- RunAzimuth(
    seurat_obj,
    reference = reference_name,
    assay = 'RNA'
  )

  # Summary
  cat("\nAzimuth annotation summary:\n")
  if ("predicted.celltype.l1" %in% colnames(seurat_obj@meta.data)) {
    cat("\nLevel 1 annotations:\n")
    print(table(seurat_obj$predicted.celltype.l1))
    
    # Generate plots
    p1 <- DimPlot(seurat_obj, reduction = "umap", 
                  group.by = "predicted.celltype.l1", 
                  label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.5) +
      ggtitle("Azimuth Cell Type Annotation (Level 1)") +
      theme_minimal() +
      theme(legend.position = "right")
    ggsave("annotation_plots/azimuth_l1_umap.png", p1, width = 12, height = 10, dpi = 300)
    
    # Composition plot
    celltype_df <- data.frame(
      cluster = seurat_obj$seurat_clusters,
      celltype = seurat_obj$predicted.celltype.l1
    )
    
    p_comp <- ggplot(celltype_df, aes(x = cluster, fill = celltype)) +
      geom_bar(position = "fill") +
      labs(title = "Cell Type Composition per Cluster (Level 1)",
           x = "Cluster", y = "Proportion", fill = "Cell Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("annotation_plots/azimuth_l1_composition.png", p_comp, 
           width = 12, height = 8, dpi = 300)
  }
  
  if ("predicted.celltype.l2" %in% colnames(seurat_obj@meta.data)) {
    cat("\nLevel 2 annotations:\n")
    print(table(seurat_obj$predicted.celltype.l2))
    
    p2 <- DimPlot(seurat_obj, reduction = "umap", 
                  group.by = "predicted.celltype.l2", 
                  label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.5) +
      ggtitle("Azimuth Cell Type Annotation (Level 2)") +
      theme_minimal() +
      theme(legend.position = "right")
    ggsave("annotation_plots/azimuth_l2_umap.png", p2, width = 14, height = 12, dpi = 300)
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
  
  # Save all markers
  write.csv(all_markers, "cell_type_markers.csv", row.names = FALSE)
  
  # Save top markers
  top_markers <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  write.csv(top_markers, "annotation_plots/top10_markers_per_cluster.csv", row.names = FALSE)
  
  # Generate heatmap of top markers
  top5 <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC)
  
  if (nrow(top5) > 0) {
    p <- DoHeatmap(seurat_obj, features = unique(top5$gene), size = 3) +
      scale_fill_gradientn(colors = c("blue", "white", "red")) +
      theme(axis.text.y = element_text(size = 8))
    ggsave("annotation_plots/marker_heatmap.png", p, width = 14, height = 12, dpi = 300)
  }
  
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
cat("\nFinding differential markers...\n")
markers <- find_cluster_markers(seurat_obj)

# Generate comparison plots
cat("\nGenerating comparison plots...\n")

# Original clusters vs annotations
if (opt$method == "singler") {
  p_comp <- DimPlot(seurat_obj, reduction = "umap", 
                    group.by = c("seurat_clusters", "singler_labels"), 
                    ncol = 2, pt.size = 0.5) +
    plot_annotation(title = "Clusters vs SingleR Annotations")
  ggsave("annotation_plots/clusters_vs_annotations.png", p_comp,
        width = 16, height = 6, dpi = 300)
  
} else if (opt$method == "azimuth" && 
          "predicted.celltype.l1" %in% colnames(seurat_obj@meta.data)) {
  p_comp <- DimPlot(seurat_obj, reduction = "umap", 
                    group.by = c("seurat_clusters", "predicted.celltype.l1"), 
                    ncol = 2, pt.size = 0.5) +
    plot_annotation(title = "Clusters vs Azimuth Annotations")
  ggsave("annotation_plots/clusters_vs_annotations.png", p_comp,
          width = 16, height = 6, dpi = 300)
}

# Split by sample if multiple samples exist
if (length(unique(seurat_obj$sample)) > 1) {
  n_samples <- length(unique(seurat_obj$sample))
  n_cols <- min(3, n_samples)
  n_rows <- ceiling(n_samples / n_cols)
  
  annotation_col <- if (opt$method == "singler") {
    "singler_labels"
  } else if ("predicted.celltype.l1" %in% colnames(seurat_obj@meta.data)) {
    "predicted.celltype.l1"
  } else {
    "seurat_clusters"
  }
  
  p_split <- DimPlot(seurat_obj, reduction = "umap", 
                      group.by = annotation_col, 
                      split.by = "sample", 
                      ncol = n_cols,
                      label = TRUE,
                      repel = TRUE,
                      pt.size = 0.3) +
    theme_minimal()
  
  ggsave("annotation_plots/umap_by_sample_split.png", p_split,
         width = n_cols * 5, height = n_rows * 4.5, dpi = 300)
}

# Save annotated object
cat("\nSaving annotated Seurat object...\n")
saveRDS(seurat_obj, file = opt$output)

# ============================================================================
# Generate summary report
# ============================================================================
sink("annotation_summary.txt")
cat(rep("=", 70), "\n", sep = "")
cat("CELL TYPE ANNOTATION SUMMARY\n")
cat(rep("=", 70), "\n", sep = "")
cat("Annotation method:", opt$method, "\n")
if (opt$method == "singler") {
  cat("Reference dataset:", opt$reference, "\n")
} else if (opt$method == "azimuth") {
  cat("Azimuth reference:", opt$azimuth_reference, "\n")
}
cat("\nData summary:\n")
cat("  Total cells:", ncol(seurat_obj), "\n")
cat("  Total features:", nrow(seurat_obj), "\n")
cat("  Number of samples:", length(unique(seurat_obj$sample)), "\n")
cat("  Number of clusters:", length(unique(seurat_obj$seurat_clusters)), "\n")

if (opt$method == "singler") {
  cat("  Number of cell types (SingleR main):", 
      length(unique(seurat_obj$singler_labels)), "\n")
  cat("\nCell type distribution:\n")
  print(table(seurat_obj$singler_labels))
  
  if ("singler_labels_fine" %in% colnames(seurat_obj@meta.data)) {
    cat("\n  Number of cell types (SingleR fine):", 
        length(unique(seurat_obj$singler_labels_fine)), "\n")
  }
  
} else if (opt$method == "azimuth") {
  if ("predicted.celltype.l1" %in% colnames(seurat_obj@meta.data)) {
    cat("  Number of cell types (Azimuth L1):", 
        length(unique(seurat_obj$predicted.celltype.l1)), "\n")
    cat("\nCell type distribution (Level 1):\n")
    print(table(seurat_obj$predicted.celltype.l1))
  }
  
  if ("predicted.celltype.l2" %in% colnames(seurat_obj@meta.data)) {
    cat("\n  Number of cell types (Azimuth L2):", 
        length(unique(seurat_obj$predicted.celltype.l2)), "\n")
    cat("\nCell type distribution (Level 2):\n")
    print(table(seurat_obj$predicted.celltype.l2))
  }
}

cat("\nCells per sample:\n")
print(table(seurat_obj$sample))

cat("\nCells per cluster:\n")
print(table(seurat_obj$seurat_clusters))

cat("\nOutput files:\n")
cat("  Annotated RDS object:", opt$output, "\n")
cat("  Marker genes: cell_type_markers.csv\n")
cat("  Plots: annotation_plots/\n")
cat(rep("=", 70), "\n", sep = "")
sink()

cat("\nCell type annotation complete!\n")
cat("Annotated object saved to:", opt$output, "\n")
cat("Plots saved to: annotation_plots/\n")
cat("Marker genes saved to: cell_type_markers.csv\n")
cat("Summary saved to: annotation_summary.txt\n")