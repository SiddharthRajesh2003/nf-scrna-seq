#!/usr/bin/env Rscript

# Single-cell RNA-seq Visualization Shiny App
# For visualizing integrated and annotated Seurat Objects

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(DT)
  library(plotly)
  library(shinythemes)
  library(shinycssloaders)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})


# ============================================================================
# UI Definition
# ============================================================================

ui <- fluidPage(
  theme = shinytheme("flatly"),

  titlePanel("Single-cell RNA-seq Data Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,

      # File Input
      fileInput("seurat_file", "Upload Seurat Object (.rds)",
                accept = c(".rds", ".RDS")),
      
      hr(),

      # Conditional Panels based on data loading
      conditionalPanel(
        condition = "output.data_loaded",

        h4("Dataset Information"),
        verbatimTextOutput("data_summary"),

        hr(),

        h4("Visualization Options"),

        selectInput("reduction",
                    "Reduction Method:",
                    choices = c("UMAP" = "umap",
                                "TSNE" = "tsne",
                                "PCA" = "pca"),
                    selected = "umap"
        ),

        selectInput("color_by",
                    "Color By:",
                    choices = NULL
        ),

        checkboxInput("show_labels", "Show Labels", value = TRUE),
        checkboxInput("split_by_sample", "Split by Sample", value = FALSE),

        sliderInput("point_size",
                    "Point Size:",
                    min = 0.5,
                    max = 5,
                    value = 1.5,
                    step = 0.1
        ),

        hr(),

        # Feature Plot Options
        h4("Feature Plot Options"),
        selectizeInput('features',
                      'Select Gene(s)',
                      choices = NULL,
                      multiple = TRUE,
                      options = list(maxItems = 6)),

        hr(),

        # Marker Genes
        h4("Cluster markers"),
        selectInput("marker_cluster",
                    "Select Cluster:",
                    choices = NULL),
        
        numericInput("top_n_markers",
                    "Top N Markers:",
                    value = 10,
                    min = 5,
                    max = 50,
                    step = 5),
          
        actionButton("find_markers", "Find Markers",
                      class = "btn-primary",
                      icon = icon("search")
                    ),
        
        hr(),
        
        h4("Download"),
        downloadButton("download_plot", "Download Current Plot"),
        downloadButton("download_data", "Download Metadata")
      )
    ),
    mainPanel(
      width = 9,

      tabsetPanel(
        id = 'main_tabs',

        # Tab 1: Overview
        tabPanel("Overview",
          fluidRow(
            column(6,
              h4("UMAP by Cell Type"),
              withSpinner(plotlyOutput("overview_celltype", height = "500px"))
            ),
            column(6,
              h4("UMAP by Sample"),
              withSpinner(plotlyOutput("overview_sample", height = "500px"))
            )
          ),
          hr(),
          fluidRow(
            column(6,
              h4("Cell Type Distribution"),
              withSpinner(plotlyOutput("celltype_barplot", height = "400px"))
            ),
            column(6,
              h4("QC Metrics"),
              withSpinner(plotlyOutput("qc_violin", height = '400px'))
            )
          )
        ),
        
        # Tab 2: Dim Reduction
        tabPanel("Dim Reduction",
          fluidRow(
            column(12,
              h4("Interactive Dimensionality Reduction Plot"),
              withSpinner(plotlyOutput("dim_reduction_plot", height = '700px'))
            )
          ),
          hr(),
          fluidRow(
            column(12,
              conditionalPanel(
                condition = 'input.split_by_sample',
                withSpinner(plotlyOutput("dim_reduction_split", height = '600px'))
              )
            )
          )
        ),

        # Tab 3: Gene Expression
        tabPanel("Gene Expression",
          fluidRow(
            column(12,
              h4("Feature Expression Plot"),
              withSpinner(plotOutput("feature_plot", height = '700px'))
            )
          ),
          hr(),
          fluidRow(
            column(6,
              h4("Violin Plot"),
              withSpinner(plotlyOutput("feature_violin", height = '500px'))
            ),
            column(6,
              h4("Dot Plot"),
              withSpinner(plotlyOutput("feature_dot", height = '500px'))
            )
          )
        ),

        # Tab 4: Cluster Markers
        tabPanel("Cluster Markers",
          fluidRow(
            column(12,
              h4("Top Marker Genes"),
              withSpinner(DTOutput("marker_table"))
            )
          ),
          hr(),
          fluidRow(
            column(12,
              h4("Marker Gene Heatmap"),
              withSpinner(plotOutput("marker_heatmap", height = "600px"))
            )
          ),
          hr(),
          fluidRow(
            column(12,
              h4("Marker Gene Expression"),
              withSpinner(plotlyOutput("marker_feature_plot", height = '600px'))
            )
          )
        ),

        # Tab 5: Cell-type Composition
        tabPanel("Cell-type Composition",
          fluidRow(
            column(6,
              h4("Cell Type by Sample"),
              withSpinner(plotlyOutput("composition_stacked", height = '500px'))
            ),
            column(6,
              h4("Sample by Cell type"),
              withSpinner(plotlyOutput("composition_stacked_flip", height = '500px'))
            )
          ),
          hr(),
          fluidRow(
            column(12,
              h4("Composition table"),
              withSpinner(DTOutput("composition_table"))
            )
          )
        ),

        # Tab 6: Metadata
        tabPanel("Metadata",
          fluidRow(
            column(12,
              h4("Cell Metadata"),
              withSpinner(DTOutput("metadata_table"))
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# Server Logic
# ============================================================================

server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 5000*1024^2)
  # Reactive values
  seurat_obj <- reactiveVal(NULL)
  current_markers <- reactiveVal(NULL)

  # Load Seurat Object
  observeEvent(input$seurat_file, {
    req(input$seurat_file)

    withProgress(message = 'Loading Seurat object...', value = 0, {
      incProgress(0.3, detail = 'Reading RDS file')
      obj <- readRDS(input$seurat_file$datapath)

      incProgress(0.7, detail = 'Processing metadata')
      seurat_obj(obj)

      # Update UI elements
      incProgress(0.9, detail = 'Updating UI')

      # Get available metdata columns
      meta_cols <- colnames(obj@meta.data)

      # Prioritize cell type columns
      celltype_cols <- grep("label|type|celltype|cluster", meta_cols,
                            value = TRUE, ignore.case = TRUE)
      
      updateSelectInput(session, "color_by",
                        choices = meta_cols,
                        selected = (if(length(celltype_cols) > 0) celltype_cols[1] else "seurat_clusters")
                      )
      
      genes <- rownames(obj)
      updateSelectizeInput(session, "features",
                            choices = genes,
                            server = TRUE
                          )
      
      # Update cluster selection
      clusters <- sort(unique(obj$seurat_clusters))
      updateSelectInput(session, "marker_cluster",
                        choices = clusters
                        )
      
      incProgress(1, detail = 'Done!')
    })
  })

  # Data loaded Flag
  output$data_loaded <- reactive({
    !is.null(seurat_obj())
  })

  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

  # Data summary
  output$data_summary <- renderText({
    req(seurat_obj())

    obj <- seurat_obj()

    paste0(
      "Cells: ", ncol(obj), "\n",
      "Features: ", nrow(obj), "\n",
      "Samples: ", length(unique(obj$Sample_ID)), "\n",
      "Clusters: ", length(unique(obj$seurat_clusters))
    )
  })

  # ============================================================================
  # Overview Tab
  # ============================================================================
  output$overview_celltype <- renderPlotly({
    req(seurat_obj())
    obj <- seurat_obj()

    # Get cell type column (prioritize annotation columns)
    celltype_col <- grep("label|type|celltype",
                          colnames(obj@meta.data),
                          value = TRUE, ignore.case = TRUE) [1]
    if(is.na(celltype_col)) celltype_col <- "seurat_clusters"

    p <- DimPlot(obj, reduction = "umap", group.by = celltype_col,
                  label = TRUE, repel = TRUE) +
          theme_minimal() +
          theme(legend.position = 'right')
    
        ggplotly(p, tooltip = c("group"))
  })

  output$overview_sample <- renderPlotly({
    req(seurat_obj())
    obj <- seurat_obj()

    p <- DimPlot(obj, reduction = 'umap', group.by = 'Sample_ID') +
      theme_minimal() +
      theme(legend.position = 'right')

    ggplotly(p, tooltip = c("Group"))
  })

  output$celltype_barplot <- renderPlotly({
    req(seurat_obj())
    obj <- seurat_obj()

    celltype_col <- input$color_by

    if (is.null(celltype_col)) celltype_col <- 'seurat_clusters'

    df <- data.frame(
      celltype = obj@meta.data[[celltype_col]],
      sample = obj$Sample_ID
    )

    ggplot(df, aes(x = celltype, fill = sample)) +
      geom_bar(position = 'dodge') +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = 'Cell Type', y = 'Number of Cells', fill = 'Sample_ID')
  })

  output$qc_violin <- renderPlotly({
    req(seurat_obj())
    obj <- seurat_obj()

    VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
              ncol=3, pt.size = 0, group.by = 'Sample_ID') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  })
  
  # ============================================================================
  # Dimensionality Reduction Tab
  # ============================================================================
  
  output$dim_reduction_plot <- renderPlotly({
      req(seurat_obj(), input$color_by)
    
    obj <- seurat_obj()

    p <- DimPlot(
        obj, reduction = input$reduction,
        group.by = input$color_by, label = input$show_labels,
        repel = TRUE, pt.size = input$point_size
    ) +
      theme_minimal() +
      theme(legend.position = 'right')

    ggplotly(p, tooltip = c("group"))
  })

  output$dim_reduction_split <- renderPlotly({
    req(seurat_obj(), input$split_by_sample)
    if(!input$split_by_sample) return(NULL)
    
    obj <- seurat_obj()

    p <- DimPlot(
      obj, reduction = input$reduction,
      group.by = input$color_by,
      split.by = 'Sample_ID',
      label = input$show_labels,
      repel = TRUE, pt.size = input$point_size,
      ncol = 2
    ) +
      theme_minimal()
    
    ggplotly(p, tooltip = c("group"))
  })
  
  # ============================================================================
  # Feature Expression Tab
  # ============================================================================
  
  output$feature_plot <- renderPlot({
    req(seurat_obj())
    req(length(input$features) > 0)
    
    obj <- seurat_obj()
    
    p <- FeaturePlot(obj,
                features = input$features,
                reduction = input$reduction,
                ncol = min(3, length(input$features)),
                pt.size = input$point_size) +
      scale_color_viridis_c() +
      theme_minimal()
  }, height = function() {
    req(input$features)
    n_features <- length(input$features)
    if(n_features == 0) return(400)
    ceiling(n_features / 3) * 250
  })

  output$feature_violin <- renderPlotly({
    req(seurat_obj(), input$features)
    if(length(input$features) == 0) return(NULL)
    
    obj <- seurat_obj()

    p <- VlnPlot(
      obj, features = input$features,
      group.by = input$color_by, ncol = min(3, length(input$features)),
      pt.size = 0
    ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggplotly(p, tooltip = c('group'))
  })

  output$feature_dot <- renderPlotly({
    req(seurat_obj(), input$features)
    if(length(input$features) == 0) return(NULL)
    
    obj <- seurat_obj()

    p <- DotPlot(
      obj, features = input$features,
      group.by = input$color_by
    ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggplotly(p, tooltip = c('group'))
  })
  
  # ============================================================================
  # Cluster Markers Tab
  # ============================================================================
  
  observeEvent(input$find_markers, {
    req(seurat_obj(), input$marker_cluster)

    withProgress(message = 'Finding marker genes...', value = 0, {
      obj <- seurat_obj()
      cluster <- input$marker_cluster

      incProgress(0.5, detail = paste("Analyzing cluster", cluster))

      # Find markers for selected cluster
      markers <- FindMarkers(
        obj, ident.1 = cluster,
        only.pos = TRUE, min.pct = 0.25,
        logfc.threshold = 0.25
      )

      markers$gene <- rownames(markers)
      markers$cluster <- cluster

      current_markers(markers)

      incProgress(1, detail = 'Done!')
    })
  })

  output$marker_table <- renderDT({
    req(current_markers())

    markers <- current_markers()

    top_markers <- markers %>%
      arrange(desc(avg_log2FC)) %>%
      head(input$top_n_markers) %>%
      select(gene, avg_log2FC, pct.1, pct.2, p_val_adj)

    datatable(top_markers,
              options = list(pageLength = 10, scrollX = TRUE),
              rownames = FALSE
    ) %>%
      formatRound(columns = c("avg_log2FC", "pct.1", "pct.2"), digits = 3) %>%
      formatSignif(columns = "p_val_adj", digits = 3)
  })

  output$marker_heatmap <- renderPlot({
    req(seurat_obj(), current_markers())
    obj <- seurat_obj()
    markers <- current_markers()

    top_genes <- markers %>%
      arrange(desc(avg_log2FC)) %>%
      head(min(20, input$top_n_markers)) %>%
      pull(gene)

    # Subset and get expression data
    obj_subset <- subset(obj, downsample = 100)
    expr_data <- GetAssayData(obj_subset, assay = 'RNA', layer = 'data')
    expr_matrix <- as.matrix(expr_data[top_genes, ])

    # Scale the data
    expr_matrix <- t(scale(t(expr_matrix)))

    # Get metadata for annotation
    cell_groups <- obj_subset@meta.data[[input$color_by]]

    # Create color mapping
    n_groups <- length(unique(cell_groups))
    group_colors <- setNames(
      scales::hue_pal()(n_groups),
      unique(cell_groups)
    )

    # Create column annotation
    col_anno <- HeatmapAnnotation(
      CellType = cell_groups,
      col = list(CellType = group_colors),
      show_legend = TRUE,
      annotation_name_side = 'left'
    )

    # Create Heatmap
    ComplexHeatmap::Heatmap(
      expr_matrix,
      name = "Expression",
      col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
      top_annotation = col_anno,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_column_names = FALSE,
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 8),
      column_title = paste("Top ", length(top_genes), " Marker genes"),
      heatmap_legend_param = list(
        title = 'Scaled\nExpression',
        direction = 'vertical'
      )
    )
  })

  output$marker_feature_plot <- renderPlotly({
    req(seurat_obj(), current_markers())

    obj <- seurat_obj()
    markers <- current_markers()

    top_genes <- markers %>%
      arrange(desc(avg_log2FC)) %>%
      head(min(6, input$top_n_markers)) %>%
      pull(gene)

    p <- FeaturePlot(
      obj, features = top_genes,
      reduction = input$reduction,
      ncol = 3
    ) +
      scale_color_viridis_c() +
      theme_minimal()

    ggplotly(p, tooltip = c('group'))
  })

  # ============================================================================
  # Composition Tab
  # ============================================================================
  
  output$composition_stacked <- renderPlotly({
    req(seurat_obj())
    obj <- seurat_obj()

    df <- data.frame(
      celltype = obj@meta.data[[input$color_by]],
      sample = obj$Sample_ID
    )

    p <- ggplot(df, aes(x = sample, fill = celltype)) +
      geom_bar(position = 'fill') +
      scale_y_continuous(labels = scales::percent) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = 'Sample_ID', y = 'Proportion', fill = 'Cell Type')

    ggplotly(p, tooltip = c('group'))
  })

  output$composition_stacked_flip <- renderPlotly({
    req(seurat_obj())
    obj <- seurat_obj()

    df <- data.frame(
      celltype = obj@meta.data[[input$color_by]],
      sample = obj$Sample_ID
    )

    p <- ggplot(df, aes(x = celltype, fill = sample)) +
      geom_bar(position = 'fill') +
      scale_y_continuous(labels = scales::percent) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = 'Cell Type', y = 'Proportion', fill = 'Sample_ID')

    ggplotly(p, tooltip = c("group"))
  })

  output$composition_table <- renderDT({
    req(seurat_obj())

    obj <- seurat_obj()

    df <- data.frame(
      celltype = obj@meta.data[[input$color_by]],
      sample = obj$Sample_ID
    )

    comp_table <- df %>%
      group_by(sample, celltype) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(sample) %>%
      mutate(percentage = round(count / sum(count) * 100, 2)) %>%
      arrange(sample, desc(count))

    datatable(
      comp_table,
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    )
  })

  # ============================================================================
  # Metadata Tab
  # ============================================================================
  

  output$metadata_table <- renderDT({
    req(seurat_obj())

    obj <- seurat_obj()

    datatable(obj@meta.data,
              options = list(pageLength = 20, scrollX = TRUE),
              filter = 'top'
    )
  })

  # ============================================================================
  # Download Handlers
  # ============================================================================
  
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      req(seurat_obj())
      obj <- seurat_obj()
      
      p <- DimPlot(obj,
                  reduction = input$reduction,
                  group.by = input$color_by,
                  label = input$show_labels,
                  repel = TRUE,
                  pt.size = input$point_size) +
        theme_minimal()
      
      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("metadata_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(seurat_obj())
      obj <- seurat_obj()
      write.csv(obj@meta.data, file, row.names = TRUE)
    }
  )
}

# ============================================================================
# Run App
# ============================================================================

shinyApp(ui = ui, server = server)