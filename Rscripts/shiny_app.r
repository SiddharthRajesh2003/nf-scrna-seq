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
})


# ============================================================================
# UI Definition
# ============================================================================

ui <- fluidPage(
  theme = shinytheme("flatly"),

  titlePanel("Single-cell RNA-seq Data Explorer"),
  
  sidebarLayout(
    sidebarPanel (
      width = 3,

      # File Input
      fileInput("seuratFile", "Upload Seurat Object (.rds)",
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
                    selected = "umap",
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
                    choices = NULL,
                    options = list(maxItems = 6)),
        
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
        downloadButtuon("download_data", "Download Metadata")
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
                        withSpinner(plotlyOutput("overview_celltype", height = "500px")))
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
      
      tabPanel("Dim Reduction",
                fluidRow(12,
                          h4("Interactive Dimnesionality Reduction Plot"),
                          withSpinner(plotlyOutput("dim_reduction_plot", height = '700px'))
                        ),
                hr(),

                fluidRow(
                  column(12,
                          conditionalPanel(
                              condition = 'input.split_by_sample',
                              withSpinner(plotlyOutput("dim_reduction_split", height = '600px'))
                          ))
                )
      ),

      tabPanel("Gene Expression",
                fluidRow(
                    column(12,
                      h4("Feature Expression Plot"),
                      withSpinner(plotlyOutput("feature_plot", height = '700px'))
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
                        withSpinner(plotlyOutput("marker_heatmap", height = "600px"))
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

# ============================================================================
# Server Logic
# ============================================================================

server <- function(input, output, session) {
  
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
                        selected = if(length(celltype_cols) > 0) celltype_cols[1]
                                    else "seurat_clusters"
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
      "Features: ", nrows(obj), "\n",
      "Samples: ", length(unique(obj$sample)), "\n",
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

    p <- DimPlot(obj, reduction = 'umap', group.by = 'sample') +
      theme_minimal() +
      theme(legend.position = 'right')

    ggplotly(p, tooltip = c("Group"))
  })

  output$celltype_barplot <- renderPlot({
    req(seurat_obj())
    obj <- seurat_obj()
  })
}