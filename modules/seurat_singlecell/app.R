# ---
# Transcriptomic Explorer (v2.0)
#
# Description:
#   Interactive R/Shiny dashboard for exploring Single-Cell RNA-seq data.
#   Features:
#     - Interactive UMAP/Feature Plots with Split/Side-by-Side views.
#     - Advanced Density & Glow visualization layers.
#     - Real-time subsetting and reprojection (with Seurat v5 compatibility).
#     - Automated Differential Expression (DEG) and Heatmap generation.
#     - High-quality vector (PDF) downloads.
#
# Usage:
#   1. Ensure all packages listed below are installed.
#   2. Run this script in RStudio: shiny::runApp()
# ---

# 1. Automatic Package Installation and Loading
# --------------------------------------------------------------------------
cat("--- [Transcriptomic Explorer] Initializing Environment ---\n")
required_packages <- c(
  "shiny", "Seurat", "ggplot2", "dplyr", "DT", "shinyjqui", "patchwork",
  "ggrepel", "RColorBrewer", "plotly", "shinycssloaders", "shinyjs",
  "scales", "colourpicker", "ggnewscale", "tibble", "viridis",
  "pheatmap", "harmony", "ggrastr", "MASS", "stringr"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
lapply(required_packages, library, character.only = TRUE)
cat("--- [Transcriptomic Explorer] Dependencies Loaded ---\n\n")

# --- Helper function for download UI ---
download_plot_ui <- function(id_prefix, width_val = 8, height_val = 6) {
  ns <- NS(id_prefix)
  tagList(
    actionLink(ns("toggle_download_link"), label = h5("Download Options ⯆"), style = "text-decoration: none; color: #337ab7; font-weight: bold;"),
    shinyjs::hidden(
      div(id = ns("download_content"),
          wellPanel(style="background-color: #f8f9fa; border-left: 3px solid #337ab7;",
                    fluidRow(
                      column(3, textInput(ns("filename"), "Filename (.pdf)", value = paste0(id_prefix, "_plot.pdf"))),
                      column(2, numericInput(ns("pdf_width"), "Width (in)", value = width_val, min = 4, max = 30, step = 0.5)),
                      column(2, numericInput(ns("pdf_height"), "Height (in)", value = height_val, min = 4, max = 30, step = 0.5)),
                      column(3, checkboxInput(ns("rasterize"), "Rasterize points", value = TRUE),
                             helpText("Recommended for large datasets.", style="font-size: smaller; color: grey;")),
                      column(2, style = "margin-top: 25px;", downloadButton(ns("download_button"), "Download PDF", class = "btn-primary"))
                    )
          )
      )
    )
  )
}

# 2. Define the User Interface (UI)
# --------------------------------------------------------------------------
ui <- fluidPage(
  shinyjs::useShinyjs(),
  theme = shinythemes::shinytheme("flatly"), # Clean modern theme
  
  # Custom CSS for better anchoring and spacing
  tags$head(
    tags$style(HTML("
      .well { background-color: #ffffff; border: 1px solid #e3e3e3; box-shadow: 0 1px 1px rgba(0,0,0,.05); }
      h4 { color: #2c3e50; font-weight: 600; margin-top: 20px; border-bottom: 1px solid #eee; padding-bottom: 5px; }
      h5 { color: #7f8c8d; font-weight: 600; }
      .action-link-anchor { cursor: pointer; color: #2980b9; font-weight: bold; }
      .sidebar-panel { background-color: #f8f9fa; }
    "))
  ),
  
  titlePanel("Transcriptomic Explorer (v2.0)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. Load Data"),
      fileInput("seurat_object_file", "Upload Seurat Object (.rds)", accept = ".rds"),
      hr(),
      uiOutput("main_controls_ui")
    ),
    
    mainPanel(
      width = 9,
      
      # --- GLOBAL SETTINGS BAR ---
      fluidRow(
        wellPanel(
          fluidRow(
            column(3,
                   h5("Global Appearance"),
                   sliderInput("pt_size", "Point Size:", min = 0.0, max = 5, value = 1.0, step = 0.1),
                   sliderInput("alpha", "Transparency:", min = 0.0, max = 1, value = 0.8, step = 0.1)
            ),
            column(3,
                   h5("Clustering"),
                   sliderInput("resolution", "Leiden Resolution:", min = 0.1, max = 2.0, value = 0.5, step = 0.1),
                   actionButton("run_clustering", "Apply Resolution", class = "btn-info btn-sm")
            ),
            column(3,
                   uiOutput("label_overlay_ui_placeholder")
            ),
            column(3,
                   actionLink("toggle_colors_link", label = "Custom Colors ⯆", class = "action-link-anchor"),
                   shinyjs::hidden(
                     div(id = "colors_content",
                         uiOutput("color_category_ui"),
                         uiOutput("color_picker_ui")
                     )
                   )
            )
          ),
          
          hr(style="margin-top: 10px; margin-bottom: 10px; border-color: #ddd;"),
          
          # --- DENSITY / GLOW SETTINGS ---
          fluidRow(
            column(12,
                   actionLink("toggle_density_link", label = "Density & Glow Overlay ⯆", class = "action-link-anchor"),
                   shinyjs::hidden(
                     div(id = "density_content",
                         wellPanel(style = "background-color: #fcfcfc;",
                                   fluidRow(
                                     column(3,
                                            h5("General"),
                                            checkboxInput("show_density", "Enable Density Layers", value = FALSE),
                                            conditionalPanel(
                                              condition = "input.show_density == true",
                                              checkboxInput("greyscale_density_mode", "Greyscale Mode", value = FALSE),
                                              conditionalPanel(
                                                condition = "input.greyscale_density_mode == true",
                                                checkboxInput("color_density_outlines", "Color Outlines by Group", value = FALSE)
                                              ),
                                              hr(style="margin: 5px 0;"),
                                              checkboxInput("normalize_density", "Normalize by Group", value = FALSE),
                                              sliderInput("padding_factor", "Padding (%):", min = 5, max = 100, value = 10, step = 1),
                                              sliderInput("smoothing", "Smoothing:", min = 0.1, max = 3.0, value = 0.7, step = 0.05)
                                            )
                                     ),
                                     column(3,
                                            h5("Peripheral Glow"),
                                            conditionalPanel(
                                              condition = "input.show_density == true",
                                              checkboxInput("show_glow", "Show Glow", value = TRUE),
                                              sliderInput("glow_opacity", "Opacity:", min = 0, max = 1, value = 0.14, step = 0.01),
                                              sliderInput("glow_extent", "Extent (Rings):", min = 1, max = 20, value = 20, step = 1),
                                              sliderInput("glow_outline_width", "Outline Width:", min = 0, max = 2, value = 0, step = 0.1)
                                            )
                                     ),
                                     column(3,
                                            h5("Main Core"),
                                            conditionalPanel(
                                              condition = "input.show_density == true",
                                              checkboxInput("show_core", "Show Core", value = TRUE),
                                              sliderInput("core_opacity", "Opacity:", min = 0, max = 1, value = 0.1, step = 0.01),
                                              sliderInput("core_size", "Size (Rings):", min = 1, max = 20, value = 18, step = 1),
                                              sliderInput("darken_factor", "Darkening:", min = 0, max = 1, value = 0.92, step = 0.01),
                                              sliderInput("core_outline_width", "Outline Width:", min = 0, max = 2, value = 0, step = 0.1)
                                            )
                                     ),
                                     column(3,
                                            h5("Dense Core"),
                                            conditionalPanel(
                                              condition = "input.show_density == true",
                                              checkboxInput("show_dense_core", "Show Dense Core", value = TRUE),
                                              sliderInput("dense_core_opacity", "Opacity:", min = 0, max = 1, value = 0.24, step = 0.01),
                                              sliderInput("dense_core_size", "Size (Rings):", min = 1, max = 10, value = 5, step = 1),
                                              sliderInput("dense_core_darkening", "Darkening:", min = 0, max = 1, value = 0.65, step = 0.01),
                                              sliderInput("dense_core_outline_width", "Outline Width:", min = 0, max = 2, value = 0, step = 0.1)
                                            )
                                     )
                                   ),
                                   # Density Colors Sub-Toggle
                                   fluidRow(
                                     column(12, 
                                            actionLink("toggle_density_color_link", label = "Density Custom Colors ⯆", style = "font-size: 0.9em; color: #555;"),
                                            shinyjs::hidden(
                                              div(id = "density_color_content",
                                                  hr(),
                                                  fluidRow(
                                                    column(4, colourpicker::colourInput("glow_outline_color", "Glow Outline", value = "#FFFFFF"), checkboxInput("glow_fill_use_group_color", "Fill w/ Group Color", value = TRUE)),
                                                    column(4, colourpicker::colourInput("core_outline_color", "Core Outline", value = "#FFFFFF"), checkboxInput("core_fill_use_group_color", "Fill w/ Group Color", value = TRUE)),
                                                    column(4, colourpicker::colourInput("dense_core_outline_color", "Dense Outline", value = "#FFFFFF"), checkboxInput("dense_core_fill_use_group_color", "Fill w/ Group Color", value = TRUE))
                                                  )
                                              )
                                            )
                                     )
                                   )
                         )
                     )
                   )
            )
          )
        )
      ),
      
      # --- MAIN TABSET ---
      tabsetPanel(
        type = "tabs",
        id = "main_tabs",
        
        # TAB 1: UMAP
        tabPanel("Explorer (UMAP)",
                 br(),
                 fluidRow(
                   column(12,
                          h4("Dimensionality Reduction"),
                          jqui_resizable(plotOutput("umap_plot", height = "600px")),
                          download_plot_ui("umap", width_val = 10, height_val = 8),
                          hr(),
                          h4("Feature Expression"),
                          jqui_resizable(plotOutput("feature_plot_grid", height = "600px")),
                          download_plot_ui("feature", width_val = 10, height_val = 8)
                   )
                 )
        ),
        
        # TAB 2: DOT PLOT
        tabPanel("Dot Plot Analysis",
                 br(),
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     h4("Settings"),
                     uiOutput("dot_plot_controls_ui")
                   ),
                   mainPanel(
                     width = 9,
                     shinycssloaders::withSpinner(
                       jqui_resizable(plotOutput("dot_plot", height="800px"))
                     ),
                     download_plot_ui("dotplot", width_val = 10, height_val = 8)
                   )
                 )
        ),
        
        # TAB 3: DISTRIBUTION
        tabPanel("Distribution Analysis",
                 br(),
                 sidebarLayout(
                   sidebarPanel(
                     width = 3,
                     h4("Settings"),
                     uiOutput("dist_analysis_controls_ui"),
                     hr(),
                     h5("Download"),
                     textInput("dist_analysis_filename", "Filename", value = "Distribution.pdf"),
                     numericInput("dist_analysis_pdf_width", "Width", value = 10),
                     numericInput("dist_analysis_pdf_height", "Height", value = 8),
                     downloadButton("download_dist_analysis_plot", "Download Plot")
                   ),
                   mainPanel(
                     width = 9,
                     tabsetPanel(
                       id = "distribution_plot_tabs",
                       tabPanel("Bar Plot",
                                shinycssloaders::withSpinner(
                                  jqui_resizable(plotOutput("dist_analysis_bar_plot", height="800px"))
                                )
                       ),
                       tabPanel("Donut Plot",
                                shinycssloaders::withSpinner(
                                  jqui_resizable(plotOutput("dist_analysis_donut_plot", height="800px"))
                                )
                       )
                     )
                   )
                 )
        ),
        
        # TAB 4: DEG HEATMAP
        tabPanel("DEG Heatmap",
                 br(),
                 sidebarLayout(
                   sidebarPanel(
                     width = 4,
                     h4("Differential Expression"),
                     selectInput("deg_group_by", "Group Cells By:", choices = NULL),
                     radioButtons("deg_mode", "Mode:",
                                  choices = c("Pairwise (Group 1 vs 2)" = "two_groups",
                                              "One-vs-Rest (All Markers)" = "one_vs_all"),
                                  selected = "two_groups"),
                     hr(),
                     conditionalPanel(
                       condition = "input.deg_mode == 'two_groups'",
                       uiOutput("deg_group1_ui"),
                       uiOutput("deg_group2_ui")
                     ),
                     conditionalPanel(
                       condition = "input.deg_mode == 'one_vs_all'",
                       uiOutput("deg_one_vs_all_groups_ui")
                     ),
                     hr(),
                     numericInput("deg_n_genes", "Top N Genes:", value = 10, min = 1, max = 50),
                     actionButton("run_deg_heatmap", "Generate Heatmap", icon = icon("fire"), class = "btn-success")
                   ),
                   mainPanel(
                     width = 8,
                     shinycssloaders::withSpinner(
                       jqui_resizable(plotOutput("deg_heatmap_plot", height = "800px"))
                     ),
                     download_plot_ui("heatmap", width_val = 8, height_val = 10)
                   )
                 )
        )
      )
    )
  )
)


# 3. Define the Server Logic
# --------------------------------------------------------------------------
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 3000 * 1024^2) # Allow up to 3GB uploads
  
  # --- DENSITY HELPER FUNCTIONS ---
  darken_color <- function(color_hex, factor=0.5) {
    tryCatch({
      rgb_vals <- col2rgb(color_hex)
      darker_rgb <- rgb_vals * factor
      rgb(darker_rgb[1,], darker_rgb[2,], darker_rgb[3,], maxColorValue=255)
    }, error = function(e) color_hex)
  }
  
  calculate_density_polygons <- function(data, x_col, y_col, group_col, smoothing, limits, normalize = FALSE, nlevels = 20) {
    all_groups <- na.omit(unique(as.character(data[[group_col]])))
    polygon_data <- lapply(all_groups, function(current_group) {
      df_subset <- data %>% filter(.data[[group_col]] == current_group)
      if (nrow(df_subset) < 5) return(NULL)
      bw_x <- smoothing * IQR(df_subset[[x_col]], na.rm=TRUE)
      bw_y <- smoothing * IQR(df_subset[[y_col]], na.rm=TRUE)
      if (is.na(bw_x) || bw_x <= 0.01) bw_x <- 0.1
      if (is.na(bw_y) || bw_y <= 0.01) bw_y <- 0.1
      kde_result <- tryCatch({
        MASS::kde2d(df_subset[[x_col]], df_subset[[y_col]], n=100, h=c(bw_x, bw_y), lims=limits)
      }, error = function(e) NULL)
      if (is.null(kde_result)) return(NULL)
      if (isTRUE(normalize)) {
        sorted_z <- sort(kde_result$z, decreasing = TRUE)
        cumulative_z <- cumsum(sorted_z)
        total_sum <- sum(kde_result$z)
        normalized_cumulative_z <- cumulative_z / total_sum
        prob_levels <- seq(0.01, 0.99, length.out = nlevels)
        level_indices <- sapply(prob_levels, function(p) which.min(abs(normalized_cumulative_z - p)))
        contour_levels <- sorted_z[level_indices]
        contour_levels <- unique(sort(contour_levels))
        contourLines(kde_result$x, kde_result$y, kde_result$z, levels = contour_levels)
      } else {
        contourLines(kde_result$x, kde_result$y, kde_result$z, nlevels = nlevels)
      }
    })
    names(polygon_data) <- all_groups
    return(polygon_data)
  }
  
  # --- REACTIVE VALUES ---
  active_object <- reactiveVal(NULL)
  object_history <- reactiveVal(list())
  deg_heatmap_data <- reactiveVal(NULL)
  color_palettes <- reactiveVal(list())
  signature_scores <- reactiveVal(character(0))
  
  # --- HELPERS: METADATA & FEATURES ---
  categorical_metadata <- reactive({
    req(active_object())
    obj <- active_object()
    meta_cols <- colnames(obj@meta.data)
    potential_cat_cols <- sapply(meta_cols, function(col) {
      col_data <- obj@meta.data[[col]]
      is.factor(col_data) || is.character(col_data) || (is.numeric(col_data) && length(unique(col_data)) < 50)
    })
    choices <- meta_cols[potential_cat_cols]
    names(choices) <- gsub("_", " ", tools::toTitleCase(choices))
    if ("seurat_clusters" %in% meta_cols) {
      choices <- choices[choices != "seurat_clusters"]
      choices <- c("Clusters" = "seurat_clusters", choices)
    }
    return(choices)
  })
  
  numerical_features_for_subset <- reactive({
    req(active_object())
    obj <- active_object()
    genes_in_use <- sapply(1:input$n_feature_plots, function(i) input[[paste0("feature_gene_", i)]])
    valid_genes <- genes_in_use[!is.null(genes_in_use) & genes_in_use != "" & genes_in_use %in% rownames(obj)]
    valid_scores <- signature_scores()[signature_scores() %in% colnames(obj@meta.data)]
    choices <- c(valid_genes, valid_scores)
    names(choices) <- c(
      if(length(valid_genes) > 0) paste0("Gene: ", valid_genes) else NULL,
      if(length(valid_scores) > 0) paste0("Score: ", valid_scores) else NULL
    )
    return(choices)
  })
  
  active_reduction_dims <- reactive({
    obj <- active_object(); req(obj)
    if ("harmony" %in% names(obj@reductions)) { "harmony" } else { "pca" }
  })
  
  active_reduction_umap <- reactive({ "umap" })
  
  # --- DATA LOADING ---
  observeEvent(input$seurat_object_file, {
    req(input$seurat_object_file)
    withProgress(message = 'Loading Seurat Object...', value = 0, {
      tryCatch({
        obj <- readRDS(input$seurat_object_file$datapath)
        if (!"pca" %in% names(obj@reductions)) {
          showNotification("Error: No 'pca' reduction found.", type = "error"); return()
        }
        reduction_for_neighbors <- if ("harmony" %in% names(obj@reductions)) "harmony" else "pca"
        dims_to_use <- 1:min(20, ncol(obj@reductions[[reduction_for_neighbors]]))
        obj <- FindNeighbors(obj, dims = dims_to_use, reduction = reduction_for_neighbors, verbose = FALSE)
        if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
          obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
        }
        active_object(obj)
        object_history(list(obj))
        showNotification("Object loaded successfully.", type = "message")
      }, error = function(e) { showNotification(paste("Error:", e$message), type = "error") })
    })
  })
  
  # --- DYNAMIC UI GENERATION (CONTROLS) ---
  output$main_controls_ui <- renderUI({
    req(active_object())
    
    all_palettes <- rownames(RColorBrewer::brewer.pal.info)
    color_choices <- c("Default (Grey->Blue)", "Viridis", all_palettes)
    
    tagList(
      h4("2. UMAP Visualization"),
      selectInput("color_by", "Color UMAP By:", choices = categorical_metadata()),
      uiOutput("subgroup_selection_ui"),
      hr(),
      
      h4("Advanced Visualization (Sided UMAP)"),
      radioButtons("split_mode", "Display Mode:",
                   choices = c("Combined" = "none", "Side-by-Side (Split View)" = "split", "Single Facet (Isolate)" = "isolate"),
                   selected = "none"),
      checkboxInput("apply_mode_to_features", "Apply Split to Feature Plots", value = FALSE),
      checkboxInput("apply_mode_to_dot_plots", "Apply Split to Dot Plots", value = FALSE),
      
      conditionalPanel(
        condition = "input.split_mode == 'split' || input.split_mode == 'isolate'",
        selectInput("split_by_var", "Variable to Split By:", choices = categorical_metadata())
      ),
      conditionalPanel(
        condition = "input.split_mode == 'split'",
        uiOutput("reorder_split_ui")
      ),
      conditionalPanel(
        condition = "input.split_mode == 'isolate'",
        uiOutput("isolate_values_ui")
      ),
      hr(),
      
      h4("3. Cell Subsetting & State"),
      uiOutput("subsetting_controls_ui"),
      uiOutput("object_state_controls_ui"),
      
      actionLink("toggle_umap_params_link", label = "3.2 UMAP Shape Parameters ⯆", class = "action-link-anchor"),
      shinyjs::hidden(
        div(id = "umap_params_content",
            wellPanel(
              sliderInput("umap_n_neighbors", "Number of Neighbors:", min = 5, max = 150, value = 30, step = 1),
              sliderInput("umap_min_dist", "Minimum Distance:", min = 0.0, max = 1.0, value = 0.3, step = 0.01),
              sliderInput("umap_spread", "Spread:", min = 0.1, max = 5.0, value = 1.0, step = 0.1),
              actionButton("rerun_umap", "Re-run UMAP", icon = icon("sync"))
            )
        )
      ),
      hr(),
      
      actionLink("toggle_annotation_link", label = h4("4. Cluster Annotation ⯆"), style = "text-decoration: none;"),
      shinyjs::hidden(
        div(id = "annotation_content",
            uiOutput("annotation_ui")
        )
      ),
      hr(),
      
      h4("5. Marker Analysis"),
      tabsetPanel(id = "marker_analysis_tabs",
                  tabPanel("Single Marker Plots",
                           numericInput("n_feature_plots", "Number of Feature Plots:", 1, 1, 12),
                           uiOutput("feature_plot_gene_inputs_ui"),
                           hr(),
                           h5("Plot Appearance"),
                           selectInput("feature_color_palette", "Color Palette:", choices = color_choices, selected = "Default (Grey->Blue)"),
                           checkboxInput("enable_custom_scale", "Enable Custom Scale", value = FALSE),
                           conditionalPanel(
                             condition = "input.enable_custom_scale == true",
                             sliderInput("feature_min_cutoff", "Min Cutoff (Quantile %)", min = 0, max = 49, value = 10, step = 1),
                             sliderInput("feature_max_cutoff", "Max Cutoff (Quantile %)", min = 50, max = 100, value = 90, step = 1)
                           )
                  ),
                  tabPanel("Composed Marker Scoring",
                           h5("Define Gene Signature"),
                           textAreaInput("positive_markers_csv", "Positive Markers (comma/space/newline separated)", rows=3, placeholder="Cd4, Cd8a, Il2"),
                           shinycssloaders::withSpinner(uiOutput("positive_marker_validation_ui"), type=6, size=0.5),
                           textAreaInput("negative_markers_csv", "Negative Markers (optional)", rows=3, placeholder="Foxp3, Tcf7"),
                           shinycssloaders::withSpinner(uiOutput("negative_marker_validation_ui"), type=6, size=0.5),
                           textInput("signature_name", "Signature Name", value="Signature1"),
                           shinyjs::disabled(actionButton("score_signature", "Calculate & Add Score", icon=icon("calculator")))
                  )
      )
    )
  })
  
  # --- SUBSETTING CONTROLS ---
  output$subsetting_controls_ui <- renderUI({
    req(active_object())
    meta_choices <- categorical_metadata()
    all_choices <- c("Numerical Feature (Gene/Score)" = "numerical_feature", meta_choices)
    
    tagList(
      selectInput("subset_by_var", "Subset Using Category:", choices = all_choices),
      conditionalPanel(
        condition = "input.subset_by_var != 'numerical_feature'",
        uiOutput("subset_by_values_ui")
      ),
      conditionalPanel(
        condition = "input.subset_by_var == 'numerical_feature'",
        uiOutput("subset_by_numerical_ui")
      ),
      actionButton("apply_subset", "Apply Subset (Filter)", icon = icon("filter"), class = "btn-info"),
      actionButton("apply_and_reproject", "Apply & Reproject", icon = icon("sync"), class = "btn-warning")
    )
  })
  
  output$subset_by_numerical_ui <- renderUI({
    req(numerical_features_for_subset())
    validate(need(length(numerical_features_for_subset()) > 0,
                  "To filter by numerical feature, first enter a gene in 'Single Marker Plots' or calculate a score."))
    tagList(
      selectInput("subset_numerical_select", "Select Feature to Filter On:", choices = numerical_features_for_subset()),
      uiOutput("subset_numerical_threshold_ui")
    )
  })
  
  # --- SEURAT V5 FIX: NUMERICAL THRESHOLD ---
  output$subset_numerical_threshold_ui <- renderUI({
    req(active_object(), input$subset_numerical_select)
    obj <- active_object()
    feature <- input$subset_numerical_select
    
    if (feature %in% rownames(obj)) {
      obj_for_fetch <- obj
      if ("RNA" %in% Assays(obj_for_fetch)) {
        if (inherits(obj_for_fetch[["RNA"]], "Assay5")) {
          obj_for_fetch[["RNA"]] <- JoinLayers(obj_for_fetch[["RNA"]], assay="RNA")
        }
      }
      vals <- FetchData(obj_for_fetch, vars = feature, layer = "data")[, 1]
    } else if (feature %in% colnames(obj@meta.data)) {
      vals <- obj@meta.data[[feature]]
    } else {
      return(NULL)
    }
    
    min_val <- min(vals, na.rm = TRUE)
    max_val <- max(vals, na.rm = TRUE)
    
    validate(need(is.finite(min_val) && is.finite(max_val), "Cannot determine valid range."),
             need(min_val < max_val, paste0("All cells have same value (", round(min_val, 2), ").")))
    
    sliderInput("subset_numerical_threshold", "Keep cells with value >",
                min = round(min_val, 2), max = round(max_val, 2),
                value = round(min_val + (max_val - min_val) * 0.25, 2), step = 0.05)
  })
  
  output$subset_by_values_ui <- renderUI({
    req(active_object(), input$subset_by_var, input$subset_by_var != 'numerical_feature')
    choices <- sort(unique(as.character(active_object()@meta.data[[input$subset_by_var]])))
    
    tagList(
      div(style="display: flex; justify-content: space-between; align-items: center;",
          tags$b("Select Subgroup(s) to Keep:"),
          div(actionLink("select_all_subsets", "All"), "|", actionLink("select_none_subsets", "None"))
      ),
      checkboxGroupInput("subset_values", label = NULL, choices = choices, selected = choices, inline = TRUE)
    )
  })
  
  output$object_state_controls_ui <- renderUI({
    req(active_object())
    is_initial_state <- length(object_history()) <= 1
    
    tagList(
      hr(style="margin: 10px 0;"),
      div(style="display: flex; justify-content: space-between;",
          actionButton("undo_subset", "Go to Previous", icon = icon("undo"), disabled = is_initial_state),
          actionButton("reset_object", "Go to Initial", icon = icon("fast-backward"))
      ),
      hr(style="margin: 10px 0;"),
      downloadButton("save_subset_rds", "Save Current View (.rds)", icon = icon("download"),
                     style = if(is_initial_state) "color: #a9a9a9; background-color: #f0f0f0; border-color: #dcdcdc; pointer-events: none;" else "")
    )
  })
  
  # --- MISC UI GENERATORS ---
  output$label_overlay_ui_placeholder <- renderUI({
    req(categorical_metadata())
    tagList(
      h5("Labels"),
      checkboxInput("overlay_labels", "Show Labels", value = FALSE),
      conditionalPanel(
        condition = "input.overlay_labels == true",
        selectInput("label_by", "Label Using:", choices = categorical_metadata()),
        sliderInput("overlay_alpha", "Label Alpha:", min = 0, max = 1, value = 0.8, step = 0.1)
      )
    )
  })
  
  output$annotation_ui <- renderUI({
    req(active_object())
    obj <- active_object()
    current_clusters <- sort(unique(obj$seurat_clusters))
    cluster_inputs <- lapply(current_clusters, function(cluster_id) {
      textInput(inputId = paste0("cluster_id_", cluster_id), label = paste("Cluster", cluster_id), value = as.character(cluster_id))
    })
    tagList(
      textInput("annotation_set_name", "New Annotation Set Name:", value = "cell_type"),
      cluster_inputs,
      actionButton("save_annotations", "Save Annotations", icon = icon("save"))
    )
  })
  
  output$isolate_values_ui <- renderUI({
    req(active_object(), input$split_mode == 'isolate', input$split_by_var)
    choices <- sort(unique(as.character(active_object()@meta.data[[input$split_by_var]])))
    checkboxGroupInput("isolate_values", "Select Subgroup(s) to Display:", choices = choices, selected = choices, inline = TRUE)
  })
  
  output$color_category_ui <- renderUI({
    req(categorical_metadata())
    tagList(
      selectInput("color_assign_by", "Coloring Category:", choices = categorical_metadata(), selected = isolate(input$color_by)),
      actionButton("apply_colors", "Apply Colors", icon = icon("paint-brush"))
    )
  })
  
  output$color_picker_ui <- renderUI({
    req(active_object(), input$color_assign_by)
    category <- input$color_assign_by
    subgroups <- sort(unique(as.character(active_object()@meta.data[[category]])))
    palettes <- color_palettes()
    current_palette <- palettes[[category]]
    if (is.null(current_palette)) {
      default_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(subgroups))
      names(default_colors) <- subgroups
      current_palette <- as.list(default_colors)
    }
    lapply(subgroups, function(subgroup) {
      safe_subgroup_name <- make.names(subgroup)
      color_val <- current_palette[[subgroup]] %||% "#808080"
      colourpicker::colourInput(inputId = paste0("col_", category, "_", safe_subgroup_name), label = subgroup, value = color_val)
    })
  })
  
  observeEvent(input$select_all_subsets, {
    req(active_object(), input$subset_by_var)
    choices <- sort(unique(as.character(active_object()@meta.data[[input$subset_by_var]])))
    updateCheckboxGroupInput(session, "subset_values", selected = choices)
  })
  
  observeEvent(input$select_none_subsets, {
    updateCheckboxGroupInput(session, "subset_values", selected = character(0))
  })
  
  output$subgroup_selection_ui <- renderUI({
    req(active_object(), input$color_by)
    meta_col <- active_object()@meta.data[[input$color_by]]
    counts <- table(factor(meta_col))
    subgroups <- names(counts)
    labels_with_counts <- paste0(subgroups, " (", format(counts, big.mark=","), " cells)")
    choices_named <- setNames(subgroups, labels_with_counts)
    checkboxGroupInput("selected_subgroups", "Highlight Subgroups:",
                       choices = choices_named, selected = subgroups, inline = TRUE)
  })
  
  output$feature_plot_gene_inputs_ui <- renderUI({
    req(active_object(), input$n_feature_plots)
    lapply(1:input$n_feature_plots, function(i) {
      selectizeInput(inputId = paste0("feature_gene_", i), label = paste("Feature for Plot", i), choices = NULL,
                     options = list(placeholder = 'Type to search...', maxOptions = 5))
    })
  })
  
  observe({
    req(active_object(), input$n_feature_plots > 0)
    obj <- active_object()
    all_genes <- rownames(obj)
    scores <- signature_scores()
    choices_list <- list()
    if (length(scores) > 0) { choices_list[['Signature Scores']] <- scores }
    choices_list[['Genes']] <- all_genes
    for (i in 1:input$n_feature_plots) {
      updateSelectizeInput(session, paste0("feature_gene_", i), choices = choices_list,
                           selected = isolate(input[[paste0("feature_gene_", i)]]), server = TRUE)
    }
  })
  
  output$reorder_split_ui <- renderUI({
    req(active_object(), input$split_mode == 'split', input$split_by_var)
    items <- sort(unique(as.character(active_object()@meta.data[[input$split_by_var]])))
    orderInput(inputId = "split_order", label = "Drag to reorder facets:", items = items)
  })
  
  observe({ req(active_object()); updateSelectInput(session, "deg_group_by", choices = categorical_metadata()) })
  
  output$deg_group1_ui <- renderUI({
    req(active_object(), input$deg_group_by); choices <- sort(unique(as.character(active_object()@meta.data[[input$deg_group_by]])))
    selectInput("deg_group1", "Group 1 (Test Group):", choices = choices, selected = choices[1])
  })
  
  output$deg_group2_ui <- renderUI({
    req(active_object(), input$deg_group_by, input$deg_group1)
    available_choices <- setdiff(sort(unique(as.character(active_object()@meta.data[[input$deg_group_by]]))), input$deg_group1)
    selectInput("deg_group2", "Group 2 (Control Group):", choices = available_choices, selected = available_choices[1])
  })
  
  output$deg_one_vs_all_groups_ui <- renderUI({
    req(active_object(), input$deg_group_by); choices <- sort(unique(as.character(active_object()@meta.data[[input$deg_group_by]])))
    checkboxGroupInput("deg_one_vs_all_groups", "Select groups to find markers for:", choices = choices, selected = choices, inline = TRUE)
  })
  
  observeEvent(input$run_deg_heatmap, {
    obj <- active_object(); req(obj, input$deg_group_by)
    withProgress(message = "Running DEG analysis...", {
      Idents(obj) <- input$deg_group_by
      if (input$deg_mode == "two_groups") {
        req(input$deg_group1, input$deg_group2)
        markers <- FindMarkers(obj, ident.1 = input$deg_group1, ident.2 = input$deg_group2, verbose = FALSE)
        top_markers <- head(rownames(markers), input$deg_n_genes)
      } else {
        req(input$deg_one_vs_all_groups)
        all_markers <- FindAllMarkers(obj, only.pos = TRUE, verbose = FALSE)
        top_markers <- all_markers %>%
          filter(cluster %in% input$deg_one_vs_all_groups) %>%
          group_by(cluster) %>%
          top_n(n = input$deg_n_genes, wt = avg_log2FC) %>%
          pull(gene) %>%
          unique()
      }
      if (length(top_markers) > 0) {
        deg_heatmap_data(DoHeatmap(obj, features = top_markers, group.by = input$deg_group_by))
      } else { showNotification("No significant markers found.", type="warning"); deg_heatmap_data(NULL) }
    })
  })
  
  # --- TOGGLE OBSERVERS (ANCHORS) ---
  observeEvent(input$toggle_colors_link, { shinyjs::toggle(id = "colors_content", anim = TRUE) })
  observeEvent(input$toggle_annotation_link, { shinyjs::toggle(id = "annotation_content", anim = TRUE) })
  observeEvent(input$toggle_umap_params_link, { shinyjs::toggle(id = "umap_params_content", anim = TRUE) })
  observeEvent(input$toggle_density_link, { shinyjs::toggle(id = "density_content", anim = TRUE) })
  observeEvent(input$toggle_density_color_link, { shinyjs::toggle(id = "density_color_content", anim = TRUE) })
  
  observeEvent(input$`umap-toggle_download_link`, { shinyjs::toggle(id = "umap-download_content", anim = TRUE) })
  observeEvent(input$`feature-toggle_download_link`, { shinyjs::toggle(id = "feature-download_content", anim = TRUE) })
  observeEvent(input$`heatmap-toggle_download_link`, { shinyjs::toggle(id = "heatmap-download_content", anim = TRUE) })
  observeEvent(input$`dotplot-toggle_download_link`, { shinyjs::toggle(id = "dotplot-download_content", anim = TRUE) })
  
  observeEvent(input$greyscale_density_mode, {
    if(isTRUE(input$greyscale_density_mode)) {
      showNotification("Applying greyscale presets...", type="message", duration = 2)
      updateSliderInput(session, "padding_factor", value = 12); updateSliderInput(session, "smoothing", value = 0.5)
      updateSliderInput(session, "glow_opacity", value = 0.23); updateSliderInput(session, "glow_extent", value = 3)
      updateSliderInput(session, "glow_outline_width", value = 0.2); updateSliderInput(session, "core_opacity", value = 0.07)
      updateSliderInput(session, "core_size", value = 18); updateSliderInput(session, "darken_factor", value = 0.67)
      updateSliderInput(session, "core_outline_width", value = 0.1); updateSliderInput(session, "dense_core_opacity", value = 0.3)
      updateSliderInput(session, "dense_core_size", value = 5); updateSliderInput(session, "dense_core_darkening", value = 0.56)
      updateSliderInput(session, "dense_core_outline_width", value = 0.2)
    }
  })
  
  observeEvent(input$rerun_umap, {
    req(active_object())
    withProgress(message = 'Re-running UMAP...', value = 0.3, {
      obj <- active_object()
      reduction_to_use <- active_reduction_dims()
      dims_to_use <- 1:min(20, ncol(obj@reductions[[reduction_to_use]]))
      obj <- RunUMAP(obj,
                     reduction = reduction_to_use, dims = dims_to_use,
                     n.neighbors = as.integer(input$umap_n_neighbors),
                     min.dist = input$umap_min_dist, spread = input$umap_spread,
                     verbose = FALSE)
      active_object(obj)
      showNotification("UMAP recalculated.", type = "message")
    })
  })
  
  parse_genes <- function(gene_string) {
    genes <- unlist(stringr::str_extract_all(gene_string, "[A-Za-z0-9.-]+"))
    return(unique(genes[genes != ""]))
  }
  
  validate_gene_list <- function(gene_vector, all_genes) {
    results <- list(); invalid_genes_found <- FALSE
    for (gene in gene_vector) {
      if (gene %in% all_genes) {
        results[[gene]] <- list(valid = TRUE, suggestion = "")
      } else {
        invalid_genes_found <- TRUE
        distances <- adist(gene, all_genes, ignore.case = TRUE)
        min_dist <- min(distances)
        if (min_dist <= 2) {
          suggestion <- all_genes[which.min(distances)]
          results[[gene]] <- list(valid = FALSE, suggestion = suggestion)
        } else {
          results[[gene]] <- list(valid = FALSE, suggestion = "")
        }
      }
    }
    return(list(results = results, all_valid = !invalid_genes_found))
  }
  
  pos_genes_status <- reactiveVal(list(all_valid = TRUE)); neg_genes_status <- reactiveVal(list(all_valid = TRUE))
  debounced_pos_genes <- reactive(input$positive_markers_csv) %>% debounce(500)
  debounced_neg_genes <- reactive(input$negative_markers_csv) %>% debounce(500)
  
  observe({
    req(active_object(), debounced_pos_genes()); all_genes <- rownames(active_object())
    parsed <- parse_genes(debounced_pos_genes())
    if (length(parsed) == 0) { pos_genes_status(list(all_valid = TRUE, results = list())); return() }
    validation_result <- validate_gene_list(parsed, all_genes); pos_genes_status(validation_result)
  })
  
  observe({
    req(active_object(), debounced_neg_genes()); all_genes <- rownames(active_object())
    parsed <- parse_genes(debounced_neg_genes())
    if (length(parsed) == 0) { neg_genes_status(list(all_valid = TRUE, results = list())); return() }
    validation_result <- validate_gene_list(parsed, all_genes); neg_genes_status(validation_result)
  })
  
  output$positive_marker_validation_ui <- renderUI({
    validation <- pos_genes_status(); if (length(validation$results) == 0) return(NULL)
    tags$div(class = "gene-validation-box", lapply(names(validation$results), function(gene) {
      res <- validation$results[[gene]]
      if (res$valid) { tags$span(class = "text-success", style="margin-right: 5px;", paste0(gene, " (✓)"))
      } else {
        suggestion_text <- if(res$suggestion != "") paste0("did you mean '", res$suggestion, "'?") else ""
        tags$span(class = "text-danger", style="margin-right: 5px;", title=suggestion_text, paste0(gene, " (✖)"))
      }
    }))
  })
  
  output$negative_marker_validation_ui <- renderUI({
    validation <- neg_genes_status(); if (length(validation$results) == 0) return(NULL)
    tags$div(class = "gene-validation-box", lapply(names(validation$results), function(gene) {
      res <- validation$results[[gene]]
      if (res$valid) { tags$span(class = "text-success", style="margin-right: 5px;", paste0(gene, " (✓)"))
      } else {
        suggestion_text <- if(res$suggestion != "") paste0("did you mean '", res$suggestion, "'?") else ""
        tags$span(class = "text-danger", style="margin-right: 5px;", title=suggestion_text, paste0(gene, " (✖)"))
      }
    }))
  })
  
  observe({
    pos_valid <- pos_genes_status()$all_valid; neg_valid <- neg_genes_status()$all_valid
    pos_genes_present <- length(pos_genes_status()$results) > 0
    if (pos_valid && neg_valid && pos_genes_present) { shinyjs::enable("score_signature")
    } else { shinyjs::disable("score_signature") }
  })
  
  observeEvent(input$score_signature, {
    req(active_object(), input$positive_markers_csv, input$signature_name)
    withProgress(message = "Calculating score...", {
      obj <- active_object(); sig_name <- make.names(input$signature_name)
      pos_genes_found <- parse_genes(input$positive_markers_csv); neg_genes_found <- parse_genes(input$negative_markers_csv)
      req(length(pos_genes_found) > 0)
      feature_list <- list(pos_genes_found)
      if (length(neg_genes_found) > 0) { feature_list <- list(Positive = pos_genes_found, Negative = neg_genes_found) }
      obj <- Seurat::AddModuleScore(object = obj, features = feature_list, name = sig_name)
      actual_score_name <- paste0(sig_name, "1")
      signature_scores(unique(c(signature_scores(), actual_score_name)))
      active_object(obj)
      showNotification(paste("Score added:", actual_score_name), type="message", duration=8)
    })
  })
  
  # --- SEURAT V5 FIX: PERFORM SUBSET ---
  perform_subset <- function() {
    obj <- active_object(); req(input$subset_by_var)
    if (input$subset_by_var == "numerical_feature") {
      req(input$subset_numerical_select, input$subset_numerical_threshold)
      feature <- input$subset_numerical_select; threshold <- input$subset_numerical_threshold
      if (feature %in% rownames(obj)) {
        obj_for_fetch <- obj
        if ("RNA" %in% Assays(obj_for_fetch)) {
          if(inherits(obj_for_fetch[["RNA"]], "Assay5")) {
            obj_for_fetch[["RNA"]] <- JoinLayers(obj_for_fetch[["RNA"]], assay="RNA")
          }
        }
        values <- FetchData(obj_for_fetch, vars = feature, layer = "data")[, 1]
      } else if (feature %in% colnames(obj@meta.data)) { values <- obj@meta.data[[feature]]
      } else { showNotification(paste("Error: Selected feature", feature, "not found."), type="error"); return(NULL) }
      cells_to_keep <- Cells(obj)[values > threshold]
    } else {
      req(input$subset_values, message = "Please select at least one subgroup to subset.")
      cells_to_keep <- Cells(obj)[obj@meta.data[[input$subset_by_var]] %in% input$subset_values]
    }
    if (length(cells_to_keep) == 0) { showNotification("Subset resulted in 0 cells.", type="warning"); return(NULL) }
    return(subset(obj, cells = cells_to_keep))
  }
  
  observeEvent(input$apply_subset, {
    new_obj <- perform_subset(); req(new_obj)
    current_history <- object_history(); object_history(c(current_history, list(new_obj)))
    active_object(new_obj)
    showNotification(paste("Subset applied.", ncol(new_obj), "cells remain."), type = "message")
  })
  
  # --- SEURAT V5 FIX: REPROJECT ---
  observeEvent(input$apply_and_reproject, {
    tryCatch({
      withProgress(message = 'Reprojecting Subset...', value = 0.1, {
        subset_obj <- perform_subset(); req(subset_obj)
        
        if ("RNA" %in% names(subset_obj@assays)) {
          if(inherits(subset_obj[["RNA"]], "Assay5")) {
            subset_obj[["RNA"]] <- JoinLayers(subset_obj[["RNA"]])
          }
        }
        
        counts <- LayerData(subset_obj, assay = "RNA", layer = "counts"); metadata <- subset_obj@meta.data
        new_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
        DefaultAssay(new_obj) <- "RNA"
        new_obj <- NormalizeData(new_obj, verbose = FALSE)
        new_obj <- FindVariableFeatures(new_obj, verbose = FALSE); new_obj <- ScaleData(new_obj, verbose = FALSE); new_obj <- RunPCA(new_obj, verbose = FALSE)
        harmony_was_used <- !is.null(object_history()[[1]]@reductions[["harmony"]])
        can_run_harmony <- "orig.ident" %in% colnames(new_obj@meta.data) && length(unique(new_obj$orig.ident)) > 1
        if (harmony_was_used && can_run_harmony) {
          new_obj <- harmony::RunHarmony(new_obj, group.by.vars = "orig.ident", reduction = "pca", assay.use = "RNA", reduction.save = "harmony", verbose = FALSE)
        }
        reduction_for_reproject <- if ("harmony" %in% names(new_obj@reductions)) "harmony" else "pca"
        dims_to_use <- 1:min(20, ncol(new_obj@reductions[[reduction_for_reproject]]))
        new_obj <- FindNeighbors(new_obj, dims = dims_to_use, reduction = reduction_for_reproject, verbose = FALSE)
        new_obj <- FindClusters(new_obj, resolution = input$resolution, verbose = FALSE)
        new_obj <- RunUMAP(new_obj, dims = dims_to_use, reduction = reduction_for_reproject, verbose = FALSE)
        current_history <- object_history(); object_history(c(current_history, list(new_obj)))
        active_object(new_obj)
        showNotification(paste("Subset reprojected.", ncol(new_obj), "cells."), type = "message")
      })
    }, error = function(e) { showNotification(paste("Reprojection Error:", e$message), type = "error", duration = 15) })
  })
  
  observeEvent(input$undo_subset, {
    hist <- object_history(); if (length(hist) > 1) {
      new_hist <- head(hist, -1); prev_obj <- new_hist[[length(new_hist)]]
      object_history(new_hist); active_object(prev_obj)
      showNotification("Undo successful.", type = "message")
    }
  })
  
  observeEvent(input$reset_object, {
    hist <- object_history(); if(length(hist) > 0) {
      initial_obj <- hist[[1]]; object_history(list(initial_obj)); active_object(initial_obj)
      showNotification("Reset to initial state.", type = "message")
    }
  })
  
  output$save_subset_rds <- downloadHandler(
    filename = function() { paste0("seurat_subset_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds") },
    content = function(file) { obj_to_save <- active_object(); req(obj_to_save); saveRDS(obj_to_save, file = file) }
  )
  
  observeEvent(input$apply_colors, {
    req(input$color_assign_by); category <- input$color_assign_by
    subgroups <- sort(unique(as.character(active_object()@meta.data[[category]])))
    current_palettes <- color_palettes()
    new_palette <- lapply(subgroups, function(subgroup) {
      safe_subgroup_name <- make.names(subgroup)
      input[[paste0("col_", category, "_", safe_subgroup_name)]]
    })
    names(new_palette) <- subgroups
    current_palettes[[category]] <- new_palette; color_palettes(current_palettes)
    showNotification("Custom colors applied.", type = "message")
  })
  
  observeEvent(input$run_clustering, {
    obj <- active_object(); req(obj, input$resolution)
    withProgress(message = 'Re-clustering...', {
      dims_to_use <- 1:min(20, ncol(obj@reductions[[active_reduction_dims()]]))
      obj <- FindNeighbors(obj, dims = dims_to_use, reduction = active_reduction_dims(), verbose = FALSE)
      obj <- FindClusters(obj, resolution = input$resolution, verbose = FALSE)
      active_object(obj)
      showNotification(paste("Clustering updated (Res:", input$resolution, ")"), type = "message")
    })
  })
  
  observeEvent(input$save_annotations, {
    req(active_object(), input$annotation_set_name)
    obj <- active_object(); set_name <- make.names(input$annotation_set_name)
    current_clusters <- sort(unique(obj$seurat_clusters))
    id_map <- sapply(current_clusters, function(cluster_id) input[[paste0("cluster_id_", cluster_id)]])
    names(id_map) <- current_clusters; obj@meta.data[[set_name]] <- id_map[as.character(obj$seurat_clusters)]
    current_history <- object_history(); current_history[[length(current_history)]] <- obj
    object_history(current_history); active_object(obj)
    showNotification(paste("Annotations saved to:", set_name), type = "message")
  })
  
  # --- PLOTTING DATA GENERATORS ---
  plot_data_reactive <- reactive({
    obj <- active_object(); req(obj, input$color_by, !is.null(input$greyscale_density_mode))
    if (isTRUE(input$greyscale_density_mode)) {
      obj$plot_group <- "background_points"; obj$plot_group <- factor(obj$plot_group)
    } else {
      req(input$selected_subgroups)
      obj$plot_group <- as.character(obj@meta.data[[input$color_by]])
      obj$plot_group[!obj$plot_group %in% input$selected_subgroups] <- "Unselected"
      obj$plot_group <- factor(obj$plot_group, levels = c(sort(input$selected_subgroups), "Unselected"))
    }
    return(obj)
  })
  
  density_data_reactive <- reactive({
    obj <- active_object(); req(obj, input$color_by, input$selected_subgroups)
    umap_coords <- Embeddings(obj, reduction = active_reduction_umap()); colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    density_df <- cbind(obj@meta.data, umap_coords)
    density_df$plot_group <- as.character(density_df[[input$color_by]])
    density_df$plot_group[!density_df$plot_group %in% input$selected_subgroups] <- "Unselected"
    density_df$plot_group <- factor(density_df$plot_group, levels = c(sort(input$selected_subgroups), "Unselected"))
    return(density_df)
  })
  
  color_palette_reactive <- reactive({
    obj <- plot_data_reactive(); req(obj); category_to_color <- input$color_by; all_palettes <- color_palettes()
    custom_palette <- all_palettes[[category_to_color]]; density_obj_for_levels <- active_object(); req(density_obj_for_levels)
    color_levels <- unique(as.character(density_obj_for_levels@meta.data[[category_to_color]]))
    final_colors <- if (!is.null(custom_palette)) { unlist(custom_palette)
    } else {
      n_colors <- length(color_levels)
      palette <- if (n_colors > 0) colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_colors) else character(0)
      names(palette) <- color_levels; palette
    }
    final_colors <- c(final_colors, `Unselected` = "#d3d3d3", `background_points` = "#d3d3d3")
    return(final_colors)
  })
  
  # --- UMAP PLOT OBJECT ---
  umap_plot_object <- reactive({
    obj_with_colors <- plot_data_reactive(); final_colors <- color_palette_reactive()
    req(obj_with_colors, input$pt_size, input$alpha, input$split_mode)
    plot_obj <- obj_with_colors; split_variable <- NULL
    if (input$split_mode == "isolate") {
      req(input$split_by_var, input$isolate_values)
      cells_to_keep <- Cells(plot_obj)[plot_obj@meta.data[[input$split_by_var]] %in% input$isolate_values]
      validate(need(length(cells_to_keep) > 0, "Isolation resulted in 0 cells."))
      plot_obj <- subset(plot_obj, cells = cells_to_keep)
    }
    if (input$split_mode == "split") {
      req(input$split_by_var, input$split_order); split_variable <- input$split_by_var
      plot_obj@meta.data[[split_variable]] <- factor(plot_obj@meta.data[[split_variable]], levels = input$split_order)
    }
    p <- DimPlot(plot_obj, group.by = "plot_group", reduction = active_reduction_umap(), label = TRUE, repel = TRUE, pt.size = input$pt_size, alpha = input$alpha, cols = final_colors, split.by = split_variable) + NoLegend()
    if (isTRUE(input$show_density)) {
      req(!is.null(input$glow_outline_color)); plot_df <- density_data_reactive()
      if (input$split_mode == "split" && !is.null(split_variable)) {
        plot_df[[split_variable]] <- factor(plot_df[[split_variable]], levels = input$split_order)
        split_levels <- input$split_order
        glow_polygons_list <- list(); core_polygons_list <- list(); dense_core_polygons_list <- list()
        for (current_split_level in split_levels) {
          df_split_subset <- plot_df %>% filter(.data[[split_variable]] == current_split_level)
          if (nrow(df_split_subset) < 10) next
          padding_val_x <- (max(df_split_subset$UMAP_1) - min(df_split_subset$UMAP_1)) * (input$padding_factor/100)
          padding_val_y <- (max(df_split_subset$UMAP_2) - min(df_split_subset$UMAP_2)) * (input$padding_factor/100)
          lims <- c(min(df_split_subset$UMAP_1)-padding_val_x, max(df_split_subset$UMAP_1)+padding_val_x, min(df_split_subset$UMAP_2)-padding_val_y, max(df_split_subset$UMAP_2)+padding_val_y)
          density_polygons_facet <- calculate_density_polygons(data=df_split_subset, x_col="UMAP_1", y_col="UMAP_2", group_col="plot_group", smoothing=input$smoothing, limits=lims, normalize=input$normalize_density)
          if (is.null(density_polygons_facet)) next
          for (group in names(density_polygons_facet)) {
            if (group == "Unselected" || is.null(density_polygons_facet[[group]])) next
            group_color <- final_colors[group]; n_levels <- length(density_polygons_facet[[group]])
            if (isTRUE(input$show_glow) && n_levels > 0) {
              fill <- if(isTRUE(input$glow_fill_use_group_color)) group_color else "#CCCCCC"
              for(l in 1:min(n_levels, input$glow_extent)) { poly_df <- as_tibble(density_polygons_facet[[group]][[l]]); poly_df$fill_color <- fill; poly_df$outline_color <- input$glow_outline_color; poly_df[[split_variable]] <- current_split_level; glow_polygons_list[[length(glow_polygons_list) + 1]] <- poly_df }
            }
            if (isTRUE(input$show_core) && n_levels > 0) {
              fill_base <- if(isTRUE(input$core_fill_use_group_color)) group_color else "#CCCCCC"
              for(l in 1:min(n_levels, input$core_size)) { poly_df <- as_tibble(density_polygons_facet[[group]][[l]]); poly_df$fill_color <- darken_color(fill_base, input$darken_factor); poly_df$outline_color <- input$core_outline_color; poly_df[[split_variable]] <- current_split_level; core_polygons_list[[length(core_polygons_list) + 1]] <- poly_df }
            }
            if (isTRUE(input$show_dense_core) && n_levels > 0) {
              fill_base <- if(isTRUE(input$dense_core_fill_use_group_color)) group_color else "#CCCCCC"
              for(l in max(1, n_levels - input$dense_core_size + 1):n_levels) { poly_df <- as_tibble(density_polygons_facet[[group]][[l]]); poly_df$fill_color <- darken_color(fill_base, input$dense_core_darkening); poly_df$outline_color <- input$dense_core_outline_color; poly_df[[split_variable]] <- current_split_level; dense_core_polygons_list[[length(dense_core_polygons_list) + 1]] <- poly_df }
            }
          }
        }
        if (length(glow_polygons_list) > 0) { glow_df <- bind_rows(glow_polygons_list, .id = "polygon_id"); glow_df[[split_variable]] <- factor(glow_df[[split_variable]], levels = input$split_order); p <- p + geom_polygon(data=glow_df, aes(x=x, y=y, group=polygon_id, fill=I(fill_color), color=I(outline_color)), alpha=input$glow_opacity, linewidth=input$glow_outline_width) }
        if (length(core_polygons_list) > 0) { core_df <- bind_rows(core_polygons_list, .id = "polygon_id"); core_df[[split_variable]] <- factor(core_df[[split_variable]], levels = input$split_order); p <- p + geom_polygon(data=core_df, aes(x=x, y=y, group=polygon_id, fill=I(fill_color), color=I(outline_color)), alpha=input$core_opacity, linewidth=input$core_outline_width) }
        if (length(dense_core_polygons_list) > 0) { dense_core_df <- bind_rows(dense_core_polygons_list, .id = "polygon_id"); dense_core_df[[split_variable]] <- factor(dense_core_df[[split_variable]], levels = input$split_order); p <- p + geom_polygon(data=dense_core_df, aes(x=x, y=y, group=polygon_id, fill=I(fill_color), color=I(outline_color)), alpha=input$dense_core_opacity, linewidth=input$dense_core_outline_width) }
      } else {
        padding_val_x <- (max(plot_df$UMAP_1) - min(plot_df$UMAP_1)) * (input$padding_factor/100)
        padding_val_y <- (max(plot_df$UMAP_2) - min(plot_df$UMAP_2)) * (input$padding_factor/100)
        lims <- c(min(plot_df$UMAP_1)-padding_val_x, max(plot_df$UMAP_1)+padding_val_x, min(plot_df$UMAP_2)-padding_val_y, max(plot_df$UMAP_2)+padding_val_y)
        density_polygons <- calculate_density_polygons(data=plot_df, x_col="UMAP_1", y_col="UMAP_2", group_col="plot_group", smoothing=input$smoothing, limits=lims, normalize=input$normalize_density)
        if (!is.null(density_polygons)) {
          for (group in names(density_polygons)) {
            if (group == "Unselected" || is.null(density_polygons[[group]])) next
            group_color <- final_colors[group]; n_levels <- length(density_polygons[[group]])
            if (isTRUE(input$show_glow) && n_levels > 0) {
              fill <- if(isTRUE(input$glow_fill_use_group_color)) group_color else "#CCCCCC"
              for(l in 1:min(n_levels, input$glow_extent)) { p <- p + geom_polygon(data=data.frame(x=density_polygons[[group]][[l]]$x, y=density_polygons[[group]][[l]]$y), aes(x=x, y=y), fill=fill, color=input$glow_outline_color, alpha=input$glow_opacity, linewidth=input$glow_outline_width) }
            }
            if (isTRUE(input$show_core) && n_levels > 0) {
              fill_base <- if(isTRUE(input$core_fill_use_group_color)) group_color else "#CCCCCC"
              for(l in 1:min(n_levels, input$core_size)) { p <- p + geom_polygon(data=data.frame(x=density_polygons[[group]][[l]]$x, y=density_polygons[[group]][[l]]$y), aes(x=x, y=y), fill=darken_color(fill_base, input$darken_factor), color=input$core_outline_color, alpha=input$core_opacity, linewidth=input$core_outline_width) }
            }
            if (isTRUE(input$show_dense_core) && n_levels > 0) {
              fill_base <- if(isTRUE(input$dense_core_fill_use_group_color)) group_color else "#CCCCCC"
              for(l in max(1, n_levels - input$dense_core_size + 1):n_levels) { p <- p + geom_polygon(data=data.frame(x=density_polygons[[group]][[l]]$x, y=density_polygons[[group]][[l]]$y), aes(x=x, y=y), fill=darken_color(fill_base, input$dense_core_darkening), color=input$dense_core_outline_color, alpha=input$dense_core_opacity, linewidth=input$dense_core_outline_width) }
            }
          }
        }
      }
    }
    if (isTRUE(input$overlay_labels)) {
      req(input$label_by); label_data <- plot_obj@meta.data
      label_data$UMAP_1 <- Embeddings(plot_obj, reduction = active_reduction_umap())[,1]
      label_data$UMAP_2 <- Embeddings(plot_obj, reduction = active_reduction_umap())[,2]
      label_centers <- label_data %>% group_by(.data[[input$label_by]]) %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = 'drop')
      p <- p + ggrepel::geom_label_repel(data = label_centers, aes(x=UMAP_1, y=UMAP_2, label=.data[[input$label_by]]), alpha=input$overlay_alpha, color="white", fill="black", inherit.aes=FALSE)
    }
    return(p)
  })
  
  feature_plot_object <- reactive({
    obj_to_use <- active_object(); req(obj_to_use, input$n_feature_plots)
    features_to_plot <- sapply(1:input$n_feature_plots, function(i) input[[paste0("feature_gene_", i)]])
    valid_features_in_object <- c(rownames(obj_to_use), colnames(obj_to_use@meta.data))
    features_to_plot <- features_to_plot[ !is.null(features_to_plot) & features_to_plot != "" & features_to_plot %in% valid_features_in_object ]
    if (length(features_to_plot) == 0) return(ggplot() + theme_void() + ggtitle("Enter a valid gene or select a score..."))
    
    plot_obj <- obj_to_use; split_variable <- NULL
    if (isTRUE(input$apply_mode_to_features)) {
      if (input$split_mode == "isolate") {
        req(input$split_by_var, input$isolate_values)
        cells_to_keep <- Cells(plot_obj)[plot_obj@meta.data[[input$split_by_var]] %in% input$isolate_values]
        plot_obj <- subset(plot_obj, cells = cells_to_keep)
      }
      if (input$split_mode == "split") {
        req(input$split_by_var, input$split_order); split_variable <- input$split_by_var
        plot_obj@meta.data[[split_variable]] <- factor(plot_obj@meta.data[[split_variable]], levels = input$split_order)
      }
    }
    
    fp_args <- list(object = plot_obj, features = features_to_plot, reduction = active_reduction_umap(),
                    pt.size = input$pt_size, order = TRUE, alpha = input$alpha, split.by = split_variable)
    
    if (isTRUE(input$enable_custom_scale)) {
      req(input$feature_min_cutoff, input$feature_max_cutoff)
      fp_args$min.cutoff <- paste0('q', input$feature_min_cutoff)
      fp_args$max.cutoff <- paste0('q', input$feature_max_cutoff)
    }
    
    p <- do.call(FeaturePlot, fp_args)
    
    if (req(input$feature_color_palette) == "Viridis") {
      p <- p & scale_color_viridis()
    } else if (input$feature_color_palette != "Default (Grey->Blue)") {
      p <- p & scale_color_distiller(palette = input$feature_color_palette)
    }
    
    return(p)
  })
  
  # --- PLOT RENDERING ---
  output$umap_plot <- renderPlot({ umap_plot_object() })
  output$feature_plot_grid <- renderPlot({ feature_plot_object() })
  output$dot_plot <- renderPlot({ dot_plot_object() })
  output$dist_analysis_bar_plot <- renderPlot({ dist_bar_plot_object() })
  output$dist_analysis_donut_plot <- renderPlot({ dist_donut_plot_object() })
  output$deg_heatmap_plot <- renderPlot({ deg_heatmap_data() })
  
  # --- DOT PLOT SECTION ---
  output$dot_plot_controls_ui <- renderUI({
    req(active_object())
    all_palettes <- rownames(RColorBrewer::brewer.pal.info)
    color_choices <- c("Viridis", "Default (RdBu)", all_palettes)
    tagList(
      selectInput("dotplot_group_by", "Group Cells Using:", choices = categorical_metadata()),
      selectInput("dotplot_color_palette", "Color Palette:", choices = color_choices, selected = "Viridis"),
      hr(),
      helpText("Features for the Dot Plot are selected from the 'Single Marker Plots' section.")
    )
  })
  
  dot_plot_object <- reactive({
    obj_to_use <- active_object(); req(obj_to_use, input$dotplot_group_by)
    features_to_plot <- sapply(1:input$n_feature_plots, function(i) input[[paste0("feature_gene_", i)]])
    valid_features_in_object <- c(rownames(obj_to_use), colnames(obj_to_use@meta.data))
    features_to_plot <- features_to_plot[ !is.null(features_to_plot) & features_to_plot != "" & features_to_plot %in% valid_features_in_object ]
    validate(need(length(features_to_plot) > 0, "Please select at least one valid feature."))
    
    plot_obj <- obj_to_use; split_variable <- NULL
    if (isTRUE(input$apply_mode_to_dot_plots)) {
      if (input$split_mode == "isolate") {
        req(input$split_by_var, input$isolate_values)
        cells_to_keep <- Cells(plot_obj)[plot_obj@meta.data[[input$split_by_var]] %in% input$isolate_values]
        plot_obj <- subset(plot_obj, cells = cells_to_keep)
      }
      if (input$split_mode == "split") {
        req(input$split_by_var, input$split_order)
        split_variable <- input$split_by_var
        plot_obj@meta.data[[split_variable]] <- factor(plot_obj@meta.data[[split_variable]], levels = input$split_order)
      }
    }
    
    p <- DotPlot(plot_obj, features = features_to_plot, group.by = input$dotplot_group_by, split.by = split_variable) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (req(input$dotplot_color_palette) == "Viridis") {
      p <- p + scale_color_viridis()
    } else if (input$dotplot_color_palette != "Default (RdBu)") {
      p <- p + scale_color_distiller(palette = input$dotplot_color_palette)
    }
    
    return(p)
  })
  
  # --- DISTRIBUTION ANALYSIS ---
  output$dist_analysis_controls_ui <- renderUI({
    req(active_object())
    choices <- categorical_metadata()[categorical_metadata() != "seurat_clusters"]
    tagList(
      radioButtons("dist_analysis_mode", "Analysis Mode:",
                   choices = c("Cluster Makeup per Group" = "clusters_in_group", "Group Makeup per Cluster" = "groups_in_cluster")),
      selectInput("dist_analysis_var", "Group data by:", choices = choices),
      uiOutput("dist_analysis_reorder_ui")
    )
  })
  
  output$dist_analysis_reorder_ui <- renderUI({
    req(active_object(), input$dist_analysis_var)
    items <- sort(unique(as.character(active_object()@meta.data[[input$dist_analysis_var]])))
    orderInput(inputId = "dist_analysis_order", label = "Drag to reorder groups:", items = items)
  })
  
  dist_analysis_data_reactive <- reactive({
    obj <- active_object(); req(obj, input$dist_analysis_var, input$dist_analysis_mode)
    comp_var <- input$dist_analysis_var; custom_order <- input$dist_analysis_order
    base_data <- obj@meta.data %>% as_tibble()
    is_order_valid <- !is.null(custom_order) && all(sort(unique(as.character(base_data[[comp_var]]))) == sort(custom_order))
    if (input$dist_analysis_mode == "clusters_in_group") {
      df <- base_data %>% count(.data[[comp_var]], seurat_clusters, name = "n") %>% group_by(.data[[comp_var]]) %>% mutate(proportion = n / sum(n)) %>% ungroup()
    } else {
      df <- base_data %>% count(seurat_clusters, .data[[comp_var]], name = "n") %>% group_by(seurat_clusters) %>% mutate(proportion = n / sum(n)) %>% ungroup()
    }
    df <- df %>% mutate(seurat_clusters = factor(seurat_clusters), !!comp_var := factor(.data[[comp_var]], levels = if(is_order_valid) custom_order else sort(unique(as.character(.data[[comp_var]])))))
    return(df)
  })
  
  dist_analysis_colors_reactive <- reactive({
    data <- dist_analysis_data_reactive(); req(data)
    custom_cluster_colors <- color_palettes()[['seurat_clusters']]; cluster_levels <- levels(data$seurat_clusters)
    if (is.null(custom_cluster_colors)) { setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(cluster_levels)), cluster_levels)
    } else { unlist(custom_cluster_colors) }
  })
  
  dist_bar_plot_object <- reactive({
    composition_data <- dist_analysis_data_reactive(); req(composition_data)
    comp_var <- isolate(input$dist_analysis_var)
    if (isolate(input$dist_analysis_mode) == "clusters_in_group") {
      ggplot(composition_data, aes(x = .data[[comp_var]], y = proportion, fill = seurat_clusters)) +
        geom_col(position = "stack", color = "white") + scale_fill_manual(values = dist_analysis_colors_reactive(), name = "Cluster") +
        scale_y_continuous(labels = scales::percent_format()) + labs(title = "Cluster Composition per Group", x = tools::toTitleCase(gsub("_", " ", comp_var)), y = "Proportion") +
        theme_minimal(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))
    } else {
      ggplot(composition_data, aes(x = seurat_clusters, y = proportion, fill = .data[[comp_var]])) +
        geom_col(position = "stack", color = "white") + scale_fill_brewer(palette = "Set3", name = tools::toTitleCase(gsub("_", " ", comp_var))) +
        scale_y_continuous(labels = scales::percent_format()) + labs(title = "Group Composition per Cluster", x = "Cluster", y = "Proportion") + theme_minimal(base_size = 14)
    }
  })
  
  dist_donut_plot_object <- reactive({
    composition_data <- dist_analysis_data_reactive(); req(composition_data)
    comp_var <- isolate(input$dist_analysis_var)
    if (isolate(input$dist_analysis_mode) == "clusters_in_group") {
      df <- composition_data %>% group_by(.data[[comp_var]]) %>% arrange(desc(seurat_clusters)) %>% mutate(lab.ypos = cumsum(proportion) - 0.5 * proportion)
      ggplot(df, aes(x=2, y=proportion, fill=seurat_clusters)) + geom_col(color="white") + coord_polar(theta="y") + facet_wrap(vars(.data[[comp_var]])) +
        geom_text(aes(y=lab.ypos, label=seurat_clusters), color="white") + scale_fill_manual(values=dist_analysis_colors_reactive()) + theme_void() + xlim(0.5, 2.5)
    } else {
      df <- composition_data %>% group_by(seurat_clusters) %>% arrange(desc(.data[[comp_var]])) %>% mutate(lab.ypos = cumsum(proportion) - 0.5 * proportion)
      ggplot(df, aes(x=2, y=proportion, fill=.data[[comp_var]])) + geom_col(color="white") + coord_polar(theta="y") + facet_wrap(vars(seurat_clusters)) +
        geom_text(aes(y=lab.ypos, label=.data[[comp_var]]), color="black", size=3) + scale_fill_brewer(palette="Set3") + theme_void() + xlim(0.5, 2.5)
    }
  })
  
  # --- PLOT RESIZING OBSERVERS ---
  observeEvent(input$umap_plot_size, {
    size <- input$umap_plot_size; req(size)
    updateNumericInput(session, "umap-pdf_width", value = round(size$width / 96, 1))
    updateNumericInput(session, "umap-pdf_height", value = round(size$height / 96, 1))
  })
  
  observeEvent(input$feature_plot_grid_size, {
    size <- input$feature_plot_grid_size; req(size)
    updateNumericInput(session, "feature-pdf_width", value = round(size$width / 96, 1))
    updateNumericInput(session, "feature-pdf_height", value = round(size$height / 96, 1))
  })
  
  observeEvent(input$deg_heatmap_plot_size, {
    size <- input$deg_heatmap_plot_size; req(size)
    updateNumericInput(session, "heatmap-pdf_width", value = round(size$width / 96, 1))
    updateNumericInput(session, "heatmap-pdf_height", value = round(size$height / 96, 1))
  })
  
  observeEvent(input$dot_plot_size, {
    size <- input$dot_plot_size; req(size)
    updateNumericInput(session, "dotplot-pdf_width", value = round(size$width / 96, 1))
    updateNumericInput(session, "dotplot-pdf_height", value = round(size$height / 96, 1))
  })
  
  observe({
    size <- if (req(input$distribution_plot_tabs) == "Bar Plot") { input$dist_analysis_bar_plot_size } else { input$dist_analysis_donut_plot_size }
    req(size)
    updateNumericInput(session, "dist_analysis_pdf_width", value = round(size$width / 96, 1))
    updateNumericInput(session, "dist_analysis_pdf_height", value = round(size$height / 96, 1))
  })
  
  # --- DOWNLOAD HANDLERS ---
  setup_download_handler <- function(input, output, session, id_prefix, plot_object_reactive) {
    output$download_button <- downloadHandler(
      filename = function() {
        req(input$filename)
        if (!endsWith(input$filename, ".pdf")) paste0(input$filename, ".pdf") else input$filename
      },
      content = function(file) {
        withProgress(message = 'Generating PDF...', {
          plot_to_save <- plot_object_reactive(); req(plot_to_save)
          if (isTRUE(input$rasterize)) { try({ plot_to_save <- ggrastr::rasterise(plot_to_save, dpi = 300) }, silent = TRUE) }
          ggsave(file, plot_to_save, width = input$pdf_width, height = input$pdf_height, device = "pdf", useDingbats = FALSE)
        })
      }
    )
  }
  
  callModule(setup_download_handler, "umap", plot_object_reactive = umap_plot_object)
  callModule(setup_download_handler, "feature", plot_object_reactive = feature_plot_object)
  callModule(setup_download_handler, "heatmap", plot_object_reactive = deg_heatmap_data)
  callModule(setup_download_handler, "dotplot", plot_object_reactive = dot_plot_object)
  
  output$download_dist_analysis_plot <- downloadHandler(
    filename = function() {
      req(input$dist_analysis_filename)
      if (!endsWith(input$dist_analysis_filename, ".pdf")) paste0(input$dist_analysis_filename, ".pdf") else input$dist_analysis_filename
    },
    content = function(file) {
      req(input$distribution_plot_tabs)
      plot_to_save <- if (input$distribution_plot_tabs == "Bar Plot") { dist_bar_plot_object() } else { dist_donut_plot_object() }
      ggsave(file, plot_to_save, width=input$dist_analysis_pdf_width, height=input$dist_analysis_pdf_height, device="pdf", useDingbats = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)