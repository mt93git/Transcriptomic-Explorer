# =============================================================================
# NeuMap Explorer :: ui.R
#
# DESCRIPTION:
# This script defines the user interface (UI) for the application. It lays
# out all the visible elements, tabs, sidebars, panels, inputs, and outputs.
# =============================================================================

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "yeti"),
  useShinyjs(),
  shiny::tags$head(
    shiny::tags$style(HTML(".rank-list-title { font-size: 0.85em; }")),
    shiny::tags$style(HTML(".btn-sm { padding: 0.2rem 0.4rem; font-size: 0.8rem; }"))
  ),
  titlePanel("NeuMap Explorer (v22.9.0)"),
  tabsetPanel(
    id = "main_tabs",
    tabPanel(
      "Interactive Analysis Dashboard",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          style = "max-height: 85vh; overflow-y: auto;",
          h3("Analysis Controls"),
          h4("1. Select Cell Population"),
          p("Use the checkboxes to filter cells.", style = "font-size: 0.9em; color: #666;"),
          create_filter_panel("filter", DASHBOARD_FILTERS, neumap_data_base@meta.data),
          checkboxInput("cumulative_filter", "Enable Combinatorial (AND) Filtering", value = TRUE),
          uiOutput("filter_mode_explanation_ui"),
          actionButton("apply_filters", "Apply Selections", icon = icon("search"), class = "btn-primary btn-block"),
          actionButton("reset_filters", "Reset Selections", icon = icon("refresh"), class = "btn-default btn-block"),
          uiOutput("cell_count_summary_ui"),
          hr(),
          h4("2. Select Feature to Visualize"),
          radioButtons("analysis_type", "Feature Type:", choices = c("Gene Expression", "Signature Score"), selected = "Gene Expression"),
          conditionalPanel(
            condition = "input.analysis_type == 'Gene Expression'",
            selectizeInput("explorer_gene_select", "Select Gene:", choices = NULL)
          ),
          conditionalPanel(
            condition = "input.analysis_type == 'Signature Score'",
            uiOutput("signature_score_ui_main")
          ),
          hr(),
          h5("Threshold for 'Positive' Cells"),
          sliderInput("explorer_threshold", NULL, min = 0, max = 5, value = 1.0, step = 0.1),
          hr(),
          h4("Global Data Capping Controls"),
          p("These controls affect all plots in all tabs.", style = "font-size: 0.9em; color: #666;"),
          checkboxInput("enable_capping", "Enable Global Data Capping", value = TRUE),
          conditionalPanel(
            condition = "input.enable_capping == true",
            radioButtons("capping_mode", "Mode:", choices = c("Cap color range (Winsorize)", "Hide outlier cells (Trim)"), selected = "Cap color range (Winsorize)"),
            fluidRow(
              column(6, numericInput("lower_percentile", "Lower %", value = 5, min = 0, max = 49, step = 1)),
              column(6, numericInput("upper_percentile", "Upper %", value = 95, min = 51, max = 100, step = 1))
            )
          ),
          create_color_menu_ui("tab1")
        ),
        mainPanel(
          width = 9,
          conditionalPanel(
            condition = "input.cumulative_filter == true",
            wellPanel(
              h4("Refine Analysis Population"),
              p("The checkbox selections on the left generate a list of all unique metadata combinations. You can deselect specific combinations here before running the analysis.", style = "font-size: 0.9em; color: #666;"),
              checkboxInput("explorer_include_hub_status", "Include 'Hub Status' in Combinations", value = FALSE),
              hr(),
              uiOutput("explorer_combination_selector_ui")
            )
          ),
          wellPanel(
            h4("UMAP Visualization Controls"),
            fluidRow(
              column(
                6,
                h5("Featured Expression UMAP"),
                radioButtons("featured_umap_mode", "Display Mode:",
                  choices = c("Show all cells" = "all", "Highlight selected population" = "highlight"),
                  selected = "highlight", inline = TRUE
                ),
                # --- Updated Labels and Added Slider ---
                sliderInput("featured_alpha", "Dot Transparency:", min = 0.1, max = 1.0, value = 0.8, step = 0.1),
                sliderInput("featured_pt_size", "Dot Size:", min = 0.1, max = 3.0, value = 1.5, step = 0.1) # Default size adjusted slightly for consistency
              ),
              column(
                6,
                h5("Reference UMAP"),
                selectInput("reference_umap_variable", "Reference Variable:", choices = ANNOTATION_CATEGORIES_RNA, selected = if ("hub_status" %in% ANNOTATION_CATEGORIES_RNA) "hub_status" else ANNOTATION_CATEGORIES_RNA[1]),
                radioButtons("reference_umap_mode", "Display Mode:",
                  choices = c("Show all cells" = "all", "Highlight selected population" = "highlight"),
                  selected = "highlight", inline = TRUE
                ),
                # --- Updated Labels and Added Slider ---
                sliderInput("reference_alpha", "Dot Transparency:", min = 0.1, max = 1.0, value = 0.7, step = 0.1),
                sliderInput("reference_pt_size", "Dot Size:", min = 0.1, max = 3.0, value = 1.0, step = 0.1)
              )
            )
          ),
          fluidRow(
            column(
              6,
              h4("Featured Expression UMAP"),
              p(style = "font-size: 0.9em; color: #6c757d;", "Expression level of the selected gene or signature score across all cells."),
              withSpinner(jqui_resizable(plotOutput("featured_expression_umap", height = "500px"))),
              div(style = "text-align: right; margin-top: 5px;", downloadButton("download_featured_umap", "Download PDF", class = "btn-outline-primary btn-sm"))
            ),
            column(
              6,
              h4("Reference UMAP"),
              p(style = "font-size: 0.9em; color: #6c757d;", "Categorical metadata variable for anatomical or functional context."),
              withSpinner(jqui_resizable(plotOutput("reference_umap", height = "500px"))),
              div(style = "text-align: right; margin-top: 5px;", downloadButton("download_reference_umap", "Download PDF", class = "btn-outline-primary btn-sm"))
            )
          ),
          hr(),
          fluidRow(
            column(
              12,
              h4("Subgroup Location UMAP"),
              p(style = "font-size: 0.9em; color: #6c757d;", "Spatial location of user-defined cell populations on the reference UMAP."),
              withSpinner(jqui_resizable(plotOutput("subgroup_location_umap", height = "500px"))),
              div(style = "text-align: right; margin-top: 5px;", downloadButton("download_subgroup_umap", "Download PDF", class = "btn-outline-primary btn-sm"))
            )
          ),
          hr(),
          fluidRow(
            column(
              12,
              fluidRow(
                column(12, checkboxInput("show_barchart_labels_tab1", "Show Value Labels on Bar Charts", value = TRUE))
              ),
              fluidRow(
                column(
                  6,
                  h4("% Positive Cells by Subgroup"),
                  p(style = "font-size: 0.9em; color: #6c757d;", "Percentage of cells in each subgroup exceeding the defined expression threshold."),
                  withSpinner(jqui_resizable(plotOutput("percent_positive_barchart", height = "350px"))),
                  div(style = "text-align: right; margin-top: 5px;", downloadButton("download_percent_barchart", "Download PDF", class = "btn-outline-primary btn-sm"))
                ),
                column(
                  6,
                  h4("Mean Expression in Positive Cells"),
                  p(style = "font-size: 0.9em; color: #6c757d;", "Average expression level for 'positive' cells within each subgroup."),
                  withSpinner(jqui_resizable(plotOutput("mean_expression_barchart", height = "350px"))),
                  div(style = "text-align: right; margin-top: 5px;", downloadButton("download_mean_barchart", "Download PDF", class = "btn-outline-primary btn-sm"))
                )
              )
            )
          ),
          hr(),
          fluidRow(
            column(
              12,
              h4("Expression Distribution by Subgroup"),
              p(style = "font-size: 0.9em; color: #6c757d;", "Full distribution of expression values for each subgroup."),
              withSpinner(jqui_resizable(plotOutput("expression_distribution_violin", height = "400px"))),
              div(style = "text-align: right; margin-top: 5px;", downloadButton("download_violin_plot_tab1", "Download PDF", class = "btn-outline-primary btn-sm"))
            )
          )
        )
      )
    ),
    tabPanel(
      "Differential Gene Discovery",
      sidebarLayout(
        sidebarPanel(
          width = 3,
          style = "max-height: 85vh; overflow-y: auto;",
          h3("Differential Expression Setup"), hr(),
          h4("1. Define Cell Population"),
          p("Use checkboxes to select the cells for analysis. Selections are combined with AND logic.", style = "font-size: 0.9em; color: #666;"),
          create_filter_panel("deg_filter", DASHBOARD_FILTERS, neumap_data_base@meta.data),
          actionButton("deg_reset_filters", "Reset All Filters", icon = icon("refresh"), class = "btn-default btn-block"),
          uiOutput("deg_cell_count_summary_ui"), hr(),
          h4("2. Define Comparison Groups"),
          selectInput("deg_comparison_variable", "Compare Groups By:",
            choices = c("Selected Combinations", ANNOTATION_CATEGORIES_RNA),
            selected = "Selected Combinations"
          ),
          hr(),
          h4("3. Configure Analysis"),
          numericInput("n_heatmap_genes", "Top N genes per group:", value = 10, min = 1, max = 50),
          numericInput("n_heatmap_cells", "Max cells per group (for heatmap):", value = 100, min = 10, max = 200),
          sliderInput("log2fc_filter", "Minimum Log2 Fold Change:", min = 0, max = 5, value = 1.5, step = 0.05),
          sliderInput("pval_filter", "Maximum Adjusted P-value:", min = 0, max = 0.1, value = 0.05, step = 0.01),
          checkboxInput("enable_downsampling", "Enable downsampling for speed", value = TRUE),
          conditionalPanel(
            condition = "input.enable_downsampling == true",
            numericInput("max_cells_deg", "Max cells per group (for DEG):", value = 500, min = 50, max = 5000, step = 50)
          ),
          hr(),
          actionButton("run_deg", "Run Differential Analysis", icon = icon("atom"), class = "btn-primary btn-block"),
          create_color_menu_ui("tab2")
        ),
        mainPanel(
          width = 9,
          wellPanel(
            h4("Refine Analysis Population"),
            p("The checkbox selections on the left generate a list of all unique metadata combinations. You can deselect specific combinations here before running the analysis.", style = "font-size: 0.9em; color: #666;"),
            checkboxInput("deg_include_hub_status", "Include 'Hub Status' in Combinations", value = FALSE),
            hr(),
            withSpinner(uiOutput("deg_combination_selector_ui"))
          ),
          wellPanel(
            h5("Heatmap Annotation Controls"),
            fluidRow(
              column(6, checkboxGroupInput("additional_annotations", "Additional Heatmap Annotations:", choices = ANNOTATION_CATEGORIES_RNA, inline = TRUE)),
              column(6, uiOutput("annotation_sorter_ui"))
            )
          ),
          hr(),
          h4("Marker Gene Heatmap"),
          p(style = "font-size: 0.9em; color: #6c757d;", "Scaled expression of top differentially expressed genes for each comparison group."),
          uiOutput("deg_capping_message_ui"),
          withSpinner(jqui_resizable(plotOutput("deg_heatmap", height = "800px"))),
          div(
            style = "text-align: right; margin-top: 5px;",
            downloadButton("download_deg_heatmap_csv", "Download CSV Data", class = "btn-outline-success btn-sm"),
            downloadButton("download_deg_heatmap", "Download PDF", class = "btn-outline-primary btn-sm")
          )
        )
      )
    ),
    tabPanel(
      "Chromatin & RNA Comparison",
      if (!HAS_ATAC_DATA) {
        fluidRow(column(12, h3("Chromatin & RNA Comparison Disabled"), div(class = "alert alert-warning", strong("Missing Prerequisite Data:"), "The required ATAC-seq data ('dogma_FINAL_COMPLETE_v5.rds') was not found.")))
      } else {
        sidebarLayout(
          sidebarPanel(
            width = 3,
            style = "max-height: 85vh; overflow-y: auto;",
            h3("Comparison Controls"), hr(),
            h4("1. Select Gene for Comparison"),
            selectizeInput("comparison_gene_select", "Select Gene:", choices = NULL), hr(),
            # --- Visualization Options sliders REMOVED from sidebar ---
            h4("2. Filter Cell Populations"), # Renumbered
            p("Filter cells for RNA plots:", style = "font-weight: bold;"),
            create_filter_panel("rna_comp_filter", ANNOTATION_CATEGORIES_RNA, neumap_data_base@meta.data),
            hr(),
            h4("3. Configure Analysis Mode (RNA)"), # Renumbered
            checkboxInput("rna_comp_combinatorial", "Enable Combinatorial (AND) Comparison", value = TRUE),
            p(
              style = "font-size:0.8em; color: #6c757d;",
              "Note: RNA filters above control the 'NeuMap RNA UMAP' visualization on the left."
            ),
            hr(style = "border-top: 2px solid #333;"),
            p("Filter cells for ATAC plots:", style = "font-weight: bold;"),
            create_filter_panel("atac_comp_filter", ANNOTATION_CATEGORIES_ATAC, dogma_data@meta.data),
            hr(),
            h4("4. Configure Analysis Mode (ATAC)"), # Renumbered
            checkboxInput("atac_comp_combinatorial", "Enable Combinatorial (AND) Comparison", value = TRUE),
            p(
              style = "font-size:0.8em; color: #6c757d;",
              "Note: ATAC filters above control the 'Promoter Accessibility' visualization on the right."
            ),
            hr(),
            actionButton("apply_comp_filters", "Apply Filters", icon = icon("search"), class = "btn-primary btn-block"),
            actionButton("reset_comp_filters", "Reset Filters", icon = icon("refresh"), class = "btn-default btn-block"),
            create_color_menu_ui("tab3")
          ),
          mainPanel(
            width = 9,
            wellPanel(
              p(strong("Methodology:"), "To enable direct visual comparison, the ATAC-seq cells are computationally projected onto the coordinate system of the NeuMap RNA reference. This projection was performed using Seurat's integration workflow, which maps the ATAC cells into the RNA UMAP space based on their gene accessibility scores.")
            ),
            wellPanel(
              h4("Refine RNA Analysis Population"),
              p("The checkbox selections on the left generate a list of all unique metadata combinations. You can deselect specific combinations here before running the analysis.", style = "font-size: 0.9em; color: #666;"),
              checkboxInput("rna_comp_include_hub_status", "Include 'Hub Status' in Combinations", value = FALSE),
              hr(),
              uiOutput("rna_comp_combination_selector_ui")
            ),
            wellPanel(
              h4("Refine ATAC Analysis Population"),
              p("The checkbox selections on the left generate a list of all unique metadata combinations. You can deselect specific combinations here before running the analysis.", style = "font-size: 0.9em; color: #666;"),
              checkboxInput("atac_comp_include_hub_status", "Include 'Predicted Hub Status' in Combinations", value = FALSE),
              hr(),
              uiOutput("atac_comp_combination_selector_ui")
            ),
            # --- New Structure for Controls Above Plots ---
            fluidRow(
              column(
                6, # RNA Controls
                wellPanel(
                  h5("NeuMap RNA UMAP Controls"),
                  radioButtons("rna_comp_display_mode", "Display Mode:",
                    choices = c("Show all cells" = "all", "Highlight selected" = "highlight"),
                    selected = "highlight", inline = TRUE
                  ),
                  sliderInput("rna_comp_alpha", "Dot Transparency:", min = 0.1, max = 1.0, value = 0.7, step = 0.1),
                  sliderInput("rna_comp_pt_size", "Dot Size:", min = 0.1, max = 3.0, value = 1.0, step = 0.1)
                )
              ),
              column(
                6, # ATAC Controls
                wellPanel(
                  h5("Promoter Accessibility UMAP Controls"),
                  radioButtons("atac_comp_display_mode", "Display Mode:",
                    choices = c("Show all cells" = "all", "Highlight selected" = "highlight"),
                    selected = "highlight", inline = TRUE
                  ),
                  sliderInput("atac_comp_alpha", "Dot Transparency:", min = 0.1, max = 1.0, value = 0.7, step = 0.1),
                  sliderInput("atac_comp_pt_size", "Dot Size:", min = 0.1, max = 3.0, value = 1.0, step = 0.1)
                )
              )
            ),
            # --- End New Structure ---
            fluidRow(
              column(
                6,
                h4(textOutput("rna_plot_title")),
                p(style = "font-size: 0.9em; color: #6c757d;", "Gene expression on the NeuMap RNA reference."),
                withSpinner(jqui_resizable(plotOutput("rna_expression_umap", height = "500px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_rna_comparison_umap", "Download PDF", class = "btn-outline-primary btn-sm"))
              ),
              column(
                6,
                h4(textOutput("atac_plot_title")),
                p(style = "font-size: 0.9em; color: #6c757d;", "Promoter accessibility scores, with ATAC cells projected onto the RNA reference coordinates."),
                withSpinner(jqui_resizable(plotOutput("atac_accessibility_umap", height = "500px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_atac_comparison_umap", "Download PDF", class = "btn-outline-primary btn-sm"))
              )
            ),
            hr(),
            fluidRow(
              column(
                6,
                h4("Expression vs. Accessibility Correlation"),
                p(style = "font-size: 0.9em; color: #6c757d;", "Correlation of mean gene expression (RNA) with mean promoter accessibility (ATAC), aggregated by hub status."),
                withSpinner(jqui_resizable(plotOutput("correlation_plot", height = "400px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_correlation_plot", "Download PDF", class = "btn-outline-primary btn-sm"))
              ),
              column(
                6,
                h4("ATAC Reference UMAP"),
                p(style = "font-size: 0.9em; color: #6c757d;", "Categorical metadata variable for the ATAC dataset."),
                selectInput("atac_reference_variable_tab3", "Reference Variable:", choices = ANNOTATION_CATEGORIES_ATAC, selected = "Predicted Hub Status (ATAC-seq)"),
                withSpinner(jqui_resizable(plotOutput("atac_reference_umap_tab3", height = "400px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_atac_ref_umap_tab3", "Download PDF", class = "btn-outline-primary btn-sm"))
              )
            )
          )
        )
      }
    ),
    tabPanel(
      "TF Motif Activity Analysis",
      if (!HAS_ATAC_DATA) {
        fluidRow(column(12, h3("Motif Activity Analysis Disabled"), div(class = "alert alert-warning", strong("Missing Prerequisite Data:"), "The required 'MOTIF' assay was not found.")))
      } else {
        sidebarLayout(
          sidebarPanel(
            width = 3,
            style = "max-height: 85vh; overflow-y: auto;",
            h3("TF Motif Activity Analysis"),
            p("Analyze differential TF motif activity scores."), hr(),
            h4("1. Define Cell Population"),
            p("Select the total population of cells for analysis.", style = "font-size:0.8em; color: #6c757d;"),
            create_filter_panel("motif_filter", MOTIF_TAB_FILTERS, dogma_data@meta.data),
            uiOutput("motif_cell_count_summary_ui"), hr(),
            h4("2. Configure Analysis Mode"),
            checkboxInput("motif_combinatorial_filter", "Enable Combinatorial (AND) Comparison", value = TRUE),
            conditionalPanel(condition = "input.motif_combinatorial_filter == true", p("Compares unique combinations of your selections against each other (AND logic).", style = "font-size:0.8em; color: #6c757d;")),
            conditionalPanel(condition = "input.motif_combinatorial_filter == false", p("Compares each of your selections as individual groups (OR logic).", style = "font-size:0.8em; color: #6c757d;")),
            hr(),
            h4("3. Analysis Options"),
            h5("Performance Options"),
            checkboxInput("motif_enable_downsampling", "Enable downsampling", value = TRUE),
            conditionalPanel(
              condition = "input.motif_enable_downsampling == true",
              numericInput("motif_max_cells", "Max cells per group:", value = 1000, min = 50, max = 5000, step = 50)
            ),
            selectInput("motif_test_use", "Statistical Test:", choices = c("Wilcoxon (Fast)" = "wilcox", "T-test" = "t"), selected = "wilcox"),
            hr(),
            actionButton("run_motif_analysis", "Run Analysis", icon = icon("atom"), class = "btn-primary btn-block"),
            actionButton("reset_motif_filters", "Reset Selections", icon = icon("refresh"), class = "btn-default btn-block"),
            hr(),
            h4("4. Visualize TF Activity"),
            p("Click a row in the results table to select a TF and auto-load its targets.", style = "font-size:0.8em; color: #6c757d;"),
            selectizeInput("viz_tf_select", "Select TF/Motif:", choices = NULL),
            h5("Threshold for 'Positive' Activity"),
            sliderInput("motif_threshold", NULL, min = 0, max = 1, value = 0.1, step = 0.01),
            create_color_menu_ui("tab4")
          ),
          mainPanel(
            width = 9,
            wellPanel(
              style = "background-color: #f0f8ff;",
              h4("Cross-Modality Analysis Framework"),
              p(style = "font-size: 0.9em;", "This tab integrates ATAC-seq (regulatory potential) and RNA-seq (gene expression). TF activity is inferred from ATAC-seq data using ChromVAR, which calculates a bias-corrected score for motif accessibility that accounts for known technical confounders such as sequencing depth and GC content. Putative TF-target relationships are sourced from the DoRothEA regulon database."),
              p(style = "font-size: 0.9em;", "To create a functional link between the distinct ATAC and RNA modalities, a multi-step contextual mapping is employed. First, ATAC sample identities were spatially resolved within the transcriptomic-defined UMAP landscape using k-Nearest Neighbor (kNN) imputation, effectively mapping chromatin accessibility profiles onto the established RNA topology. Second, a functional Predicted Hub Status was transferred from the RNA reference using a direct nearest-neighbor projection in this shared space. This ensures that the ATAC hub identities spatially mirror the transcriptomic topology, enabling precise downstream analysis.")
            ),
            h4("Identify Key Regulators via Differential Motif Accessibility"),
            wellPanel(h5("Analysis Parameters Summary"), textOutput("motif_analysis_summary_text")),
            conditionalPanel(
              condition = "input.motif_combinatorial_filter == false",
              div(
                class = "alert alert-info",
                p(strong("Individual Comparison Mode:"), " The plots below will show all selected groups (even if they overlap), while this statistical table uses prioritized, mutually exclusive groups for accuracy.")
              )
            ),
            withSpinner(DT::DTOutput("motif_results_table")),
            hr(),
            h4("Visualize Regulator Activity (ATAC-seq)"),
            fluidRow(
              column(
                6,
                h5("Activity Score UMAP"),
                p(style = "font-size: 0.9em; color: #6c757d;", "ChromVAR activity score of the selected TF projected onto the UMAP."),
                withSpinner(jqui_resizable(plotOutput("motif_activity_plot", height = "400px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_motif_activity_umap", "Download PDF", class = "btn-outline-primary btn-sm"))
              ),
              column(
                6,
                h5("ATAC Reference UMAP"),
                p(style = "font-size: 0.9em; color: #6c757d;", "UMAP colored by the selected ATAC metadata for contextual reference."),
                selectInput("atac_reference_variable_tab4", "Reference Variable:", choices = ANNOTATION_CATEGORIES_ATAC, selected = "Predicted Hub Status (ATAC-seq)"),
                withSpinner(jqui_resizable(plotOutput("atac_reference_umap_tab4", height = "400px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_atac_ref_umap_tab4", "Download PDF", class = "btn-outline-primary btn-sm"))
              )
            ),
            hr(),
            h4("Quantify Regulator Activity Across ATAC Subgroups"),
            fluidRow(
              column(12, checkboxInput("show_barchart_labels_tab4", "Show Value Labels on Bar Charts", value = TRUE))
            ),
            fluidRow(
              column(
                6,
                h5("% Positive Cells by Subgroup"),
                p(style = "font-size: 0.9em; color: #6c757d;", "Percentage of cells per subgroup exceeding the activity threshold."),
                withSpinner(jqui_resizable(plotOutput("motif_percent_positive_barchart", height = "350px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_motif_percent_barchart", "Download PDF", class = "btn-outline-primary btn-sm"))
              ),
              column(
                6,
                h5("Mean Activity in Positive Cells"),
                p(style = "font-size: 0.9em; color: #6c757d;", "Average activity score in cells defined as 'positive'."),
                withSpinner(jqui_resizable(plotOutput("motif_activity_barchart", height = "350px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_motif_barchart", "Download PDF", class = "btn-outline-primary btn-sm"))
              )
            ),
            hr(),
            fluidRow(
              column(
                12,
                h5("Activity Score Distribution"),
                p(style = "font-size: 0.9em; color: #6c757d;", "Distribution of ChromVAR activity scores across the selected ATAC subgroups."),
                withSpinner(jqui_resizable(plotOutput("motif_activity_distribution_plot", height = "400px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_motif_dist_plot", "Download PDF", class = "btn-outline-primary btn-sm"))
              )
            ),
            hr(),
            h4("Cross-Validate with TF Gene Expression (in RNA-seq)"),
            p(style = "font-size: 0.9em; color: #6c757d;", "Distribution of the TF's own mRNA expression levels in the corresponding RNA populations. (Cross-validation check for the inferred activity shown above)"),
            fluidRow(
              column(
                12,
                withSpinner(jqui_resizable(plotOutput("tf_expression_violin_plot", height = "400px"))),
                div(style = "text-align: right; margin-top: 5px;", downloadButton("download_tf_expression_violin", "Download PDF", class = "btn-outline-primary btn-sm"))
              )
            ),
            hr(),
            h4("Analyze Expression of Putative Downstream Targets (in RNA-seq)"),
            p(style = "font-size: 0.9em; color: #6c757d;", "This section analyzes the expression of the selected TF's putative target genes (retrieved from the DoRothEA database) within the corresponding RNA cell populations."),
            wellPanel(
              if (is.null(dorothea_regulon_mm)) {
                div(class = "alert alert-danger", "DoRothEA database could not be loaded. This feature is disabled.")
              } else {
                tagList(
                  fluidRow(
                    column(4, selectizeInput("target_analysis_tf_select", "1. Select TF:", choices = NULL)),
                    column(
                      4,
                      selectInput("dorothea_confidence_level", "2. DoRothEA Confidence:",
                        choices = c("A", "B", "C", "D", "E"),
                        selected = c("A", "B", "C", "D", "E"),
                        multiple = TRUE
                      ),
                      actionButton("load_suggested_targets", "3. Load Targets", icon = icon("download"), class = "btn-info btn-sm btn-block")
                    ),
                    column(4, actionButton("analyze_target_expression", "4. Analyze Expression", icon = icon("fire"), class = "btn-success", style = "margin-top: 25px;"))
                  ),
                  p(
                    style = "font-size: 0.8em; color: #6c757d;",
                    strong("Note:"), "Confidence levels (A-E) reflect descending evidence quality, from curated/validated (A) to co-expression-based (E)."
                  ),
                  hr(),
                  h5("Heatmap Display Options"),
                  fluidRow(
                    column(4, checkboxGroupInput("target_heatmap_annotations", "Annotations:", choices = TARGET_HEATMAP_ANNOTATIONS, inline = TRUE)),
                    column(8, uiOutput("target_heatmap_sorter_ui"))
                  ),
                  hr(),
                  h5("Heatmap Filtering & Sizing"),
                  fluidRow(
                    column(4, numericInput("n_target_heatmap_genes", "Top N Genes:", value = 5, min = 1, max = 50)),
                    column(4, sliderInput("log2fc_target_heatmap_filter", "Min. Log2FC:", min = 0, max = 4, value = 0.25, step = 0.1)),
                    column(4, numericInput("max_cells_target_heatmap", "Max Cells / Group:", value = 50, min = 10, max = 200))
                  )
                )
              }
            ),
            wellPanel(
              h4("Refine RNA Population for Heatmap"),
              p("Deselect specific combinations to exclude them from the heatmap analysis.", style = "font-size: 0.9em; color: #666;"),
              checkboxInput("target_heatmap_include_hub_status", "Include 'Hub Status' in Combinations", value = TRUE),
              hr(),
              uiOutput("target_combination_selector_ui")
            ),
            uiOutput("bridge_logic_summary_ui"),
            textAreaInput("target_gene_display", "Loaded Target Genes:", rows = 5, placeholder = "Genes will appear here after clicking 'Load Targets'..."),
            uiOutput("target_capping_message_ui"),
            withSpinner(jqui_resizable(plotOutput("target_gene_heatmap", height = "600px"))),
            div(
              style = "text-align: right; margin-top: 5px;",
              downloadButton("download_target_heatmap_csv", "Download CSV Data", class = "btn-outline-success btn-sm"),
              downloadButton("download_target_heatmap", "Download PDF", class = "btn-outline-primary btn-sm")
            ),
            hr(),
            h4("Exploratory: Find Upstream Regulators from a Gene List"),
            p(style = "font-size: 0.9em; color: #6c757d;", "Perform a reverse analysis: input a list of genes to identify their potential upstream TF regulators from the DoRothEA database."),
            wellPanel(
              textAreaInput("upstream_gene_input", "Enter Target Gene(s):", rows = 4, placeholder = "S100a8\nS100a9\nIl1b..."),
              actionButton("run_upstream_analysis", "Find Regulators", icon = icon("search"), class = "btn-info")
            ),
            withSpinner(DT::DTOutput("upstream_tf_table")),
            hr(),
            wellPanel(
              h5("Regulator Expression Visualization"),
              p(style = "font-size: 0.9em; color: #6c757d;", "Generates RNA expression plots for selected upstream regulators from the table."),
              sliderInput("upstream_threshold", "Threshold for 'Positive' Expression:", min = 0, max = 4, value = 1.0, step = 0.1),
              actionButton("visualize_selected_regulators", "Visualize Selected Regulator(s)", icon = icon("images"), class = "btn-success btn-block")
            ),
            withSpinner(jqui_resizable(uiOutput("upstream_plots_ui"))),
            div(style = "text-align: right; margin-top: 5px;", downloadButton("download_upstream_plot", "Download PDF", class = "btn-outline-primary btn-sm"))
          )
        )
      }
    ),
    tabPanel(
      "Contact",
      h4("Contact Information"),
      hr(),
      p("For questions or support, please contact:"),
      p(tags$a(href = "mailto:maxence.benjamin@gmail.com", "maxence.benjamin@gmail.com"), style = "font-size: 1.2em; font-weight: bold;")
    )
  )
)
