# =============================================================================
# NeuMap Explorer :: server.R
#
# DESCRIPTION:
# This script contains the server-side logic for the application. It defines
# all the reactive expressions, observers, and outputs (plots, tables, etc.)
# that respond to user input.
# =============================================================================

server <- function(input, output, session) {
  rv <- reactiveValues(
    seurat_object = neumap_data_base,
    dogma_object = dogma_data,
    deg_results = NULL,
    calculated_scores = c(PRECOMPUTED_SCORES),
    motif_results_raw = NULL,
    motif_results_processed = NULL,
    motif_analysis_summary = "Analysis has not been run yet.",
    motif_analysis_groups = NULL,
    motif_analysis_group_by = NULL,
    target_gene_list = NULL,
    upstream_tf_results = NULL,
    bridge_logic_message = NULL,
    color_maps = reactiveValues(
      population = list(),
      heatmap_annotations = list(),
      density_palette = "viridis"
    )
  )

  # --- Remap ATAC Hub Status (User Request) ---
  observe({
    req(rv$dogma_object)

    # Check if remapping is needed (Lung + LLC -> IS-II)
    # Using flexible matching to catch variations like "Lung" or "lung", "LLC" or "Lewis Lung Carcinoma"
    meta <- rv$dogma_object@meta.data

    # Identify cells to remap: Tissue contains "Lung" AND Disease contains "LLC" (case-insensitive)
    # Adjust column names if they differ slightly, but based on previous code they seem to be:
    # "Tissue (ATAC-seq)" and "Disease Category (ATAC-seq)"

    tissue_col <- "Tissue (ATAC-seq)"
    disease_col <- "Disease Category (ATAC-seq)"
    hub_col <- "Predicted Hub Status (ATAC-seq)"
    score_col <- "predicted.id.score"

    if (all(c(tissue_col, disease_col, hub_col) %in% colnames(meta))) {
      # Base criteria: Lung + LLC
      base_criteria <- grepl("Lung", meta[[tissue_col]], ignore.case = TRUE) &
        grepl("LLC", meta[[disease_col]], ignore.case = TRUE)

      # Spatial Refinement: Use UMAP_2 coordinate to define the "IS-II Island"
      # Analysis showed a clear gap around UMAP_2 = 2.83
      # Cells ABOVE this cutoff are the true IS-II tip.
      # Cells BELOW are likely IS-I or other states.

      umap_2_vals <- rv$dogma_object@reductions$umap@cell.embeddings[, "UMAP_2"]
      spatial_criteria <- umap_2_vals > 2.83

      # Define the two groups
      is_ii_cells <- which(base_criteria & spatial_criteria)
      is_i_cells <- which(base_criteria & !spatial_criteria)

      current_hubs <- meta[[hub_col]]
      if (is.factor(current_hubs)) {
        levels(current_hubs) <- unique(c(levels(current_hubs), "IS-II", "IS-I"))
      }

      # Apply the spatial override
      if (length(is_ii_cells) > 0) {
        current_hubs[is_ii_cells] <- "IS-II"
      }

      if (length(is_i_cells) > 0) {
        current_hubs[is_i_cells] <- "IS-I"
      }

      # --- DIRECT TRANSCRIPTOMIC TRANSFER (k-NN from RNA) ---
      # Assign each ATAC cell the hub label of its nearest RNA neighbor in UMAP space.
      # This ensures the ATAC hub map perfectly mirrors the RNA topology.

      # 1. Get RNA UMAP coordinates and Hub labels
      rna_coords <- rv$seurat_object@reductions$umap@cell.embeddings
      rna_hubs <- rv$seurat_object@meta.data$hub_status

      # 2. Get ATAC UMAP coordinates
      atac_coords <- rv$dogma_object@reductions$umap@cell.embeddings

      # 3. Find nearest RNA neighbor for each ATAC cell (k=1)
      # This effectively "projects" the RNA labels onto the ATAC cells based on location

      # Filter NAs from training data (knn does not support NAs in cl)
      valid_train_idx <- !is.na(rna_hubs)
      if (sum(valid_train_idx) == 0) {
        warning("No valid RNA hub labels found for KNN training.")
        return()
      }

      transferred_hubs <- class::knn(
        train = rna_coords[valid_train_idx, ],
        test = atac_coords,
        cl = rna_hubs[valid_train_idx],
        k = 1
      )

      # 4. Apply the transferred labels
      current_hubs <- as.character(transferred_hubs)

      # --- SPATIAL OVERRIDE (IS-II Tip) ---
      # Keep the strict definition for the IS-II tip as requested
      if (length(is_ii_cells) > 0) {
        current_hubs[is_ii_cells] <- "IS-II"
      }

      # Ensure factor levels are consistent
      if (is.factor(rv$dogma_object@meta.data[[hub_col]])) {
        current_hubs <- factor(current_hubs, levels = levels(rv$dogma_object@meta.data[[hub_col]]))
      }

      rv$dogma_object@meta.data[[hub_col]] <- current_hubs

      # --- UPDATE UI FILTERS ---
      # Ensure the sidebar filter reflects the new Hub categories (e.g., APC, IS-II, etc.)
      # We need to update the choices for "Predicted Hub Status (ATAC-seq)"

      new_hub_levels <- sort(unique(as.character(current_hubs)))

      # Sanitize the category name to match the input ID format used in create_filter_panel
      # The ID format is "filter_" + sanitized_category_name
      # But wait, the ATAC tab uses "atac_comp_filter_" + sanitized_name
      # And the Motif tab uses "motif_filter_" + sanitized_name
      # We should update BOTH if possible, or at least the one in Tab 3 (ATAC Comparison)

      sanitized_hub_col <- gsub("[^a-zA-Z0-9_]", "_", hub_col)

      # Update Tab 3 Filter
      updateCheckboxGroupInput(session, paste0("atac_comp_filter_", sanitized_hub_col),
        choices = new_hub_levels,
        selected = new_hub_levels,
        inline = TRUE
      )

      # Update Tab 4 Filter (Motif) - assuming it uses the same column
      updateCheckboxGroupInput(session, paste0("motif_filter_", sanitized_hub_col),
        choices = new_hub_levels,
        selected = new_hub_levels,
        inline = TRUE
      )
    }
  })

  target_rna_subset <- reactiveVal(NULL)
  target_deg_results <- reactiveVal(NULL)
  deg_subset_object <- reactiveVal(NULL)
  motif_analysis_subset_object <- reactiveVal(NULL)
  motif_plot_data <- reactiveVal(NULL)

  selections_for_tab1 <- reactiveVal(NULL)

  updateSelectizeInput(session, "explorer_gene_select", choices = ALL_RNA_GENES, selected = DEFAULT_GENE, server = TRUE)
  updateSelectizeInput(session, "explorer_sig_genes", choices = ALL_RNA_GENES, server = TRUE)
  updateSelectizeInput(session, "comparison_gene_select", choices = COMMON_GENES, selected = DEFAULT_GENE, server = TRUE)
  updateSelectizeInput(session, "viz_tf_select", choices = ALL_TF_NAMES, server = TRUE)
  shinyjs::disable("target_gene_display")

  available_meta_main <- reactiveVal(neumap_data_base@meta.data)
  available_meta_deg <- reactiveVal(neumap_data_base@meta.data)
  available_meta_rna_comp <- reactiveVal(neumap_data_base@meta.data)

  create_select_all_observers <- function(prefix, categories_list, metadata) {
    lapply(categories_list, function(cat) {
      local({
        sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
        id <- paste0(prefix, "_", sanitized_cat)
        choices <- get_initial_choices_ui(metadata, cat)
        if (cat %in% c("hub_status", "predicted_hub_status", "Predicted Hub Status (ATAC-seq)")) choices <- choices[choices != "NA"]

        observeEvent(input[[paste0("select_all_", id)]], {
          updateCheckboxGroupInput(session, id, selected = choices)
        })
        observeEvent(input[[paste0("deselect_all_", id)]], {
          updateCheckboxGroupInput(session, id, selected = character(0))
        })
      })
    })
  }

  create_cascading_observers <- function(prefix, categories_list, source_metadata, available_meta_rv) {
    lapply(seq_along(categories_list), function(i) {
      sanitized_current_cat <- gsub("[^a-zA-Z0-9_]", "_", categories_list[i])
      observeEvent(input[[paste0(prefix, sanitized_current_cat)]],
        {
          meta_subset <- source_metadata
          for (j in 1:i) {
            sanitized_j_cat <- gsub("[^a-zA-Z0-9_]", "_", categories_list[j])
            selection <- input[[paste0(prefix, sanitized_j_cat)]]
            if (!is.null(selection) && length(selection) > 0) {
              meta_subset <- meta_subset %>% dplyr::filter(.data[[categories_list[j]]] %in% selection)
            }
          }
          if (i < length(categories_list)) {
            for (k in (i + 1):length(categories_list)) {
              downstream_filter_name_orig <- categories_list[k]
              sanitized_downstream_name <- gsub("[^a-zA-Z0-9_]", "_", downstream_filter_name_orig)
              choices_df <- source_metadata
              for (m in 1:(k - 1)) {
                sanitized_m_cat <- gsub("[^a-zA-Z0-9_]", "_", categories_list[m])
                selection_for_choice <- input[[paste0(prefix, sanitized_m_cat)]]
                if (!is.null(selection_for_choice) && length(selection_for_choice) > 0) {
                  choices_df <- choices_df %>% dplyr::filter(.data[[categories_list[m]]] %in% selection_for_choice)
                }
              }
              current_downstream_selection <- input[[paste0(prefix, sanitized_downstream_name)]]
              new_choices <- get_initial_choices_ui(choices_df, downstream_filter_name_orig)
              if (downstream_filter_name_orig %in% c("hub_status", "predicted_hub_status", "Predicted Hub Status (ATAC-seq)")) new_choices <- new_choices[new_choices != "NA"]
              valid_selection <- intersect(current_downstream_selection, new_choices)
              updateCheckboxGroupInput(session, inputId = paste0(prefix, sanitized_downstream_name), choices = new_choices, selected = valid_selection)
            }
          }
          available_meta_rv(meta_subset)
        },
        ignoreNULL = FALSE,
        ignoreInit = TRUE
      )
    })
  }

  observe_reset_button <- function(button_id, prefix, categories_list) {
    observeEvent(input[[button_id]], {
      for (cat in categories_list) {
        sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
        id <- paste0(prefix, "_", sanitized_cat)
        updateCheckboxGroupInput(session, id, selected = character(0))
      }
    })
  }

  observe({
    create_select_all_observers("filter", DASHBOARD_FILTERS, neumap_data_base@meta.data)
    create_select_all_observers("deg_filter", DASHBOARD_FILTERS, neumap_data_base@meta.data)
    if (HAS_ATAC_DATA) {
      create_select_all_observers("motif_filter", MOTIF_TAB_FILTERS, dogma_data@meta.data)
      create_select_all_observers("atac_comp_filter", ANNOTATION_CATEGORIES_ATAC, dogma_data@meta.data)
    }
    create_select_all_observers("rna_comp_filter", ANNOTATION_CATEGORIES_RNA, neumap_data_base@meta.data)

    create_cascading_observers("filter_", DASHBOARD_FILTERS, neumap_data_base@meta.data, available_meta_main)
    create_cascading_observers("deg_filter_", DASHBOARD_FILTERS, neumap_data_base@meta.data, available_meta_deg)
    create_cascading_observers("rna_comp_filter_", ANNOTATION_CATEGORIES_RNA, neumap_data_base@meta.data, available_meta_rna_comp)

    observe_reset_button("deg_reset_filters", "deg_filter", DASHBOARD_FILTERS)
    observe_reset_button("reset_motif_filters", "motif_filter", MOTIF_TAB_FILTERS)
    observeEvent(input$reset_comp_filters, {
      for (cat in ANNOTATION_CATEGORIES_RNA) {
        sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
        updateCheckboxGroupInput(session, paste0("rna_comp_filter_", sanitized_cat), selected = character(0))
      }
      for (cat in ANNOTATION_CATEGORIES_ATAC) {
        sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
        updateCheckboxGroupInput(session, paste0("atac_comp_filter_", sanitized_cat), selected = character(0))
      }
    })
  })

  apply_capping <- function(df, value_col = NULL) {
    if (!isTRUE(input$enable_capping)) {
      if (is.matrix(df)) {
        return(df)
      }
      if (!is.null(value_col) && value_col %in% names(df)) df$capped_value <- df[[value_col]]
      return(df)
    }

    req(input$capping_mode, input$lower_percentile, input$upper_percentile)

    is_matrix <- is.matrix(df)
    original_df <- df

    if (is_matrix) {
      if (ncol(df) == 0) {
        return(df)
      }
      df <- as.data.frame(as.matrix(df)) %>%
        rownames_to_column("feature_id") %>%
        pivot_longer(-feature_id, names_to = "cell_id", values_to = "value")
      value_col <- "value"
    }

    if (is.null(value_col) || !value_col %in% names(df) || nrow(df) == 0) {
      if (is_matrix) {
        return(original_df)
      } else {
        return(df)
      }
    }

    valid_vals <- df[[value_col]][is.finite(df[[value_col]])]
    if (length(valid_vals) < 2) {
      if (is_matrix) {
        return(original_df)
      }
      df$capped_value <- df[[value_col]]
      return(df)
    }

    positive_vals <- valid_vals[valid_vals > 0]

    use_positive_subset <- (length(positive_vals) > 20 && length(positive_vals) / length(valid_vals) < 0.3)

    vals_for_quantile <- if (use_positive_subset) positive_vals else valid_vals

    bounds <- quantile(vals_for_quantile, probs = c(input$lower_percentile / 100, input$upper_percentile / 100), na.rm = TRUE)

    if (any(is.na(bounds))) {
      if (is_matrix) {
        return(original_df)
      }
      df$capped_value <- df[[value_col]]
      return(df)
    }

    lower_bound <- bounds[1]
    upper_bound <- bounds[2]

    if (upper_bound <= lower_bound) {
      upper_bound <- lower_bound + 1e-9
    }

    if (input$capping_mode == "Hide outlier cells (Trim)") {
      if (is_matrix) {
        cells_to_keep <- df %>%
          dplyr::filter(dplyr::between(!!sym(value_col), lower_bound, upper_bound)) %>%
          pull(cell_id) %>%
          unique()
        return(original_df[, intersect(colnames(original_df), cells_to_keep), drop = FALSE])
      } else {
        df_filtered <- df %>% dplyr::filter(dplyr::between(!!sym(value_col), lower_bound, upper_bound))
        df_filtered$capped_value <- df_filtered[[value_col]]
        return(df_filtered)
      }
    } else { # "Cap color range (Winsorize)"
      if (is_matrix) {
        mat <- original_df
        mat[mat < lower_bound] <- lower_bound
        mat[mat > upper_bound] <- upper_bound
        return(mat)
      } else {
        df$capped_value <- pmax(lower_bound, pmin(upper_bound, df[[value_col]]))
        return(df)
      }
    }
  }

  plot_data_base <- reactive({
    req(rv$seurat_object)
    if (!"umap" %in% names(rv$seurat_object@reductions)) {
      showNotification("UMAP reduction not found in RNA object.", type = "error")
      return(NULL)
    }
    umap_coords <- as.data.frame(rv$seurat_object@reductions$umap@cell.embeddings)
    if (ncol(umap_coords) >= 2) {
      colnames(umap_coords)[1:2] <- c("UMAP_1", "UMAP_2")
    } else {
      showNotification("UMAP reduction has fewer than 2 dimensions.", type = "error")
      return(NULL)
    }
    cbind(umap_coords[, c("UMAP_1", "UMAP_2")], rv$seurat_object@meta.data) %>%
      rownames_to_column("cell_id")
  })

  observe({
    req(plot_data_base())
    if (is.null(selections_for_tab1())) {
      selections_for_tab1(plot_data_base() %>% mutate(plot_group = "All Cells"))
    }
  })

  observeEvent(input$apply_filters, {
    base_df <- plot_data_base()
    req(base_df)

    filter_inputs <- list()
    for (cat in DASHBOARD_FILTERS) {
      sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
      filter_inputs[[cat]] <- input[[paste0("filter_", sanitized_cat)]]
    }
    active_filters <- Filter(function(x) !is.null(x) && length(x) > 0, filter_inputs)

    if (length(active_filters) == 0) {
      final_plot_df <- base_df %>% mutate(plot_group = "All Cells")
    } else {
      if (isTRUE(input$cumulative_filter)) {
        # Stage 1: Broad filtering based on all active sidebar filters
        filtered_df_stage1 <- base_df
        for (cat_name in names(active_filters)) {
          selected_values <- active_filters[[cat_name]]
          filtered_df_stage1 <- filtered_df_stage1 %>% dplyr::filter(.data[[cat_name]] %in% selected_values)
        }

        # Add combination column for refinement
        cols_to_unite <- c("disease_category", "tissue")
        if (isTRUE(input$explorer_include_hub_status)) {
          cols_to_unite <- c(cols_to_unite, "hub_status")
        }
        cols_to_unite <- intersect(cols_to_unite, colnames(filtered_df_stage1))

        filtered_df_with_comb <- filtered_df_stage1 %>%
          tidyr::unite("combination", all_of(cols_to_unite), sep = " & ", remove = FALSE, na.rm = TRUE)

        # Stage 2: Refine by selected combinations
        req(input$explorer_refined_selection)
        cells_to_keep_df <- filtered_df_with_comb %>%
          dplyr::filter(combination %in% input$explorer_refined_selection)

        # Join back to base_df to create final plot df with groups
        if (nrow(cells_to_keep_df) > 0) {
          group_labels <- cells_to_keep_df %>% dplyr::select(cell_id, plot_group = combination)
          final_plot_df <- base_df %>%
            dplyr::left_join(group_labels, by = "cell_id") %>%
            mutate(plot_group = replace_na(plot_group, "Not Selected"))
        } else {
          final_plot_df <- base_df %>% mutate(plot_group = "Not Selected")
        }
      } else { # Individual (OR) mode
        plot_df <- base_df %>% mutate(plot_group = "Not Selected")
        filter_hierarchy <- rev(DASHBOARD_FILTERS)
        active_filter_names <- intersect(filter_hierarchy, names(active_filters))

        for (cat_name in active_filter_names) {
          plot_df <- plot_df %>%
            mutate(plot_group = if_else(.data[[cat_name]] %in% active_filters[[cat_name]],
              as.character(.data[[cat_name]]),
              plot_group
            ))
        }
        final_plot_df <- plot_df
      }
    }
    selections_for_tab1(final_plot_df)
  })

  observeEvent(input$reset_filters, {
    for (cat in DASHBOARD_FILTERS) {
      sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
      updateCheckboxGroupInput(session, paste0("filter_", sanitized_cat), selected = character(0))
    }
    req(plot_data_base())
    selections_for_tab1(plot_data_base() %>% mutate(plot_group = "All Cells"))
  })

  output$filter_mode_explanation_ui <- renderUI({
    text <- if (isTRUE(input$cumulative_filter)) {
      "Combinatorial (AND) Mode: Finds cells that match ALL selected criteria across different categories."
    } else {
      "Individual (OR) Mode: Shows all cells that match ANY selected criterion."
    }
    p(style = "font-size: 0.8em; color: #6c757d; font-style: italic; margin-top: -10px;", text)
  })

  explorer_selected_combinations <- reactive({
    meta_subset <- rv$seurat_object@meta.data

    any_filter_active <- FALSE
    for (cat in DASHBOARD_FILTERS) {
      sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
      selected_vals <- input[[paste0("filter_", sanitized_cat)]]
      if (!is.null(selected_vals) && length(selected_vals) > 0) {
        any_filter_active <- TRUE
        meta_subset <- meta_subset %>% dplyr::filter(.data[[cat]] %in% selected_vals)
      }
    }

    if (!any_filter_active || nrow(meta_subset) == 0) {
      return(character(0))
    }

    cols_to_unite <- c("disease_category", "tissue")
    if (isTRUE(input$explorer_include_hub_status)) {
      cols_to_unite <- c(cols_to_unite, "hub_status")
    }
    cols_to_unite <- intersect(cols_to_unite, colnames(meta_subset))

    if (length(cols_to_unite) == 0) {
      return(character(0))
    }

    meta_subset %>%
      dplyr::distinct(across(all_of(cols_to_unite))) %>%
      tidyr::unite("combination", all_of(cols_to_unite), sep = " & ", remove = TRUE, na.rm = TRUE) %>%
      pull(combination) %>%
      sort()
  })

  output$explorer_combination_selector_ui <- renderUI({
    combinations <- explorer_selected_combinations()
    if (length(combinations) == 0) {
      return(p("Select filters above to generate population combinations.", style = "color: #666;"))
    }

    tagList(
      fluidRow(
        column(6, actionLink("select_all_explorer_combinations", "Select All", style = "font-size: 9pt;")),
        column(6, actionLink("deselect_all_explorer_combinations", "Deselect All", style = "font-size: 9pt;"))
      ),
      checkboxGroupInput(
        inputId = "explorer_refined_selection",
        label = NULL,
        choices = combinations,
        selected = combinations,
        inline = TRUE
      )
    )
  })

  observeEvent(input$select_all_explorer_combinations, {
    updateCheckboxGroupInput(session, "explorer_refined_selection", selected = explorer_selected_combinations())
  })
  observeEvent(input$deselect_all_explorer_combinations, {
    updateCheckboxGroupInput(session, "explorer_refined_selection", selected = character(0))
  })

  output$cell_count_summary_ui <- renderUI({
    shiny::tags$p(
      style = "font-size: 0.9em; text-align: center; margin-top: 10px;",
      strong("Available in UI (Live):"),
      em(paste(scales::comma(nrow(available_meta_main())), "cells"))
    )
  })

  active_feature <- reactive({
    if (input$analysis_type == "Gene Expression") {
      req(input$explorer_gene_select)
      return(input$explorer_gene_select)
    } else {
      req(input$explorer_score_select)
      return(input$explorer_score_select)
    }
  })

  feature_expression_data <- reactive({
    req(input$analysis_type)
    feature <- active_feature()
    req(feature)
    base_df <- plot_data_base()
    req(base_df)

    feature_data <- tryCatch(
      {
        if (input$analysis_type == "Gene Expression") {
          if (!feature %in% rownames(rv$seurat_object[["RNA"]])) {
            stop(paste("Gene '", feature, "' not found in RNA assay rownames."))
          }
          expr_values <- GetAssayData(rv$seurat_object, assay = "RNA", layer = "data")[feature, ]
          data.frame(
            cell_id = names(expr_values),
            expression_value = as.numeric(expr_values),
            stringsAsFactors = FALSE
          )
        } else {
          if (!feature %in% colnames(rv$seurat_object@meta.data)) {
            stop(paste("Signature score '", feature, "' not found in metadata."))
          }
          FetchData(rv$seurat_object, vars = feature) %>%
            rownames_to_column("cell_id") %>%
            setNames(c("cell_id", "expression_value"))
        }
      },
      error = function(e) {
        showNotification(paste("Error fetching data for", feature, ":", e$message), type = "error")
        return(NULL)
      }
    )

    req(feature_data)

    plot_df <- base_df %>% left_join(feature_data, by = "cell_id")

    if (!"expression_value" %in% names(plot_df) || !is.numeric(plot_df$expression_value)) {
      showNotification("Could not create valid plot data.", type = "error")
      return(NULL)
    }
    plot_df$is_positive <- plot_df$expression_value > input$explorer_threshold

    return(plot_df)
  })

  summary_stats_data <- reactive({
    primary_selection <- selections_for_tab1()
    req(primary_selection)
    expression_data <- feature_expression_data()
    req(expression_data)

    plot_df <- expression_data %>%
      inner_join(primary_selection %>% dplyr::select(cell_id, plot_group), by = "cell_id") %>%
      dplyr::filter(plot_group != "Not Selected" & plot_group != "All Cells")

    req(nrow(plot_df) > 0)

    data_for_stats <- apply_capping(plot_df, "expression_value")

    if (nrow(data_for_stats) == 0) {
      return(NULL)
    }

    percent_positive <- data_for_stats %>%
      group_by(plot_group) %>%
      summarise(percent_positive = (sum(is_positive) / n()) * 100, .groups = "drop")

    mean_expression <- data_for_stats %>%
      dplyr::filter(is_positive) %>%
      group_by(plot_group) %>%
      summarise(mean_expr = mean(expression_value, na.rm = TRUE), .groups = "drop")

    return(list(percent_positive = percent_positive, mean_expression = mean_expression))
  })

  tab1_color_map <- reactive({
    plot_data <- selections_for_tab1()
    req(plot_data)

    group_levels <- plot_data %>%
      dplyr::filter(plot_group != "Not Selected" & plot_group != "All Cells") %>%
      pull(plot_group) %>%
      unique() %>%
      sort()

    if (length(group_levels) == 0) {
      return(NULL)
    }

    default_palette <- generate_palette(length(group_levels))
    default_color_map <- setNames(default_palette, group_levels)

    # Use global HUB_STATUS_COLORS if applicable
    if (all(group_levels %in% names(HUB_STATUS_COLORS))) {
      default_color_map <- HUB_STATUS_COLORS[group_levels]
    }

    custom_colors <- rv$color_maps$population

    final_color_map <- utils::modifyList(as.list(default_color_map), custom_colors)

    return(final_color_map)
  })

  featured_expression_umap_plot <- reactive({
    feature <- active_feature()
    req(feature)
    plot_df_full <- feature_expression_data()
    req(plot_df_full)
    base_df <- plot_data_base()
    req(base_df)
    selection_data <- selections_for_tab1()
    req(selection_data)

    p_base <- ggplot() +
      labs(title = paste("Expression of", feature)) +
      theme_minimal(base_size = 16) +
      theme(panel.grid = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.position = "right", plot.title = element_text(hjust = 0.5))

    if (input$featured_umap_mode == "highlight") {
      selected_cell_ids <- selection_data %>%
        dplyr::filter(plot_group != "Not Selected" & plot_group != "All Cells") %>%
        pull(cell_id)
      p_base <- p_base + geom_point(data = base_df, aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5, alpha = 0.5)
      if (length(selected_cell_ids) == 0) {
        return(p_base + annotate("text", x = mean(base_df$UMAP_1), y = mean(base_df$UMAP_2), label = "No cells match current selections.", size = 5))
      }
      plot_df <- plot_df_full %>% dplyr::filter(cell_id %in% selected_cell_ids)
    } else {
      plot_df <- plot_df_full
    }

    plot_df_capped <- apply_capping(plot_df, "expression_value")

    # --- Dynamic Color Scale Logic ---
    palette_name <- rv$color_maps$density_palette
    viridis_palettes <- c("viridis", "magma", "plasma", "inferno", "cividis")

    if (palette_name %in% viridis_palettes) {
      color_scale <- scale_color_viridis_c(name = str_to_title(gsub("_", " ", feature)), option = palette_name)
    } else {
      color_scale <- scale_color_distiller(name = str_to_title(gsub("_", " ", feature)), palette = palette_name)
    }
    # --- End Dynamic Logic ---

    p <- p_base +
      geom_point(
        data = plot_df_capped %>% arrange(capped_value),
        aes(x = UMAP_1, y = UMAP_2, color = capped_value),
        size = input$featured_pt_size, alpha = input$featured_alpha
      ) +
      color_scale

    return(p)
  })

  reference_umap_plot <- reactive({
    req(input$reference_umap_variable)
    ref_var <- input$reference_umap_variable
    plot_df_full <- plot_data_base()
    req(plot_df_full)
    base_df <- plot_data_base()
    req(base_df)
    selection_data <- selections_for_tab1()
    req(selection_data)

    p_base <- ggplot() +
      theme_void(base_size = 14) +
      labs(title = str_to_title(gsub("_", " ", ref_var))) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    if (input$reference_umap_mode == "highlight") {
      selected_cell_ids <- selection_data %>%
        dplyr::filter(plot_group != "Not Selected" & plot_group != "All Cells") %>%
        pull(cell_id)
      p_base <- p_base + geom_point(data = base_df, aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5, alpha = 0.5)
      if (length(selected_cell_ids) == 0) {
        return(p_base + annotate("text", x = mean(base_df$UMAP_1), y = mean(base_df$UMAP_2), label = "No cells match current selections.", size = 5))
      }
      plot_df <- plot_df_full %>% dplyr::filter(cell_id %in% selected_cell_ids)
    } else {
      plot_df <- plot_df_full
    }

    if (nrow(plot_df) > 1) plot_df <- plot_df[sample(nrow(plot_df)), ]
    plot_df[[ref_var]] <- droplevels(as.factor(plot_df[[ref_var]]))

    p <- p_base +
      geom_point(data = plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[ref_var]]), size = input$reference_pt_size, alpha = input$reference_alpha)

    if (ref_var == "hub_status") {
      plot_levels <- levels(plot_df[[ref_var]])
      levels_to_color <- plot_levels[plot_levels != "NA"]
      palette <- generate_palette(length(levels_to_color))
      color_map <- setNames(palette, levels_to_color)
      final_color_map <- c(color_map, "NA" = "grey80")

      centroids <- plot_df %>%
        dplyr::filter(.data[[ref_var]] != "NA") %>%
        group_by(.data[[ref_var]]) %>%
        summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")

      p <- p +
        ggrepel::geom_label_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = .data[[ref_var]], fill = .data[[ref_var]]), fontface = "bold", color = "white", size = 4) +
        scale_color_manual(values = final_color_map, breaks = levels_to_color) +
        scale_fill_manual(values = final_color_map) +
        guides(color = "none", fill = "none")
    } else {
      centroids <- plot_df %>%
        group_by(.data[[ref_var]]) %>%
        summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")
      p <- p +
        ggrepel::geom_label_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = .data[[ref_var]], fill = .data[[ref_var]]), fontface = "bold", color = "white", size = 4) +
        guides(color = "none", fill = "none")
    }

    return(p)
  })

  subgroup_location_umap_plot <- reactive({
    plot_df <- selections_for_tab1()
    req(plot_df)
    base_df <- plot_data_base()
    req(base_df)
    color_map <- tab1_color_map()
    req(color_map)

    highlighted_df <- plot_df %>% dplyr::filter(plot_group != "Not Selected" & plot_group != "All Cells")

    p <- ggplot(base_df, aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(color = "grey85", size = 1.0, alpha = 0.5) +
      theme_void(base_size = 14) +
      theme(legend.position = "none")

    if (nrow(highlighted_df) > 0) {
      centroids <- highlighted_df %>%
        group_by(plot_group) %>%
        summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")

      p <- p + geom_point(data = highlighted_df, aes(color = plot_group), size = 1.5, alpha = 0.5) +
        ggrepel::geom_label_repel(data = centroids, aes(label = plot_group, fill = plot_group), color = "white", fontface = "bold") +
        scale_color_manual(values = color_map) +
        scale_fill_manual(values = color_map)
    }
    return(p)
  })

  percent_positive_barchart_plot <- reactive({
    plot_data <- summary_stats_data()
    req(plot_data, !is.null(plot_data$percent_positive), nrow(plot_data$percent_positive) > 0)
    color_map <- tab1_color_map()
    req(color_map)

    p <- ggplot(plot_data$percent_positive, aes(x = plot_group, y = percent_positive, fill = plot_group)) +
      geom_col(show.legend = FALSE)

    if (isTRUE(input$show_barchart_labels_tab1)) {
      p <- p + geom_text(aes(label = paste0(round(percent_positive, 1), "%")), vjust = -0.5)
    }

    p + scale_fill_manual(values = color_map) +
      scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.15))) + # Increased expansion for labels
      labs(x = NULL, y = "Percent (%)") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  })

  mean_expression_barchart_plot <- reactive({
    plot_data <- summary_stats_data()
    req(plot_data, !is.null(plot_data$mean_expression), nrow(plot_data$mean_expression) > 0)
    color_map <- tab1_color_map()
    req(color_map)

    p <- ggplot(plot_data$mean_expression, aes(x = plot_group, y = mean_expr, fill = plot_group)) +
      geom_col(show.legend = FALSE)

    if (isTRUE(input$show_barchart_labels_tab1)) {
      p <- p + geom_text(aes(label = round(mean_expr, 2)), vjust = -0.5)
    }

    p + scale_fill_manual(values = color_map) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + # Increased expansion for labels
      labs(x = NULL, y = "Mean Expression") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  })

  expression_distribution_violin_plot <- reactive({
    primary_selection <- selections_for_tab1()
    req(primary_selection)
    expression_data <- feature_expression_data()
    req(expression_data)
    feature <- active_feature()
    color_map <- tab1_color_map()
    req(color_map)

    plot_df <- expression_data %>%
      inner_join(primary_selection %>% dplyr::select(cell_id, plot_group), by = "cell_id") %>%
      dplyr::filter(plot_group != "Not Selected" & plot_group != "All Cells")

    req(nrow(plot_df) > 0)

    plot_df_capped <- apply_capping(plot_df, "expression_value")

    ggplot(plot_df_capped, aes(x = plot_group, y = expression_value, fill = plot_group)) +
      geom_violin(trim = FALSE, alpha = 0.8, show.legend = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, show.legend = FALSE) +
      scale_fill_manual(values = color_map) +
      labs(
        title = paste("Distribution of", feature),
        x = "Group",
        y = "Expression Level"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()
      )
  })

  output$featured_expression_umap <- renderPlot({
    featured_expression_umap_plot()
  })
  output$reference_umap <- renderPlot({
    reference_umap_plot()
  })
  output$subgroup_location_umap <- renderPlot({
    subgroup_location_umap_plot()
  })
  output$percent_positive_barchart <- renderPlot({
    percent_positive_barchart_plot()
  })
  output$mean_expression_barchart <- renderPlot({
    mean_expression_barchart_plot()
  })
  output$expression_distribution_violin <- renderPlot({
    expression_distribution_violin_plot()
  })

  create_download_handler <- function(plot_id, plot_reactive, filename_prefix, feature_name_reactive = active_feature) {
    downloadHandler(
      filename = function() {
        feature_name <- tryCatch(feature_name_reactive(), error = function(e) "plot")
        clean_feature_name <- gsub("[^a-zA-Z0-9_.-]", "_", paste(feature_name, collapse = "_"))
        paste0(filename_prefix, "_", clean_feature_name, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        plot_size <- input[[paste0(plot_id, "_size")]]
        width_px <- if (!is.null(plot_size$width)) as.numeric(plot_size$width) else 800
        height_px <- if (!is.null(plot_size$height)) as.numeric(plot_size$height) else 600
        if (plot_id == "upstream_composite_plot") {
          num_tfs <- length(visualize_selected_tfs())
          if (num_tfs > 0) height_px <- max(600, 500 * num_tfs)
        }

        width_in <- width_px / 96
        height_in <- height_px / 96
        p <- plot_reactive()
        ggsave(file, plot = p, device = "pdf", width = width_in, height = height_in, units = "in", dpi = 300, limitsize = FALSE)
      }
    )
  }

  output$download_featured_umap <- create_download_handler("featured_expression_umap", featured_expression_umap_plot, "Featured_UMAP")
  output$download_reference_umap <- create_download_handler("reference_umap", reference_umap_plot, "Reference_UMAP", feature_name_reactive = reactive(input$reference_umap_variable))
  output$download_subgroup_umap <- create_download_handler("subgroup_location_umap", subgroup_location_umap_plot, "Subgroup_Location_UMAP", feature_name_reactive = reactive("Subgroups"))
  output$download_percent_barchart <- create_download_handler("percent_positive_barchart", percent_positive_barchart_plot, "Percent_Positive")
  output$download_mean_barchart <- create_download_handler("mean_expression_barchart", mean_expression_barchart_plot, "Mean_Expression")
  output$download_violin_plot_tab1 <- create_download_handler("expression_distribution_violin", expression_distribution_violin_plot, "Expression_Distribution")

  output$signature_score_ui_main <- renderUI({
    tagList(
      selectInput("explorer_score_select", "Select Score:", choices = rv$calculated_scores),
      wellPanel(
        h5("Calculate a New Score"),
        textInput("explorer_sig_name", "Name:", placeholder = "e.g., Inflammation"),
        selectizeInput("explorer_sig_genes", "Select Genes:", choices = NULL, multiple = TRUE),
        textAreaInput("explorer_sig_genes_paste", "Or Paste Genes (Column/CSV):", rows = 3, placeholder = "GeneA\nGeneB, GeneC"),
        actionButton("explorer_run_sig_score", "Save & Use Score", icon = icon("save"), class = "btn-success btn-xs btn-block")
      )
    )
  })

  observeEvent(input$explorer_run_sig_score, {
    req(input$explorer_sig_name)

    # Combine genes from selectize and text area
    genes_from_select <- input$explorer_sig_genes
    genes_from_paste <- if (!is.null(input$explorer_sig_genes_paste) && nzchar(input$explorer_sig_genes_paste)) {
      unlist(strsplit(input$explorer_sig_genes_paste, "[\\s,]+", perl = TRUE))
    } else {
      character(0)
    }

    # ongoing cleaning
    all_genes <- unique(c(genes_from_select, genes_from_paste))
    all_genes <- all_genes[all_genes != ""]

    # Filter valid genes
    valid_genes <- intersect(all_genes, rownames(rv$seurat_object[["RNA"]]))

    sig_name <- make.names(input$explorer_sig_name)

    if (sig_name %in% names(rv$seurat_object@meta.data)) {
      showNotification("Score name already exists.", type = "error")
      return()
    }

    ignored_genes <- setdiff(all_genes, valid_genes)

    if (length(valid_genes) < 2) {
      msg <- paste("Error: At least two VALID genes are required. Found:", length(valid_genes))
      if (length(ignored_genes) > 0) {
        msg <- paste0(msg, "\nIgnored (not found): ", paste(head(ignored_genes, 5), collapse = ", "), ifelse(length(ignored_genes) > 5, "...", ""))
      }
      showNotification(msg, type = "error")
      return()
    }

    # Notify about ignored genes but proceed
    if (length(ignored_genes) > 0) {
      warning_msg <- paste(
        "Warning: The following", length(ignored_genes), "genes were not found and ignored:",
        paste(head(ignored_genes, 10), collapse = ", "), ifelse(length(ignored_genes) > 10, "...", "")
      )
      showNotification(warning_msg, type = "warning", duration = 8)
    }

    # Notify regarding valid genes
    success_msg <- paste("Calculating score with", length(valid_genes), "valid genes.")
    showNotification(success_msg, type = "message", duration = 4)

    withProgress(message = "Calculating Signature Score", value = 0.5, {
      tryCatch(
        {
          rv$seurat_object <- AddModuleScore(rv$seurat_object, features = list(valid_genes), name = sig_name, assay = "RNA")
          names(rv$seurat_object@meta.data)[names(rv$seurat_object@meta.data) == paste0(sig_name, 1)] <- sig_name
          rv$calculated_scores <- c(rv$calculated_scores, sig_name)
          updateSelectInput(session, "explorer_score_select", choices = rv$calculated_scores, selected = sig_name)
          updateSelectizeInput(session, "explorer_sig_genes", selected = "")
          updateTextAreaInput(session, "explorer_sig_genes_paste", value = "")
        },
        error = function(e) {
          showNotification(paste("Error during score calculation:", e$message), type = "error")
        }
      )
    })
  })

  output$deg_cell_count_summary_ui <- renderUI({
    shiny::tags$p(
      style = "font-size: 0.9em; text-align: center; margin-top: 10px; padding: 5px; border: 1px solid #ddd;",
      strong("Available in UI (Live):"),
      em(paste(scales::comma(nrow(available_meta_deg())), "cells"))
    )
  })

  deg_selected_combinations <- reactive({
    meta_subset <- rv$seurat_object@meta.data

    any_filter_active <- FALSE
    for (cat in DASHBOARD_FILTERS) {
      sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
      selected_vals <- input[[paste0("deg_filter_", sanitized_cat)]]
      if (!is.null(selected_vals) && length(selected_vals) > 0) {
        any_filter_active <- TRUE
        meta_subset <- meta_subset %>% dplyr::filter(.data[[cat]] %in% selected_vals)
      }
    }

    if (!any_filter_active || nrow(meta_subset) == 0) {
      return(character(0))
    }

    cols_to_unite <- DEG_COMBINATION_COLS
    if (!isTRUE(input$deg_include_hub_status)) {
      cols_to_unite <- setdiff(cols_to_unite, "hub_status")
    }
    cols_to_unite <- intersect(cols_to_unite, colnames(meta_subset))

    if (length(cols_to_unite) == 0) {
      return(character(0))
    }

    meta_subset %>%
      dplyr::distinct(across(all_of(cols_to_unite))) %>%
      tidyr::unite("combination", all_of(cols_to_unite), sep = "_", remove = TRUE, na.rm = TRUE) %>%
      pull(combination) %>%
      sort()
  })

  output$deg_combination_selector_ui <- renderUI({
    combinations <- deg_selected_combinations()
    if (length(combinations) == 0) {
      return(p("Select filters on the left to generate population combinations.", style = "color: #666;"))
    }

    tagList(
      fluidRow(
        column(6, actionLink("select_all_deg_combinations", "Select All", style = "font-size: 9pt;")),
        column(6, actionLink("deselect_all_deg_combinations", "Deselect All", style = "font-size: 9pt;"))
      ),
      checkboxGroupInput(
        inputId = "deg_refined_selection",
        label = "Select combinations to include in DEG analysis:",
        choices = combinations,
        selected = combinations
      )
    )
  })

  observeEvent(input$select_all_deg_combinations, {
    updateCheckboxGroupInput(session, "deg_refined_selection", selected = deg_selected_combinations())
  })
  observeEvent(input$deselect_all_deg_combinations, {
    updateCheckboxGroupInput(session, "deg_refined_selection", selected = character(0))
  })

  observeEvent(input$run_deg,
    {
      withProgress(message = "Finding Markers...", value = 0.1, {
        req(input$deg_refined_selection)

        setProgress(0.2, detail = "Resolving refined cell selection...")

        cols_to_unite <- DEG_COMBINATION_COLS
        if (!isTRUE(input$deg_include_hub_status)) {
          cols_to_unite <- setdiff(cols_to_unite, "hub_status")
        }
        cols_to_unite <- intersect(cols_to_unite, colnames(rv$seurat_object@meta.data))

        meta_with_combinations <- rv$seurat_object@meta.data %>%
          tidyr::unite("combination", all_of(cols_to_unite), sep = "_", remove = FALSE, na.rm = TRUE)

        filtered_meta <- meta_with_combinations %>%
          dplyr::filter(combination %in% input$deg_refined_selection)

        if (nrow(filtered_meta) == 0) {
          showNotification("CRITICAL: No cells match the refined selection of combinations.", type = "error", duration = 10)
          return()
        }

        if (nrow(filtered_meta) < 20) {
          showNotification("Fewer than 20 cells in the selected population. Please select more cells.", type = "error", duration = 8)
          return()
        }

        setProgress(0.3, detail = "Subsetting data...")
        deg_obj <- subset(rv$seurat_object, cells = rownames(filtered_meta))

        deg_obj <- AddMetaData(deg_obj, metadata = filtered_meta[Cells(deg_obj), "combination", drop = FALSE], col.name = "deg_combination")

        comp_var <- req(input$deg_comparison_variable)

        if (comp_var == "Selected Combinations") {
          Idents(deg_obj) <- "deg_combination"
        } else {
          Idents(deg_obj) <- comp_var
        }

        groups_to_compare <- levels(Idents(deg_obj))
        if (length(groups_to_compare) < 2) {
          comparison_text <- if (comp_var == "Selected Combinations") {
            "selected combinations"
          } else {
            paste0("groups for the comparison category '", comp_var, "'")
          }

          error_message <- paste0(
            "Analysis failed: The final cell population has fewer than two ",
            comparison_text,
            ". This can happen after using the 'Refine Analysis Population' panel. ",
            "Please adjust your selections to include at least two distinct groups or choose a different comparison variable."
          )

          showNotification(ui = error_message, type = "error", duration = 15)
          return()
        }

        deg_subset_object(deg_obj)

        setProgress(0.6, detail = paste("Comparing", length(groups_to_compare), "groups..."))

        deg_params <- list(object = deg_obj, only.pos = FALSE, logfc.threshold = 0.1, assay = "RNA")
        if (isTRUE(input$enable_downsampling)) {
          deg_params$max.cells.per.ident <- input$max_cells_deg
        }

        # Safety check: Filter out groups with too few cells (< 3) to prevent Seurat crash
        cells_per_group <- table(Idents(deg_obj))
        valid_groups <- names(cells_per_group[cells_per_group >= 3])

        if (length(valid_groups) < 2) {
          showNotification("Analysis failed: Fewer than 2 groups have enough cells (>= 3) for comparison.", type = "error")
          return()
        }

        if (length(valid_groups) < length(cells_per_group)) {
          deg_obj <- subset(deg_obj, idents = valid_groups)
          deg_params$object <- deg_obj
          showNotification(paste("Excluded", length(cells_per_group) - length(valid_groups), "small groups (< 3 cells) from analysis."), type = "warning")
        }

        deg_data <- tryCatch(
          {
            do.call(FindAllMarkers, deg_params)
          },
          error = function(e) {
            showNotification(paste("Error during DEG analysis:", e$message), type = "error")
            return(NULL)
          }
        )

        rv$deg_results <- deg_data

        if (is.null(rv$deg_results) || nrow(rv$deg_results) == 0) {
          showNotification("No significant markers found with the current settings.", type = "warning")
          return()
        }

        setProgress(1, detail = "Done!")
        showNotification("Differential gene analysis complete.", type = "message")
      })
    },
    ignoreInit = TRUE
  )

  observeEvent(input$deg_comparison_variable,
    {
      req(input$deg_comparison_variable)

      if (input$deg_comparison_variable != "Selected Combinations") {
        current_annots <- input$additional_annotations
        new_selection <- union(current_annots, input$deg_comparison_variable)
        updateCheckboxGroupInput(session, "additional_annotations", choices = ANNOTATION_CATEGORIES_RNA, selected = new_selection, inline = TRUE)
      }
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
  )

  output$annotation_sorter_ui <- renderUI({
    req(input$deg_comparison_variable)

    cols_for_sorter <- if (input$deg_comparison_variable == "Selected Combinations") {
      unique(c("Combination Group", input$additional_annotations))
    } else {
      unique(c(input$deg_comparison_variable, input$additional_annotations))
    }

    if (length(cols_for_sorter) > 0) {
      rank_list(
        text = "Set Annotation Priority (Top = Closest to Heatmap):",
        labels = cols_for_sorter,
        input_id = "annotation_order"
      )
    }
  })

  filtered_deg <- reactive({
    req(rv$deg_results)
    rv$deg_results %>%
      dplyr::filter(abs(avg_log2FC) >= input$log2fc_filter & p_val_adj <= input$pval_filter)
  })

  deg_heatmap_args <- reactive({
    req(input$enable_capping)
    if (isTRUE(input$enable_capping)) req(input$capping_mode)
    req(filtered_deg(), deg_subset_object())

    deg_obj <- deg_subset_object()

    top_genes_df <- filtered_deg() %>%
      group_by(cluster) %>%
      slice_max(order_by = abs(avg_log2FC), n = input$n_heatmap_genes, with_ties = FALSE) %>%
      ungroup()

    unique_genes_to_plot <- unique(top_genes_df$gene)
    if (length(unique_genes_to_plot) == 0) {
      return("No genes pass the current filters.")
    }

    cells_to_plot <- tryCatch(
      {
        Seurat::WhichCells(deg_obj, idents = levels(Idents(deg_obj)), downsample = input$n_heatmap_cells)
      },
      error = function(e) {
        return(Cells(deg_obj))
      }
    )

    sorting_priority <- input$annotation_order

    if (!is.null(sorting_priority)) {
      sorting_priority[sorting_priority == "Combination Group"] <- "deg_combination"
    }

    main_annot <- if (input$deg_comparison_variable == "Selected Combinations") "deg_combination" else input$deg_comparison_variable

    cell_sort_order <- unique(c(sorting_priority, main_annot))
    cell_sort_order <- intersect(cell_sort_order, colnames(deg_obj@meta.data))


    if (length(cell_sort_order) == 0) {
      annotation_df <- NULL
      ordered_cells <- cells_to_plot
    } else {
      full_annotation_df <- deg_obj@meta.data[cells_to_plot, cell_sort_order, drop = FALSE]
      ordered_annotation_df <- full_annotation_df %>% arrange(across(all_of(cell_sort_order)))
      ordered_cells <- rownames(ordered_annotation_df)

      annotation_df <- ordered_annotation_df[, rev(cell_sort_order), drop = FALSE]
    }

    data_for_heatmap <- GetAssayData(deg_obj, assay = "RNA", layer = "data")[unique_genes_to_plot, ordered_cells, drop = FALSE]

    data_for_heatmap <- apply_capping(data_for_heatmap, NULL)

    if (ncol(data_for_heatmap) < 2) {
      return("Fewer than 2 cells remain after trimming.")
    }

    scaled_data <- t(scale(t(as.matrix(data_for_heatmap))))

    annotation_df <- annotation_df[colnames(scaled_data), , drop = FALSE]

    cap_value <- quantile(abs(scaled_data), 0.95, na.rm = TRUE)
    if (!is.finite(cap_value) || cap_value < 2.5) cap_value <- 2.5
    scaled_data[is.na(scaled_data)] <- 0
    scaled_data[scaled_data > cap_value] <- cap_value
    scaled_data[scaled_data < -cap_value] <- -cap_value

    return(list(mat = scaled_data, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = (nrow(scaled_data) > 1), annotation_col = annotation_df))
  })

  output$deg_heatmap <- renderPlot({
    args <- deg_heatmap_args()
    if (is.character(args)) {
      return(ggplot() +
        annotate("text", x = 1, y = 1, label = args) +
        theme_void())
    }

    # Build custom color list for pheatmap
    annot_colors <- list()
    if (!is.null(args$annotation_col) && is.data.frame(args$annotation_col)) {
      for (col_name in names(args$annotation_col)) {
        levels <- levels(factor(args$annotation_col[[col_name]]))
        default_palette <- setNames(generate_palette(length(levels)), levels)

        custom_palette <- rv$color_maps$heatmap_annotations[[col_name]]
        if (!is.null(custom_palette)) {
          final_palette <- unlist(utils::modifyList(as.list(default_palette), custom_palette))
        } else {
          final_palette <- default_palette
        }
        annot_colors[[col_name]] <- final_palette
      }
    }

    final_args <- c(args, list(annotation_colors = annot_colors))
    do.call(pheatmap::pheatmap, final_args)
  })

  capping_message_dynamic_ui <- reactive({
    req(isTRUE(input$enable_capping))

    mode_text <- if (grepl("Trim", input$capping_mode, fixed = TRUE)) {
      "Trim"
    } else {
      "Winsorize"
    }

    message_str <- sprintf(
      "Data Capping Applied (Mode: %s). Lower: %s%%, Upper: %s%%. (See Tab 1 to change)",
      mode_text,
      input$lower_percentile,
      input$upper_percentile
    )

    p(style = "font-size: 0.9em; color: #666; font-style: italic;", message_str)
  })

  output$deg_capping_message_ui <- renderUI({
    capping_message_dynamic_ui()
  })

  output$target_capping_message_ui <- renderUI({
    capping_message_dynamic_ui()
  })

  output$download_deg_heatmap_csv <- downloadHandler(
    filename = function() {
      paste0("DEG_Heatmap_Data_", input$deg_comparison_variable, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      args <- deg_heatmap_args()
      if (!is.list(args) || is.null(args$mat)) {
        return()
      }
      final_data <- as.data.frame(as.matrix(args$mat)) %>% rownames_to_column("gene")
      write.csv(final_data, file, row.names = FALSE)
    }
  )

  output$download_deg_heatmap <- downloadHandler(
    filename = function() {
      paste0("DEG_Heatmap_", input$deg_comparison_variable, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      args <- deg_heatmap_args()
      if (!is.list(args) || is.null(args$mat)) {
        return()
      }
      plot_size <- input$deg_heatmap_size
      width_in <- if (!is.null(plot_size$width)) as.numeric(plot_size$width) / 96 else 8
      height_in <- if (!is.null(plot_size$height)) as.numeric(plot_size$height) / 96 else 8
      pdf(file, width = width_in, height = height_in)
      do.call(pheatmap::pheatmap, args)
      dev.off()
    }
  )

  selected_comparison_gene <- reactive({
    req(input$comparison_gene_select, nchar(input$comparison_gene_select) > 0)
    input$comparison_gene_select
  })

  filtered_comparison_cells <- eventReactive(input$apply_comp_filters,
    {
      # --- RNA Data Filtering ---
      rna_meta_subset <- rv$seurat_object@meta.data
      rna_cell_ids <- rownames(rna_meta_subset) # Default to all cells

      # Stage 1: Broad filtering based on sidebar selections
      rna_filter_inputs <- list()
      for (cat in ANNOTATION_CATEGORIES_RNA) {
        sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
        selection <- input[[paste0("rna_comp_filter_", sanitized_cat)]]
        if (!is.null(selection) && length(selection) > 0) {
          rna_filter_inputs[[cat]] <- selection
          rna_meta_subset <- rna_meta_subset %>% dplyr::filter(.data[[cat]] %in% selection)
        }
      }

      # Stage 2: Refine by selected combinations (if any filters were applied)
      if (length(rna_filter_inputs) > 0) {
        rna_cols_to_unite <- c("disease_category", "tissue")
        if (isTRUE(input$rna_comp_include_hub_status)) {
          rna_cols_to_unite <- c(rna_cols_to_unite, "hub_status")
        }
        rna_cols_to_unite <- intersect(rna_cols_to_unite, colnames(rna_meta_subset))

        if (length(rna_cols_to_unite) > 0) {
          rna_meta_with_comb <- rna_meta_subset %>%
            tidyr::unite("combination", all_of(rna_cols_to_unite), sep = " & ", remove = FALSE, na.rm = TRUE)

          # FIX: If refined selection is NULL (e.g. during update), default to all available combinations
          # instead of filtering to nothing.
          if (!is.null(input$rna_refined_selection)) {
            rna_cell_ids <- rna_meta_with_comb %>%
              dplyr::filter(combination %in% input$rna_refined_selection) %>%
              rownames()
          } else {
            rna_cell_ids <- rownames(rna_meta_subset)
          }
        } else {
          rna_cell_ids <- rownames(rna_meta_subset)
        }
      }

      # --- ATAC Data Filtering ---
      req(rv$dogma_object)
      atac_meta_subset <- rv$dogma_object@meta.data
      atac_cell_ids <- rownames(atac_meta_subset) # Default to all cells

      # Stage 1: Broad filtering
      atac_filter_inputs <- list()
      for (cat in ANNOTATION_CATEGORIES_ATAC) {
        sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
        selection <- input[[paste0("atac_comp_filter_", sanitized_cat)]]
        if (!is.null(selection) && length(selection) > 0) {
          atac_filter_inputs[[cat]] <- selection
          atac_meta_subset <- atac_meta_subset %>% dplyr::filter(.data[[cat]] %in% selection)
        }
      }

      # Stage 2: Refine by selected combinations
      if (length(atac_filter_inputs) > 0) {
        req(input$atac_refined_selection)

        atac_cols_to_unite <- c("Disease Category (ATAC-seq)", "Tissue (ATAC-seq)")
        if (isTRUE(input$atac_comp_include_hub_status)) {
          atac_cols_to_unite <- c(atac_cols_to_unite, "Predicted Hub Status (ATAC-seq)")
        }
        atac_cols_to_unite <- intersect(atac_cols_to_unite, colnames(atac_meta_subset))

        if (length(atac_cols_to_unite) > 0) {
          atac_meta_with_comb <- atac_meta_subset %>%
            tidyr::unite("combination", all_of(atac_cols_to_unite), sep = " & ", remove = FALSE, na.rm = TRUE)

          atac_cell_ids <- atac_meta_with_comb %>%
            dplyr::filter(combination %in% input$atac_refined_selection) %>%
            rownames()
        } else {
          atac_cell_ids <- rownames(atac_meta_subset)
        }
      }

      list(rna_cells = rna_cell_ids, atac_cells = atac_cell_ids)
    },
    ignoreNULL = FALSE
  )

  # --- Combination Selector Logic for Tab 3: RNA ---
  rna_comp_combinations <- reactive({
    meta_subset <- rv$seurat_object@meta.data
    any_filter_active <- FALSE
    for (cat in ANNOTATION_CATEGORIES_RNA) {
      sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
      selection <- input[[paste0("rna_comp_filter_", sanitized_cat)]]
      if (!is.null(selection) && length(selection) > 0) {
        any_filter_active <- TRUE
        meta_subset <- meta_subset %>% dplyr::filter(.data[[cat]] %in% selection)
      }
    }
    if (!any_filter_active) {
      return(character(0))
    }

    cols_to_unite <- c("disease_category", "tissue")
    if (isTRUE(input$rna_comp_include_hub_status)) {
      cols_to_unite <- c(cols_to_unite, "hub_status")
    }
    cols_to_unite <- intersect(cols_to_unite, colnames(meta_subset))
    if (length(cols_to_unite) == 0) {
      return(character(0))
    }

    meta_subset %>%
      dplyr::distinct(across(all_of(cols_to_unite))) %>%
      tidyr::unite("combination", all_of(cols_to_unite), sep = " & ", remove = TRUE, na.rm = TRUE) %>%
      pull(combination) %>%
      sort()
  })

  output$rna_comp_combination_selector_ui <- renderUI({
    combinations <- rna_comp_combinations()
    if (length(combinations) == 0) {
      return(NULL)
    }
    tagList(
      fluidRow(
        column(6, actionLink("select_all_rna_comp_combinations", "Select All", style = "font-size: 9pt;")),
        column(6, actionLink("deselect_all_rna_comp_combinations", "Deselect All", style = "font-size: 9pt;"))
      ),
      checkboxGroupInput(
        inputId = "rna_refined_selection",
        label = NULL,
        choices = combinations,
        selected = combinations,
        inline = TRUE
      )
    )
  })

  observeEvent(input$select_all_rna_comp_combinations, {
    updateCheckboxGroupInput(session, "rna_refined_selection", selected = rna_comp_combinations())
  })
  observeEvent(input$deselect_all_rna_comp_combinations, {
    updateCheckboxGroupInput(session, "rna_refined_selection", selected = character(0))
  })

  # --- Combination Selector Logic for Tab 3: ATAC ---
  if (HAS_ATAC_DATA) {
    atac_comp_combinations <- reactive({
      req(rv$dogma_object)
      meta_subset <- rv$dogma_object@meta.data
      any_filter_active <- FALSE
      for (cat in ANNOTATION_CATEGORIES_ATAC) {
        sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
        selection <- input[[paste0("atac_comp_filter_", sanitized_cat)]]
        if (!is.null(selection) && length(selection) > 0) {
          any_filter_active <- TRUE
          meta_subset <- meta_subset %>% dplyr::filter(.data[[cat]] %in% selection)
        }
      }
      if (!any_filter_active) {
        return(character(0))
      }

      cols_to_unite <- c("Disease Category (ATAC-seq)", "Tissue (ATAC-seq)")
      if (isTRUE(input$atac_comp_include_hub_status)) {
        cols_to_unite <- c(cols_to_unite, "Predicted Hub Status (ATAC-seq)")
      }
      cols_to_unite <- intersect(cols_to_unite, colnames(meta_subset))
      if (length(cols_to_unite) == 0) {
        return(character(0))
      }

      meta_subset %>%
        dplyr::distinct(across(all_of(cols_to_unite))) %>%
        tidyr::unite("combination", all_of(cols_to_unite), sep = " & ", remove = TRUE, na.rm = TRUE) %>%
        pull(combination) %>%
        sort()
    })

    output$atac_comp_combination_selector_ui <- renderUI({
      combinations <- atac_comp_combinations()
      if (length(combinations) == 0) {
        return(NULL)
      }
      tagList(
        fluidRow(
          column(6, actionLink("select_all_atac_comp_combinations", "Select All", style = "font-size: 9pt;")),
          column(6, actionLink("deselect_all_atac_comp_combinations", "Deselect All", style = "font-size: 9pt;"))
        ),
        checkboxGroupInput(
          inputId = "atac_refined_selection",
          label = NULL,
          choices = combinations,
          selected = combinations,
          inline = TRUE
        )
      )
    })

    observeEvent(input$select_all_atac_comp_combinations, {
      updateCheckboxGroupInput(session, "atac_refined_selection", selected = atac_comp_combinations())
    })
    observeEvent(input$deselect_all_atac_comp_combinations, {
      updateCheckboxGroupInput(session, "atac_refined_selection", selected = character(0))
    })
  }

  if (HAS_ATAC_DATA) {
    atac_plot_data <- reactive({
      gene <- selected_comparison_gene()
      req(gene)
      req(rv$dogma_object)
      if (!"umap" %in% names(rv$dogma_object@reductions)) {
        return(NULL)
      }

      plot_df <- as.data.frame(rv$dogma_object@reductions$umap@cell.embeddings) %>% rownames_to_column("cell_id")

      if ("GENEACCESS" %in% names(rv$dogma_object@assays)) {
        accessibility_data <- GetAssayData(rv$dogma_object, assay = "GENEACCESS", layer = "data")[gene, , drop = FALSE]
        acc_df <- data.frame(cell_id = colnames(accessibility_data), value = as.numeric(accessibility_data[1, ]))
        plot_df <- plot_df %>% left_join(acc_df, by = "cell_id")
        plot_df <- plot_df[is.finite(plot_df$value), ]
      } else {
        return(NULL)
      }
      return(plot_df)
    })
  }

  output$rna_plot_title <- renderText({
    req(input$rna_comp_display_mode) # FIX: Wait for input

    rna_ids <- filtered_comparison_cells()$rna_cells
    total_ids <- Cells(rv$seurat_object)

    num_shown <- if (input$rna_comp_display_mode == "highlight") {
      length(total_ids)
    } else {
      length(total_ids) # "all" mode plots all cells
    }

    paste0("NeuMap RNA UMAP (", scales::comma(num_shown), " cells shown)")
  })
  output$atac_plot_title <- renderText({
    if (!HAS_ATAC_DATA) {
      return("ATAC data not available.")
    }
    req(input$atac_comp_display_mode) # FIX: Wait for input

    atac_ids <- filtered_comparison_cells()$atac_cells
    total_ids <- Cells(rv$dogma_object)

    num_shown <- if (input$atac_comp_display_mode == "highlight") {
      length(total_ids) # Highlight mode shows all ATAC cells (on RNA bg)
    } else {
      length(total_ids) # "all" mode also shows all ATAC cells
    }
    paste0("Promoter Accessibility (", scales::comma(num_shown), " cells shown)")
  })

  rna_expression_umap_plot <- reactive({
    req(input$rna_comp_display_mode) # FIX: Wait for input
    gene <- selected_comparison_gene()
    req(gene)
    base_df <- plot_data_base()
    req(base_df)
    filtered_ids <- filtered_comparison_cells()$rna_cells

    plot_df_full <- base_df[, c("UMAP_1", "UMAP_2", "cell_id")]

    expression_data <- data.frame(
      cell_id = colnames(rv$seurat_object),
      value = as.numeric(GetAssayData(rv$seurat_object, assay = "RNA", layer = "data")[gene, ]),
      stringsAsFactors = FALSE
    )

    plot_df_full <- plot_df_full %>% left_join(expression_data, by = "cell_id")

    # --- New Display Mode Logic (matches Tab 1) ---
    p_base <- ggplot() +
      labs(title = gene) +
      theme_void(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      scale_color_viridis_c(name = "Expression")

    if (input$rna_comp_display_mode == "highlight") {
      p_base <- p_base +
        geom_point(data = base_df, aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = input$rna_comp_pt_size * 0.8)

      if (length(filtered_ids) == 0) {
        return(p_base + annotate("text", x = mean(base_df$UMAP_1), y = mean(base_df$UMAP_2), label = "No cells match current selections.", size = 5))
      }
      plot_df <- plot_df_full %>% dplyr::filter(cell_id %in% filtered_ids)
    } else { # "all" mode (plots all cells colored by value)
      plot_df <- plot_df_full
    }
    # --- End New Logic ---

    plot_df_capped <- apply_capping(plot_df, "value")

    p <- p_base +
      geom_point(
        data = plot_df_capped %>% arrange(capped_value),
        aes(x = UMAP_1, y = UMAP_2, color = capped_value),
        size = input$rna_comp_pt_size,
        alpha = input$rna_comp_alpha
      )

    p
  })
  output$rna_expression_umap <- renderPlot({
    rna_expression_umap_plot()
  })

  if (HAS_ATAC_DATA) {
    atac_accessibility_umap_plot <- reactive({
      req(input$atac_comp_display_mode) # FIX: Wait for input
      background_df <- plot_data_base()
      req(background_df)
      plot_df_full <- atac_plot_data()
      filtered_ids <- filtered_comparison_cells()$atac_cells

      if (is.null(plot_df_full)) {
        return(ggplot() +
          annotate("text", x = 1, y = 1, label = "GENEACCESS assay missing or data unavailable.") +
          theme_void())
      }

      # --- New Display Mode Logic ---
      # This plot *always* needs the RNA background for context
      p_base <- ggplot() +
        geom_point(data = background_df, aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = input$rna_comp_pt_size * 0.8) +
        labs(title = selected_comparison_gene()) +
        theme_void(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        scale_color_viridis_c(name = "Accessibility")

      if (input$atac_comp_display_mode == "highlight") {
        if (length(filtered_ids) == 0) {
          return(p_base + annotate("text", x = mean(background_df$UMAP_1), y = mean(background_df$UMAP_2), label = "No cells match current selections.", size = 5))
        }
        plot_df <- plot_df_full %>% dplyr::filter(cell_id %in% filtered_ids)
      } else { # "all" mode (plots all ATAC cells colored by value)
        plot_df <- plot_df_full
      }
      # --- End New Logic ---

      plot_df_capped <- apply_capping(plot_df, "value")

      p <- p_base +
        geom_point(
          data = plot_df_capped %>% arrange(capped_value),
          aes(x = UMAP_1, y = UMAP_2, color = capped_value),
          size = input$atac_comp_pt_size,
          alpha = input$atac_comp_alpha
        )

      p
    })
    output$atac_accessibility_umap <- renderPlot({
      atac_accessibility_umap_plot()
    })
  }

  correlation_plot_object <- reactive({
    gene <- selected_comparison_gene()
    req(gene)
    rna_ids <- filtered_comparison_cells()$rna_cells

    rna_subset_obj <- subset(rv$seurat_object, cells = rna_ids)

    rna_data_for_corr <- FetchData(rna_subset_obj, vars = c(gene, "hub_status")) %>%
      setNames(c("expression", "hub")) %>%
      filter(hub != "NA", !is.na(hub))

    rna_data_trimmed <- apply_capping(rna_data_for_corr, "expression")

    rna_summary <- rna_data_trimmed %>%
      group_by(hub) %>%
      summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop")

    if (nrow(rna_summary) == 0) {
      return(ggplot() +
        annotate("text", x = 1, y = 1, label = "Could not calculate mean RNA expression in filtered subset.") +
        theme_void())
    }

    if (HAS_ATAC_DATA) {
      atac_ids <- filtered_comparison_cells()$atac_cells
      req(rv$dogma_object)
      atac_subset_obj <- subset(rv$dogma_object, cells = atac_ids)

      atac_data_for_corr <- tryCatch(
        {
          accessibility_scores <- GetAssayData(atac_subset_obj, assay = "GENEACCESS", layer = "data")[gene, , drop = TRUE]

          # Robustly extract metadata column
          hub_col_name <- "Predicted Hub Status (ATAC-seq)"
          if (!hub_col_name %in% colnames(atac_subset_obj@meta.data)) {
            # Fallback or error if column missing
            stop("Hub column not found")
          }

          hub_values <- atac_subset_obj@meta.data[[hub_col_name]]

          data.frame(
            accessibility = as.numeric(accessibility_scores),
            predicted_hub_status = hub_values
          ) %>%
            filter(!is.na(predicted_hub_status), predicted_hub_status != "<NA>") %>%
            rename(hub = predicted_hub_status)
        },
        error = function(e) {
          NULL
        }
      )

      if (is.null(atac_data_for_corr) || nrow(atac_data_for_corr) == 0) {
        return(ggplot() +
          annotate("text", x = 1, y = 1, label = "Could not get ATAC data for correlation in filtered subset.") +
          theme_void())
      }

      atac_data_trimmed <- apply_capping(atac_data_for_corr, "accessibility")

      atac_summary <- atac_data_trimmed %>%
        group_by(hub) %>%
        summarise(mean_accessibility = mean(accessibility, na.rm = TRUE), .groups = "drop")

      if (nrow(atac_summary) == 0) {
        return(ggplot() +
          annotate("text", x = 1, y = 1, label = "Could not calculate mean promoter accessibility in filtered subset.") +
          theme_void())
      }

      merged_data <- inner_join(rna_summary, atac_summary, by = "hub")

      if (nrow(merged_data) < 3) {
        return(ggplot() +
          annotate("text", x = 1, y = 1, label = "Fewer than 3 common groups to correlate in filtered subset.") +
          theme_void())
      }

      cor_test_result <- cor.test(merged_data$mean_accessibility, merged_data$mean_expression, method = "spearman")
      rho <- round(cor_test_result$estimate, 2)
      p_val <- format(cor_test_result$p.value, scientific = TRUE, digits = 2)

      # Use global colors if available
      corr_colors <- if (all(merged_data$hub %in% names(HUB_STATUS_COLORS))) {
        HUB_STATUS_COLORS[merged_data$hub]
      } else {
        generate_palette(length(unique(merged_data$hub)))
      }

      ggplot(merged_data, aes(x = mean_accessibility, y = mean_expression)) +
        geom_point(aes(color = hub), size = 4) +
        scale_color_manual(values = corr_colors) +
        geom_point(aes(color = hub), size = 4, show.legend = FALSE) +
        geom_smooth(method = "lm", se = FALSE, color = "grey50", linetype = "dashed", fullrange = TRUE) +
        ggrepel::geom_text_repel(aes(label = hub), size = 4) +
        labs(
          title = paste("Correlation for", gene),
          subtitle = paste0("Spearman's ρ = ", rho, " (p = ", p_val, ")"),
          x = "Mean Promoter Accessibility (Aggregated by Hub)",
          y = "Mean Gene Expression (Aggregated by Hub)"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11)
        )
    } else {
      return(ggplot() +
        annotate("text", x = 1, y = 1, label = "ATAC data not available for correlation analysis.") +
        theme_void())
    }
  })
  output$correlation_plot <- renderPlot({
    correlation_plot_object()
  })

  if (HAS_ATAC_DATA) {
    create_atac_reference_umap <- function(ref_var) {
      req(ref_var)
      base_df <- plot_data_base()
      req(base_df)
      req(rv$dogma_object)

      if (!ref_var %in% colnames(rv$dogma_object@meta.data)) {
        return(ggplot() +
          annotate("text", x = 1, y = 1, label = paste("Variable '", ref_var, "' not in ATAC data.")) +
          theme_void())
      }
      plot_df <- rv$dogma_object@meta.data %>% rownames_to_column("cell_id")

      umap_coords_atac <- as.data.frame(rv$dogma_object@reductions$umap@cell.embeddings) %>%
        rownames_to_column("cell_id")
      plot_df <- plot_df %>%
        left_join(umap_coords_atac, by = "cell_id")

      plot_df[[ref_var]] <- droplevels(as.factor(plot_df[[ref_var]]))
      plot_df <- plot_df[sample(nrow(plot_df)), ]

      p <- ggplot() +
        geom_point(data = base_df, aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 1.0, alpha = 0.5) +
        geom_point(data = plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[ref_var]]), size = 1.0, alpha = 0.8) +
        theme_void(base_size = 14) +
        labs(title = ref_var) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right") +
        guides(color = guide_legend(override.aes = list(size = 4), title = NULL))

      # Apply global colors if the variable matches Hub Status categories
      if (all(levels(plot_df[[ref_var]]) %in% names(HUB_STATUS_COLORS))) {
        p <- p + scale_color_manual(values = HUB_STATUS_COLORS)
      }

      return(p)
    }

    atac_reference_umap_tab3_plot <- reactive({
      create_atac_reference_umap(input$atac_reference_variable_tab3)
    })
    atac_reference_umap_tab4_plot <- reactive({
      create_atac_reference_umap(input$atac_reference_variable_tab4)
    })
    output$atac_reference_umap_tab3 <- renderPlot({
      atac_reference_umap_tab3_plot()
    })
    output$atac_reference_umap_tab4 <- renderPlot({
      atac_reference_umap_tab4_plot()
    })
  }

  if (!is.null(MOTIF_ASSAY_NAME)) {
    motif_filtered_cells <- reactive({
      meta_subset <- rv$dogma_object@meta.data
      any_filter_active <- FALSE
      for (cat in MOTIF_TAB_FILTERS) {
        sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
        selection <- input[[paste0("motif_filter_", sanitized_cat)]]
        if (!is.null(selection) && length(selection) > 0) {
          any_filter_active <- TRUE
          meta_subset <- meta_subset %>% dplyr::filter(.data[[cat]] %in% selection)
        }
      }
      if (!any_filter_active) {
        return(character(0))
      }
      return(rownames(meta_subset))
    })

    output$motif_cell_count_summary_ui <- renderUI({
      shiny::tags$div(
        style = "font-size: 0.9em; text-align: center; margin-top: 10px; padding: 5px; border: 1px solid #ddd;",
        strong("ATAC Cells in Selection:"),
        em(paste(scales::comma(length(motif_filtered_cells())), "cells"))
      )
    })

    observeEvent(input$run_motif_analysis, {
      withProgress(message = "Analyzing TF Motif Activity...", value = 0.1, {
        rv$bridge_logic_message <- NULL

        filter_inputs <- list()
        for (cat in MOTIF_TAB_FILTERS) {
          sanitized_cat <- gsub("[^a-zA-Z0-9_]", "_", cat)
          filter_inputs[[cat]] <- input[[paste0("motif_filter_", sanitized_cat)]]
        }
        active_filters <- filter_inputs[sapply(filter_inputs, function(x) !is.null(x) && length(x) > 0)]
        if (length(active_filters) < 1) {
          showNotification("Please select at least one filter category to define the overall cell population.", type = "warning", duration = 8)
          return()
        }

        setProgress(value = 0.2, detail = "Filtering cells...")
        meta_for_analysis <- rv$dogma_object@meta.data %>%
          tibble::rownames_to_column("cell_id")

        for (cat_name in names(active_filters)) {
          meta_for_analysis <- meta_for_analysis %>%
            dplyr::filter(!is.na(.data[[cat_name]]) & .data[[cat_name]] %in% active_filters[[cat_name]])
        }

        meta_for_plotting <- meta_for_analysis

        setProgress(value = 0.3, detail = "Defining cell groups...")
        if (isTRUE(input$motif_combinatorial_filter)) {
          meta_for_analysis <- meta_for_analysis %>%
            tidyr::unite("analysis_group", all_of(names(active_filters)), sep = " & ", remove = FALSE, na.rm = TRUE)
          summary_text_prefix <- "Combinatorial analysis of "
        } else {
          meta_for_analysis$analysis_group <- NA_character_

          filter_priority <- intersect(MOTIF_TAB_FILTERS, names(active_filters))

          for (cat_name in filter_priority) {
            is_match <- meta_for_analysis[[cat_name]] %in% active_filters[[cat_name]]
            is_unassigned <- is.na(meta_for_analysis$analysis_group)
            meta_for_analysis$analysis_group[is_match & is_unassigned] <- as.character(meta_for_analysis[[cat_name]][is_match & is_unassigned])
          }

          meta_for_analysis <- meta_for_analysis %>% dplyr::filter(!is.na(analysis_group))
          summary_text_prefix <- "Individual analysis of "
        }

        cells_for_analysis_df <- meta_for_analysis %>%
          group_by(analysis_group) %>%
          filter(n() > 10, analysis_group != "Not Selected", analysis_group != "") %>%
          ungroup()

        dogma_subset <- subset(rv$dogma_object, cells = cells_for_analysis_df$cell_id)

        group_map <- setNames(cells_for_analysis_df$analysis_group, cells_for_analysis_df$cell_id)
        dogma_subset$analysis_group <- group_map[Cells(dogma_subset)]
        dogma_subset <- subset(dogma_subset, subset = !is.na(analysis_group))

        Idents(dogma_subset) <- "analysis_group"

        if (nrow(cells_for_analysis_df) == 0) {
          showNotification("Error: No cells match the current filter selections.", type = "error", duration = 8)
          return()
        }
        group_counts <- table(Idents(dogma_subset))
        if (length(group_counts) < 2) {
          showNotification("Error: Fewer than two valid groups to compare. Please check your filters or choose a different comparison variable.", type = "error", duration = 10)
          return()
        }
        if (any(group_counts < 3)) {
          bad_groups <- names(group_counts[group_counts < 3])
          showNotification(paste("Error: The group(s) '", paste(bad_groups, collapse = ", "), "' have too few cells (< 3) for analysis. Adjust filters."), type = "error", duration = 10)
          return()
        }

        if (isTRUE(input$motif_combinatorial_filter)) {
          plot_df <- dogma_subset@meta.data %>%
            tibble::rownames_to_column("cell_id") %>%
            dplyr::select(cell_id, group = analysis_group)
        } else {
          plot_df_list <- list()
          for (cat_name in names(active_filters)) {
            for (item in active_filters[[cat_name]]) {
              cells_in_group <- meta_for_plotting %>%
                dplyr::filter(.data[[cat_name]] == item) %>%
                pull(cell_id)
              if (length(cells_in_group) > 0) {
                plot_df_list[[paste0(cat_name, "_", item)]] <- tibble(cell_id = cells_in_group, group = item)
              }
            }
          }
          plot_df <- dplyr::bind_rows(plot_df_list)
        }
        motif_plot_data(plot_df)

        setProgress(value = 0.6, detail = "Running FindAllMarkers...")
        test_to_use <- input$motif_test_use
        max_cells <- if (isTRUE(input$motif_enable_downsampling)) input$motif_max_cells else NULL

        motif_markers_raw <- FindAllMarkers(object = dogma_subset, assay = MOTIF_ASSAY_NAME, test.use = test_to_use, max.cells.per.ident = max_cells, logfc.threshold = 0.01, min.pct = 0.01, only.pos = FALSE)

        setProgress(value = 0.9, detail = "Formatting results...")
        rv$motif_analysis_summary <- paste0(summary_text_prefix, " ", length(group_counts), " groups across ", length(active_filters), " filter categories (", scales::comma(ncol(dogma_subset)), " total cells).")
        motif_analysis_subset_object(dogma_subset)
        rv$motif_analysis_groups <- levels(Idents(dogma_subset))
        rv$motif_results_raw <- motif_markers_raw

        if (!is.null(motif_markers_raw) && nrow(motif_markers_raw) > 0) {
          processed_results <- motif_markers_raw %>%
            dplyr::left_join(TF_LOOKUP_TABLE, by = c("gene" = "MotifID")) %>%
            mutate(p_val_adj = format(p_val_adj, scientific = TRUE, digits = 3), avg_log2FC = round(avg_log2FC, 3)) %>%
            dplyr::select(
              Group = cluster,
              `TF Name` = TF_Name,
              `Motif ID` = gene,
              `Avg. Activity Change (log2FC)` = avg_log2FC,
              `Adj. p-value` = p_val_adj
            )

          rv$motif_results_processed <- processed_results

          found_tfs <- unique(na.omit(processed_results$`TF Name`))
          if (length(found_tfs) > 0) {
            updateSelectizeInput(session, "viz_tf_select", choices = found_tfs, selected = found_tfs[1])
            updateSelectizeInput(session, "target_analysis_tf_select", choices = found_tfs, selected = found_tfs[1])
          }
        } else {
          rv$motif_results_processed <- data.frame(Info = "No differential motif activities found.")
        }
      })
    })

    output$motif_analysis_summary_text <- renderText({
      rv$motif_analysis_summary
    })

    output$motif_results_table <- DT::renderDT({
      req(rv$motif_results_processed)
      DT::datatable(
        rv$motif_results_processed,
        options = list(pageLength = 10, filter = "top", scrollX = TRUE),
        rownames = FALSE,
        selection = "single"
      )
    })

    observeEvent(input$motif_results_table_rows_selected, {
      req(input$motif_results_table_rows_selected, rv$motif_results_processed)

      selected_row_index <- input$motif_results_table_rows_selected

      if (selected_row_index > nrow(rv$motif_results_processed)) {
        return()
      }

      selected_tf <- rv$motif_results_processed$`TF Name`[selected_row_index]

      req(!is.na(selected_tf))

      updateSelectizeInput(session, "viz_tf_select", selected = selected_tf)

      updateSelectizeInput(session, "target_analysis_tf_select", selected = selected_tf)

      shinyjs::click("load_suggested_targets")
    })

    selected_viz_tf <- reactive({
      req(input$viz_tf_select, nchar(input$viz_tf_select) > 0, !is.na(input$viz_tf_select))
      input$viz_tf_select
    })

    visualization_plot_data <- reactive({
      req(motif_plot_data(), selected_viz_tf())
      plot_groups_df <- motif_plot_data()

      tf_name <- selected_viz_tf()
      motif_id <- get_motif_id_from_name_SAFE(tf_name)
      req(motif_id, rv$dogma_object)

      all_activity_scores <- GetAssayData(rv$dogma_object, assay = MOTIF_ASSAY_NAME, layer = "data")[motif_id, , drop = TRUE]

      if (is.null(all_activity_scores)) {
        showNotification(paste("Motif ID", motif_id, "not found in object."), type = "error")
        return(NULL)
      }

      activity_df <- data.frame(
        cell_id = names(all_activity_scores),
        activity = as.numeric(all_activity_scores)
      )

      plot_df <- plot_groups_df %>%
        left_join(activity_df, by = "cell_id") %>%
        mutate(is_positive = activity > input$motif_threshold) %>%
        filter(!is.na(group), !is.na(activity))

      req(nrow(plot_df) > 0)
      return(plot_df)
    })

    motif_activity_data <- reactive({
      req(motif_analysis_subset_object())
      dogma_subset <- motif_analysis_subset_object()

      tf_name <- selected_viz_tf()
      req(tf_name)
      motif_id <- get_motif_id_from_name_SAFE(tf_name)

      # --- Robust Data Retrieval (Bypass Reductions Slot) ---
      # Pull UMAP coordinates directly from metadata, which is guaranteed
      # to be present by global.R and correctly subsetted.
      if (!all(c("refUMAP_1", "refUMAP_2") %in% colnames(dogma_subset@meta.data))) {
        showNotification("CRITICAL: refUMAP coordinates not found in subsetted metadata.", type = "error")
        return(NULL)
      }

      plot_df <- dogma_subset@meta.data %>%
        dplyr::select(UMAP_1 = refUMAP_1, UMAP_2 = refUMAP_2) %>%
        rownames_to_column("cell_id")
      # --- End Robust Data Retrieval ---

      tryCatch(
        {
          activity_data <- GetAssayData(dogma_subset, assay = MOTIF_ASSAY_NAME, layer = "data")[motif_id, , drop = FALSE]
          if (ncol(activity_data) == 0) {
            return(NULL)
          }
          activity_df <- as.data.frame(t(as.matrix(activity_data))) %>%
            rownames_to_column("cell_id")
          colnames(activity_df)[2] <- "value"

          plot_df <- plot_df %>% left_join(activity_df, by = "cell_id")
        },
        error = function(e) {
          showNotification(paste("Error fetching motif activity data for ID:", motif_id), type = "error")
          return(NULL)
        }
      )

      req(plot_df, "value" %in% names(plot_df))
      plot_df$value <- as.numeric(plot_df$value)
      plot_df <- plot_df[is.finite(plot_df$value), ]
      return(plot_df)
    })

    motif_activity_plot_object <- reactive({
      plot_df <- motif_activity_data()
      req(plot_df)
      tf_name <- selected_viz_tf()

      background_df <- plot_data_base()
      req(background_df)

      plot_df_capped <- apply_capping(plot_df, "value")

      p <- ggplot() +
        geom_point(data = background_df, aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 1, alpha = 0.5) +
        geom_point(data = plot_df_capped %>% arrange(capped_value), aes(x = UMAP_1, y = UMAP_2, color = capped_value), size = 1.5, alpha = 0.7) +
        labs(title = tf_name) +
        theme_void(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)) +
        scale_color_viridis_c(name = "Activity Score")

      p
    })
    output$motif_activity_plot <- renderPlot({
      motif_activity_plot_object()
    })

    motif_subgroup_location_umap_plot <- reactive({
      plot_data <- motif_plot_data()
      req(plot_data)

      umap_coords <- as.data.frame(rv$dogma_object@reductions$umap@cell.embeddings) %>%
        rownames_to_column("cell_id")
      plot_df <- plot_data %>%
        left_join(umap_coords, by = "cell_id") %>%
        filter(!is.na(group))

      background_df <- plot_data_base()
      req(background_df)

      p <- ggplot(background_df, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(color = "grey85", size = 1.0, alpha = 0.5) +
        theme_void(base_size = 14) +
        theme(legend.position = "none")

      if (nrow(plot_df) > 0) {
        group_levels <- unique(plot_df$group)
        palette <- generate_palette(length(group_levels))
        color_map <- setNames(palette, group_levels)

        centroids <- plot_df %>%
          group_by(group) %>%
          summarise(UMAP_1 = mean(UMAP_1, na.rm = TRUE), UMAP_2 = mean(UMAP_2, na.rm = TRUE), .groups = "drop")

        p <- p + geom_point(data = plot_df, aes(color = group), size = 1.5, alpha = 0.7) +
          ggrepel::geom_label_repel(data = centroids, aes(label = group, fill = group), color = "white", fontface = "bold") +
          scale_color_manual(values = color_map) +
          scale_fill_manual(values = color_map)
      }
      p
    })
    output$motif_subgroup_location_umap <- renderPlot({
      motif_subgroup_location_umap_plot()
    })

    motif_group_color_map <- reactive({
      plot_data <- motif_plot_data()
      req(plot_data)

      group_levels <- sort(unique(plot_data$group))
      palette <- generate_palette(length(group_levels))
      color_map <- setNames(palette, group_levels)
      return(color_map)
    })

    motif_summary_stats_data <- reactive({
      plot_data <- visualization_plot_data()
      req(plot_data)

      data_for_stats <- apply_capping(plot_data, "activity")

      if (nrow(data_for_stats) == 0) {
        return(NULL)
      }

      percent_positive <- data_for_stats %>%
        group_by(group) %>%
        summarise(percent_positive = (sum(is_positive) / n()) * 100, .groups = "drop")

      mean_activity <- data_for_stats %>%
        dplyr::filter(is_positive) %>%
        group_by(group) %>%
        summarise(mean_activity = mean(activity, na.rm = TRUE), .groups = "drop")

      return(list(percent_positive = percent_positive, mean_activity = mean_activity))
    })

    motif_percent_positive_barchart_plot <- reactive({
      plot_data_list <- motif_summary_stats_data()
      req(plot_data_list, !is.null(plot_data_list$percent_positive), nrow(plot_data_list$percent_positive) > 0)

      color_map <- motif_group_color_map()

      plot_df <- plot_data_list$percent_positive %>%
        mutate(group = factor(group, levels = names(color_map)))

      p <- ggplot(plot_df, aes(x = group, y = percent_positive, fill = group)) +
        geom_col(show.legend = FALSE)

      if (isTRUE(input$show_barchart_labels_tab4)) {
        p <- p + geom_text(aes(label = paste0(round(percent_positive, 1), "%")), vjust = -0.5)
      }

      p + scale_fill_manual(values = color_map) +
        scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1))) +
        labs(x = NULL, y = "Percent (%)") +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
    })
    output$motif_percent_positive_barchart <- renderPlot({
      motif_percent_positive_barchart_plot()
    })

    motif_activity_barchart_object <- reactive({
      plot_data_list <- motif_summary_stats_data()
      req(plot_data_list, !is.null(plot_data_list$mean_activity), nrow(plot_data_list$mean_activity) > 0)

      color_map <- motif_group_color_map()

      plot_df <- plot_data_list$mean_activity %>%
        mutate(group = factor(group, levels = names(color_map)))

      p <- ggplot(plot_df, aes(x = group, y = mean_activity, fill = group)) +
        geom_col(show.legend = FALSE)

      if (isTRUE(input$show_barchart_labels_tab4)) {
        p <- p + geom_text(aes(label = round(mean_activity, 3)), vjust = -0.5)
      }

      p + scale_fill_manual(values = color_map) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(x = NULL, y = "Mean Score") +
        theme_minimal(base_size = 14) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
    })
    output$motif_activity_barchart <- renderPlot({
      motif_activity_barchart_object()
    })

    motif_activity_distribution_plot_object <- reactive({
      plot_df <- visualization_plot_data()
      req(plot_df)

      color_map <- motif_group_color_map()

      plot_df_capped <- apply_capping(plot_df, "activity") %>%
        mutate(group = factor(group, levels = names(color_map)))

      ggplot(plot_df_capped, aes(x = group, y = activity, fill = group)) +
        geom_violin(trim = FALSE, alpha = 0.8, show.legend = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, show.legend = FALSE) +
        scale_fill_manual(values = color_map) +
        labs(title = NULL, x = "Group", y = "Activity Score") +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank()
        )
    })
    output$motif_activity_distribution_plot <- renderPlot({
      motif_activity_distribution_plot_object()
    })

    upstream_plot_data_base <- reactive({
      req(motif_plot_data())

      plot_data_groups <- motif_plot_data()
      req(nrow(plot_data_groups) > 0)
      plot_groups <- unique(plot_data_groups$group)

      rna_cells_for_plot <- list()

      tissue_label_map <- setNames(c("bm", "lung", "spleen", "blood"), c("Bone Marrow", "Lung", "Spleen", "Blood"))
      disease_label_map <- setNames(c("Lung Cancer", "Acute inflammation"), c("LLC", "LPS"))

      for (group_label in plot_groups) {
        atac_cells_in_group <- plot_data_groups %>%
          dplyr::filter(group == group_label) %>%
          pull(cell_id)
        group_context_df <- rv$dogma_object@meta.data[atac_cells_in_group, , drop = FALSE]

        context_hubs <- unique(na.omit(as.character(group_context_df[["Predicted Hub Status (ATAC-seq)"]])))
        context_diseases <- unique(na.omit(as.character(group_context_df[["Disease Category (ATAC-seq)"]])))
        context_tissues <- unique(na.omit(as.character(group_context_df[["Tissue (ATAC-seq)"]])))

        rna_meta_subset_dh <- rv$seurat_object@meta.data
        if (length(context_hubs) > 0) {
          rna_meta_subset_dh <- rna_meta_subset_dh %>% dplyr::filter(hub_status %in% context_hubs)
        }
        if (length(context_diseases) > 0) {
          rna_equivalent_diseases <- disease_label_map[context_diseases]
          rna_equivalent_diseases[is.na(rna_equivalent_diseases)] <- context_diseases[is.na(rna_equivalent_diseases)]
          rna_meta_subset_dh <- rna_meta_subset_dh %>% dplyr::filter(disease_category %in% rna_equivalent_diseases)
        }

        final_rna_meta_subset <- rna_meta_subset_dh
        if (length(context_tissues) > 0 && nrow(rna_meta_subset_dh) > 0) {
          rna_equivalent_tissues <- tissue_label_map[context_tissues]
          rna_equivalent_tissues <- na.omit(rna_equivalent_tissues)
          if (length(rna_equivalent_tissues) > 0) {
            tissue_filtered_subset <- rna_meta_subset_dh %>% dplyr::filter(tissue %in% rna_equivalent_tissues)
            if (nrow(tissue_filtered_subset) > 0) {
              final_rna_meta_subset <- tissue_filtered_subset
            }
          }
        }

        rna_cells <- rownames(final_rna_meta_subset)
        if (length(rna_cells) > 0) {
          rna_cells_for_plot[[group_label]] <- rna_cells
        }
      }

      if (length(unlist(rna_cells_for_plot)) == 0) {
        return(NULL)
      }

      all_rna_cells <- unique(unlist(rna_cells_for_plot))
      rna_subset <- subset(rv$seurat_object, cells = all_rna_cells)

      cell_group_df <- purrr::map_dfr(names(rna_cells_for_plot), function(grp) {
        tibble::tibble(cell_id = rna_cells_for_plot[[grp]], analysis_group = grp)
      })
      unique_cell_group_df <- cell_group_df %>% dplyr::distinct(cell_id, .keep_all = TRUE)

      cell_to_group_map <- setNames(unique_cell_group_df$analysis_group, unique_cell_group_df$cell_id)

      rna_subset$analysis_group <- cell_to_group_map[Cells(rna_subset)]
      rna_subset <- subset(rna_subset, cells = which(!is.na(rna_subset$analysis_group)))

      return(rna_subset)
    })

    tf_expression_violin_plot_object <- reactive({
      tf_name_upper <- selected_viz_tf()
      req(tf_name_upper)

      # Ensure consistent color mapping is available
      req(motif_group_color_map())

      rna_subset <- upstream_plot_data_base()
      req(rna_subset, "analysis_group" %in% colnames(rna_subset@meta.data))

      all_rna_genes <- rownames(rna_subset)
      matched_gene <- all_rna_genes[toupper(all_rna_genes) == toupper(tf_name_upper)]

      if (length(matched_gene) == 0) {
        return(ggplot() +
          annotate("text", x = 1, y = 1, label = paste("Gene '", tf_name_upper, "' not found in RNA assay (case-insensitive search)."), size = 5) +
          labs(title = paste("RNA Expression for", tf_name_upper)) +
          theme_void())
      }

      tf_gene_to_plot <- matched_gene[1]

      plot_df <- FetchData(rna_subset, vars = c(tf_gene_to_plot, "analysis_group"))
      colnames(plot_df) <- c("expression", "group")

      plot_df <- plot_df %>%
        filter(!is.na(group)) %>%
        mutate(group = factor(group))

      if (nrow(plot_df) == 0) {
        return(ggplot() +
          annotate("text", x = 1, y = 1, label = "No corresponding RNA cells to plot.", size = 5) +
          theme_void())
      }

      plot_df_capped <- apply_capping(plot_df, "expression")

      plt <- ggplot(plot_df_capped, aes(x = group, y = expression, fill = group)) +
        geom_violin(trim = FALSE, alpha = 0.8, show.legend = FALSE) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, show.legend = FALSE)

      # Use same color map as Motif Activity Plot for consistency
      color_map <- motif_group_color_map()
      if (!is.null(color_map)) {
        plt <- plt + scale_fill_manual(values = color_map)
      }

      plt + labs(
        title = paste("RNA Expression of", tf_gene_to_plot, "in Corresponding Cell Groups"),
        x = "Analysis Group (from ATAC context)",
        y = "RNA Expression Level"
      ) +
        theme_minimal(base_size = 14) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid = element_blank()
        )
    })

    output$tf_expression_violin_plot <- renderPlot({
      tf_expression_violin_plot_object()
    })

    observeEvent(input$load_suggested_targets, {
      req(input$target_analysis_tf_select, input$dorothea_confidence_level, !is.null(dorothea_regulon_mm))

      tf_for_dorothea <- stringr::str_to_title(input$target_analysis_tf_select)

      withProgress(message = "Fetching DoRothEA Targets...", value = 0.5, {
        suggested_targets <- dorothea_regulon_mm %>%
          dplyr::filter(tf == !!tf_for_dorothea, confidence %in% input$dorothea_confidence_level) %>%
          pull(target) %>%
          unique()
        if (length(suggested_targets) > 0) {
          rv$target_gene_list <- suggested_targets
          updateTextAreaInput(session, "target_gene_display", value = paste(suggested_targets, collapse = "\n"))
          showNotification(paste(length(suggested_targets), "target genes loaded for TF", input$target_analysis_tf_select), type = "message")
        } else {
          rv$target_gene_list <- NULL
          updateTextAreaInput(session, "target_gene_display", value = "")
          showNotification(paste("No target genes found for TF", input$target_analysis_tf_select, "with selected confidence levels."), type = "warning")
        }
      })
    })

    output$target_heatmap_sorter_ui <- renderUI({
      req(input$target_heatmap_annotations)
      cols_to_display <- unique(c("Group", input$target_heatmap_annotations))
      if (length(cols_to_display) > 1) {
        rank_list(
          text = "Set Annotation Sorting Priority:",
          labels = cols_to_display,
          input_id = "target_annotation_order"
        )
      }
    })

    observeEvent(input$analyze_target_expression, {
      withProgress(message = "Analyzing Target Gene Expression...", value = 0.1, {
        req(input$target_analysis_tf_select, motif_plot_data())

        plot_data_groups <- motif_plot_data()
        plot_groups <- unique(plot_data_groups$group)

        setProgress(0.2, detail = "Mapping ATAC groups to RNA cells...")

        bridge_results <- list()
        tissue_label_map <- setNames(c("bm", "lung", "spleen", "blood"), c("Bone Marrow", "Lung", "Spleen", "Blood"))
        disease_label_map <- setNames(c("Lung Cancer", "Acute inflammation"), c("LLC", "LPS"))

        for (group_label in plot_groups) {
          atac_cells_in_group <- plot_data_groups %>%
            dplyr::filter(group == group_label) %>%
            pull(cell_id)
          req(rv$dogma_object)
          group_context_df <- rv$dogma_object@meta.data[atac_cells_in_group, , drop = FALSE]

          context_hubs <- unique(na.omit(as.character(group_context_df[["Predicted Hub Status (ATAC-seq)"]])))
          context_diseases <- unique(na.omit(as.character(group_context_df[["Disease Category (ATAC-seq)"]])))
          context_tissues <- unique(na.omit(as.character(group_context_df[["Tissue (ATAC-seq)"]])))

          rna_meta_subset_dh <- rv$seurat_object@meta.data
          if (length(context_hubs) > 0) {
            rna_meta_subset_dh <- rna_meta_subset_dh %>% dplyr::filter(hub_status %in% context_hubs)
          }
          if (length(context_diseases) > 0) {
            rna_equivalent_diseases <- disease_label_map[context_diseases]
            rna_equivalent_diseases[is.na(rna_equivalent_diseases)] <- context_diseases[is.na(rna_equivalent_diseases)]
            rna_meta_subset_dh <- rna_meta_subset_dh %>% dplyr::filter(disease_category %in% rna_equivalent_diseases)
          }

          final_rna_meta_subset <- rna_meta_subset_dh
          fallback_triggered <- TRUE
          if (length(context_tissues) > 0 && nrow(rna_meta_subset_dh) > 0) {
            rna_equivalent_tissues <- tissue_label_map[context_tissues]
            rna_equivalent_tissues <- na.omit(rna_equivalent_tissues)
            if (length(rna_equivalent_tissues) > 0) {
              tissue_filtered_subset <- rna_meta_subset_dh %>% dplyr::filter(tissue %in% rna_equivalent_tissues)
              if (nrow(tissue_filtered_subset) > 0) {
                final_rna_meta_subset <- tissue_filtered_subset
                fallback_triggered <- FALSE
              }
            }
          }

          bridge_results[[group_label]] <- list(cells = rownames(final_rna_meta_subset), fallback = fallback_triggered)
        }

        rna_cells_for_heatmap <- lapply(bridge_results, `[[`, "cells")

        if (length(unlist(rna_cells_for_heatmap)) == 0) {
          showNotification("Could not map any analysis groups to RNA cells.", type = "warning", duration = 10)
          return()
        }

        all_rna_cells <- unique(unlist(rna_cells_for_heatmap))
        rna_subset <- subset(rv$seurat_object, cells = all_rna_cells)

        any_fallback <- any(sapply(bridge_results, `[[`, "fallback"))
        if (any_fallback) {
          tissue_counts <- rna_subset@meta.data %>%
            dplyr::count(tissue, sort = TRUE) %>%
            mutate(percent = round(n / sum(n) * 100, 1))

          tissue_summary_str <- paste(paste0(tissue_counts$tissue, " (", tissue_counts$percent, "%)"), collapse = ", ")

          rv$bridge_logic_message <- paste0(
            "Displaying RNA cells that match the ATAC context by Disease and Hub, ",
            "sourced from the following tissues: ", tissue_summary_str
          )
        } else {
          rv$bridge_logic_message <- "<b>SUCCESS:</b> Displaying RNA cells matching the ATAC context (Disease, Tissue, and Hub)."
        }

        cell_group_df <- purrr::map_dfr(
          names(rna_cells_for_heatmap),
          function(group_name) {
            tibble::tibble(
              cell_id = rna_cells_for_heatmap[[group_name]],
              Group = group_name
            )
          }
        )

        unique_cell_group_df <- cell_group_df %>% dplyr::distinct(cell_id, .keep_all = TRUE)
        cell_to_group_map <- setNames(unique_cell_group_df$Group, unique_cell_group_df$cell_id)

        rna_subset$analysis_group <- cell_to_group_map[Cells(rna_subset)]
        rna_subset <- subset(rna_subset, cells = which(!is.na(rna_subset$analysis_group)))

        target_rna_subset(rna_subset) # This triggers the refinement UI

        # Clear previous DEG results to prevent showing a stale heatmap
        target_deg_results(NULL)

        setProgress(value = 1, detail = "RNA population identified. Refine and plot.")
        showNotification("Corresponding RNA population is ready. Refine selections and view heatmap.", type = "message")
      })
    })

    output$bridge_logic_summary_ui <- renderUI({
      req(rv$bridge_logic_message)

      is_warning <- stringr::str_detect(rv$bridge_logic_message, "WARNING")

      div(
        class = if (is_warning) "alert alert-warning" else "alert alert-success",
        style = "margin-bottom: 15px; padding: 10px; font-size: 0.9em;",
        HTML(rv$bridge_logic_message)
      )
    })

    refined_target_rna_subset <- reactive({
      req(target_rna_subset(), input$target_refined_selection)

      full_subset <- target_rna_subset()
      meta_subset <- full_subset@meta.data

      cols_to_unite <- c("disease_category", "tissue", "hub_status")
      if (!isTRUE(input$target_heatmap_include_hub_status)) {
        cols_to_unite <- setdiff(cols_to_unite, "hub_status")
      }
      cols_to_unite <- intersect(cols_to_unite, colnames(meta_subset))

      meta_with_comb <- meta_subset %>%
        tidyr::unite("combination", all_of(cols_to_unite), sep = " & ", remove = FALSE, na.rm = TRUE)

      cells_to_keep <- meta_with_comb %>%
        filter(combination %in% input$target_refined_selection) %>%
        rownames()

      if (length(cells_to_keep) == 0) {
        return(NULL)
      } # Prevent subset crash

      subset(full_subset, cells = cells_to_keep)
    })

    observe({
      # Gracefully handle NULL subset without crash
      subset_obj <- refined_target_rna_subset()
      if (is.null(subset_obj)) {
        target_deg_results(NULL)
        return()
      }
      req(target_deg_results, rv$target_gene_list)

      withProgress(message = "Running DE analysis on targets...", value = 0.5, {
        rna_subset <- subset_obj

        target_genes <- rv$target_gene_list
        genes_in_rna <- intersect(target_genes, rownames(rna_subset))

        if (length(genes_in_rna) == 0) {
          showNotification("None of the target genes were found in the refined RNA dataset.", type = "error")
          target_deg_results(data.frame())
          return()
        }

        if (nlevels(as.factor(rna_subset$analysis_group)) < 2) {
          showNotification("Fewer than two groups in refined selection. Cannot run DE analysis.", type = "warning")
          target_deg_results(data.frame())
          return()
        }

        Idents(rna_subset) <- "analysis_group"

        degs <- FindAllMarkers(
          rna_subset,
          features = genes_in_rna, logfc.threshold = 0,
          min.pct = 0.01, test.use = "wilcox", only.pos = FALSE
        )

        target_deg_results(degs)
      })
    })

    target_heatmap_args <- reactive({
      # Handle case where subset is NULL (empty selection)
      subset_obj <- refined_target_rna_subset()
      if (is.null(subset_obj)) {
        return("No cells selected.")
      }

      req(target_deg_results())

      deg_results <- target_deg_results()
      rna_subset <- subset_obj

      if (nrow(deg_results) == 0) {
        return("No differential genes found in refined selection.")
      }

      top_genes_df <- deg_results %>%
        dplyr::filter(abs(avg_log2FC) >= input$log2fc_target_heatmap_filter) %>%
        group_by(cluster) %>%
        slice_max(order_by = abs(avg_log2FC), n = input$n_target_heatmap_genes, with_ties = FALSE) %>%
        ungroup()

      unique_genes_to_plot <- unique(top_genes_df$gene)
      if (length(unique_genes_to_plot) == 0) {
        return("No genes pass the current filters.")
      }

      cells_to_plot <- tryCatch(
        {
          WhichCells(rna_subset, idents = levels(as.factor(rna_subset$analysis_group)), downsample = input$max_cells_target_heatmap)
        },
        error = function(e) {
          Cells(rna_subset)
        }
      )

      if (length(cells_to_plot) == 0) {
        return("No cells in the refined selection.")
      }

      primary_annot <- "analysis_group"
      extra_annots <- input$target_heatmap_annotations
      sorting_priority <- input$target_annotation_order

      cols_to_show <- unique(c(primary_annot, extra_annots))
      cols_to_show <- intersect(cols_to_show, colnames(rna_subset@meta.data))

      full_annotation_df <- rna_subset@meta.data[cells_to_plot, , drop = FALSE]

      final_sorting_order <- if (!is.null(sorting_priority) && all(sorting_priority %in% cols_to_show)) {
        sorting_priority
      } else {
        cols_to_show
      }

      ordered_annotation_df_with_id <- full_annotation_df %>%
        tibble::rownames_to_column("cell_id") %>%
        arrange(across(all_of(final_sorting_order)))

      ordered_cells <- ordered_annotation_df_with_id$cell_id
      if (length(ordered_cells) == 0) {
        return("No cells left after sorting.")
      }

      annotation_for_pheatmap <- ordered_annotation_df_with_id %>%
        dplyr::select(all_of(rev(final_sorting_order))) %>% # Reverse for pheatmap visual top-to-bottom match
        as.data.frame()
      rownames(annotation_for_pheatmap) <- ordered_cells

      data_for_heatmap <- GetAssayData(rna_subset, assay = "RNA", layer = "data")[unique_genes_to_plot, ordered_cells, drop = FALSE]
      scaled_data <- t(scale(t(as.matrix(data_for_heatmap))))
      cap_value <- quantile(abs(scaled_data), 0.98, na.rm = TRUE)
      if (!is.finite(cap_value) || cap_value < 1.5) cap_value <- 1.5
      scaled_data[is.na(scaled_data)] <- 0
      scaled_data[scaled_data > cap_value] <- cap_value
      scaled_data[scaled_data < -cap_value] <- -cap_value

      return(list(
        mat = scaled_data, show_colnames = FALSE, cluster_cols = FALSE,
        cluster_rows = (nrow(scaled_data) > 1), annotation_col = annotation_for_pheatmap,
        main = paste("Expression of Top", input$target_analysis_tf_select, "Target Genes")
      ))
    })

    output$target_gene_heatmap <- renderPlot({
      args <- target_heatmap_args()
      if (is.character(args)) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = stringr::str_wrap(args, width = 30), size = 6) +
          theme_void() +
          xlim(0, 1) +
          ylim(0, 1))
      }
      if (is.null(args) || is.null(args$mat) || ncol(args$mat) == 0) {
        return(ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "Select TF and targets, then click 'Analyze Expression'.", size = 6) +
          theme_void() +
          xlim(0, 1) +
          ylim(0, 1))
      }

      # Build custom color list for pheatmap
      annot_colors <- list()
      if (!is.null(args$annotation_col) && is.data.frame(args$annotation_col)) {
        for (col_name in names(args$annotation_col)) {
          levels <- levels(factor(args$annotation_col[[col_name]]))
          default_palette <- setNames(generate_palette(length(levels)), levels)
          custom_palette <- rv$color_maps$heatmap_annotations[[col_name]]
          if (!is.null(custom_palette)) {
            final_palette <- unlist(utils::modifyList(as.list(default_palette), custom_palette))
          } else {
            final_palette <- default_palette
          }
          annot_colors[[col_name]] <- final_palette
        }
      }

      final_args <- c(args, list(annotation_colors = annot_colors))

      do.call(pheatmap::pheatmap, final_args)
    })

    target_rna_combinations <- reactive({
      req(target_rna_subset())
      meta_subset <- target_rna_subset()@meta.data

      cols_to_unite <- c("disease_category", "tissue", "hub_status")
      if (!isTRUE(input$target_heatmap_include_hub_status)) {
        cols_to_unite <- setdiff(cols_to_unite, "hub_status")
      }
      cols_to_unite <- intersect(cols_to_unite, colnames(meta_subset))

      if (length(cols_to_unite) < 1) {
        return(character(0))
      }

      meta_subset %>%
        dplyr::distinct(across(all_of(cols_to_unite))) %>%
        tidyr::unite("combination", all_of(cols_to_unite), sep = " & ", remove = TRUE, na.rm = TRUE) %>%
        pull(combination) %>%
        sort()
    })

    output$target_combination_selector_ui <- renderUI({
      combinations <- target_rna_combinations()
      req(length(combinations) > 0)

      tagList(
        fluidRow(
          column(6, actionLink("select_all_target_combinations", "Select All", style = "font-size: 9pt;")),
          column(6, actionLink("deselect_all_target_combinations", "Deselect All", style = "font-size: 9pt;"))
        ),
        checkboxGroupInput(
          inputId = "target_refined_selection",
          label = NULL,
          choices = combinations,
          selected = combinations
        )
      )
    })

    observeEvent(input$select_all_target_combinations, {
      updateCheckboxGroupInput(session, "target_refined_selection", selected = target_rna_combinations())
    })
    observeEvent(input$deselect_all_target_combinations, {
      updateCheckboxGroupInput(session, "target_refined_selection", selected = character(0))
    })

    output$download_target_heatmap_csv <- downloadHandler(
      filename = function() {
        paste0("Target_Gene_Heatmap_Data_", input$target_analysis_tf_select, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        args <- target_heatmap_args()

        # Handle error messages by writing them to the CSV
        if (is.character(args)) {
          write.csv(data.frame(Error = args), file, row.names = FALSE)
          return()
        }

        if (!is.list(args) || is.null(args$mat)) {
          write.csv(data.frame(Error = "No data available for heatmap."), file, row.names = FALSE)
          return()
        }
        final_data <- as.data.frame(as.matrix(args$mat)) %>% rownames_to_column("gene")
        write.csv(final_data, file, row.names = FALSE)
      }
    )

    observeEvent(input$run_upstream_analysis, {
      req(input$upstream_gene_input, dorothea_regulon_mm)

      target_genes <- str_trim(unlist(str_split(input$upstream_gene_input, "\\s+")))
      target_genes <- target_genes[target_genes != ""]

      if (length(target_genes) == 0) {
        showNotification("Please enter at least one gene symbol.", type = "warning")
        rv$upstream_tf_results <- NULL
        return()
      }

      results <- dorothea_regulon_mm %>%
        dplyr::filter(target %in% target_genes) %>%
        dplyr::select(Regulator = tf, Target = target, Confidence = confidence) %>%
        distinct()

      if (nrow(results) == 0) {
        showNotification("No known regulators found for the entered gene(s) in the DoRothEA database.", type = "warning")
        rv$upstream_tf_results <- NULL
        return()
      }

      showNotification(paste("Found", n_distinct(results$Regulator), "potential regulators."), type = "message")
      rv$upstream_tf_results <- results
    })

    output$upstream_tf_table <- DT::renderDT({
      req(rv$upstream_tf_results)
      DT::datatable(rv$upstream_tf_results,
        options = list(pageLength = 5),
        selection = "multiple",
        rownames = FALSE
      )
    })

    visualize_selected_tfs <- eventReactive(input$visualize_selected_regulators, {
      req(rv$upstream_tf_results, input$upstream_tf_table_rows_selected)
      selected_indices <- input$upstream_tf_table_rows_selected
      rv$upstream_tf_results %>%
        slice(selected_indices) %>%
        pull(Regulator) %>%
        unique()
    })

    upstream_plot_data_base <- reactive({
      req(motif_plot_data())

      plot_data_groups <- motif_plot_data()
      req(nrow(plot_data_groups) > 0)
      plot_groups <- unique(plot_data_groups$group)

      rna_cells_for_plot <- list()

      tissue_label_map <- setNames(c("bm", "lung", "spleen", "blood"), c("Bone Marrow", "Lung", "Spleen", "Blood"))
      disease_label_map <- setNames(c("Lung Cancer", "Acute inflammation"), c("LLC", "LPS"))

      for (group_label in plot_groups) {
        atac_cells_in_group <- plot_data_groups %>%
          dplyr::filter(group == group_label) %>%
          pull(cell_id)
        req(rv$dogma_object)
        group_context_df <- rv$dogma_object@meta.data[atac_cells_in_group, , drop = FALSE]

        context_hubs <- unique(na.omit(as.character(group_context_df[["Predicted Hub Status (ATAC-seq)"]])))
        context_diseases <- unique(na.omit(as.character(group_context_df[["Disease Category (ATAC-seq)"]])))
        context_tissues <- unique(na.omit(as.character(group_context_df[["Tissue (ATAC-seq)"]])))

        rna_meta_subset_dh <- rv$seurat_object@meta.data
        if (length(context_hubs) > 0) {
          rna_meta_subset_dh <- rna_meta_subset_dh %>% dplyr::filter(hub_status %in% context_hubs)
        }
        if (length(context_diseases) > 0) {
          rna_equivalent_diseases <- disease_label_map[context_diseases]
          rna_equivalent_diseases[is.na(rna_equivalent_diseases)] <- context_diseases[is.na(rna_equivalent_diseases)]
          rna_meta_subset_dh <- rna_meta_subset_dh %>% dplyr::filter(disease_category %in% rna_equivalent_diseases)
        }

        final_rna_meta_subset <- rna_meta_subset_dh
        if (length(context_tissues) > 0 && nrow(rna_meta_subset_dh) > 0) {
          rna_equivalent_tissues <- tissue_label_map[context_tissues]
          rna_equivalent_tissues <- na.omit(rna_equivalent_tissues)
          if (length(rna_equivalent_tissues) > 0) {
            tissue_filtered_subset <- rna_meta_subset_dh %>% dplyr::filter(tissue %in% rna_equivalent_tissues)
            if (nrow(tissue_filtered_subset) > 0) {
              final_rna_meta_subset <- tissue_filtered_subset
            }
          }
        }

        rna_cells <- rownames(final_rna_meta_subset)
        if (length(rna_cells) > 0) {
          rna_cells_for_plot[[group_label]] <- rna_cells
        }
      }

      if (length(unlist(rna_cells_for_plot)) == 0) {
        return(NULL)
      }

      all_rna_cells <- unique(unlist(rna_cells_for_plot))
      rna_subset <- subset(rv$seurat_object, cells = all_rna_cells)

      cell_group_df <- purrr::map_dfr(names(rna_cells_for_plot), function(grp) {
        tibble::tibble(cell_id = rna_cells_for_plot[[grp]], analysis_group = grp)
      })
      unique_cell_group_df <- cell_group_df %>% dplyr::distinct(cell_id, .keep_all = TRUE)

      cell_to_group_map <- setNames(unique_cell_group_df$analysis_group, unique_cell_group_df$cell_id)

      rna_subset$analysis_group <- cell_to_group_map[Cells(rna_subset)]
      rna_subset <- subset(rna_subset, cells = which(!is.na(rna_subset$analysis_group)))

      return(rna_subset)
    })

    upstream_composite_plot_object <- reactive({
      tfs_to_plot <- visualize_selected_tfs()
      req(length(tfs_to_plot) > 0)

      rna_subset <- upstream_plot_data_base()
      req(rna_subset, "analysis_group" %in% colnames(rna_subset@meta.data))

      plot_list <- lapply(tfs_to_plot, function(tf_name) {
        if (!tf_name %in% rownames(rna_subset)) {
          return(ggplot() +
            annotate("text", x = 1, y = 1, label = paste("TF '", tf_name, "' not found in RNA data.")) +
            theme_void())
        }

        plot_df <- FetchData(rna_subset, vars = c(tf_name, "analysis_group"))
        colnames(plot_df) <- c("expression", "group")

        plot_df <- plot_df %>% filter(!is.na(group))
        if (nrow(plot_df) == 0) {
          return(ggplot() +
            annotate("text", x = 1, y = 1, label = paste("No cells to plot for", tf_name)) +
            theme_void())
        }

        plot_df_capped <- apply_capping(plot_df, "expression")

        plot_df_capped <- plot_df_capped %>% mutate(is_positive = expression > input$upstream_threshold)

        percent_positive_data <- plot_df_capped %>%
          group_by(group) %>%
          summarise(percent_positive = (sum(is_positive, na.rm = TRUE) / n()) * 100, .groups = "drop")

        p_percent <- ggplot(percent_positive_data, aes(x = group, y = percent_positive, fill = group)) +
          geom_col(show.legend = FALSE) +
          geom_text(aes(label = paste0(round(percent_positive, 1), "%")), vjust = -0.5) +
          scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1))) +
          labs(x = NULL, y = "Percent (%)", title = "% Positive Cells") +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_blank()
          )

        mean_expr_data <- plot_df_capped %>%
          filter(is_positive == TRUE) %>%
          group_by(group) %>%
          summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop")

        p_mean <- ggplot(mean_expr_data, aes(x = group, y = mean_expr, fill = group)) +
          geom_col(show.legend = FALSE) +
          geom_text(aes(label = round(mean_expr, 2)), vjust = -0.5) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
          labs(x = NULL, y = "Mean Expression", title = "Mean Expression (Positive Cells)") +
          theme_minimal(base_size = 12) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_blank()
          )

        p_violin <- ggplot(plot_df_capped, aes(x = group, y = expression, fill = group)) +
          geom_violin(trim = FALSE, alpha = 0.8, show.legend = FALSE) +
          geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, show.legend = FALSE) +
          labs(
            title = paste("Expression Distribution for", tf_name)
          ) +
          theme_minimal(base_size = 14) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.grid = element_blank()
          )

        (p_percent + p_mean) / p_violin
      })

      patchwork::wrap_plots(plot_list, ncol = 1)
    })

    output$upstream_plots_ui <- renderUI({
      tfs <- visualize_selected_tfs()
      req(tfs)
      plot_height <- max(400, length(tfs) * 500)
      plotOutput("upstream_composite_plot", height = paste0(plot_height, "px"))
    })

    output$upstream_composite_plot <- renderPlot({
      upstream_composite_plot_object()
    })

    output$download_upstream_plot <- create_download_handler(
      "upstream_composite_plot",
      upstream_composite_plot_object,
      "Upstream_Regulator_Expression",
      feature_name_reactive = visualize_selected_tfs
    )

    output$download_rna_comparison_umap <- create_download_handler("rna_expression_umap", rna_expression_umap_plot, "RNA_Comparison_UMAP", feature_name_reactive = selected_comparison_gene)
    output$download_atac_comparison_umap <- create_download_handler("atac_accessibility_umap", atac_accessibility_umap_plot, "ATAC_Comparison_UMAP", feature_name_reactive = selected_comparison_gene)
    output$download_correlation_plot <- create_download_handler("correlation_plot", correlation_plot_object, "RNA_ATAC_Correlation", feature_name_reactive = selected_comparison_gene)
    output$download_atac_ref_umap_tab3 <- create_download_handler("atac_reference_umap_tab3", atac_reference_umap_tab3_plot, "ATAC_Reference_UMAP", feature_name_reactive = reactive(input$atac_reference_variable_tab3))

    output$download_motif_activity_umap <- create_download_handler("motif_activity_plot", motif_activity_plot_object, "Motif_Activity_UMAP", feature_name_reactive = selected_viz_tf)
    output$download_motif_subgroup_umap <- create_download_handler("motif_subgroup_location_umap", motif_subgroup_location_umap_plot, "Motif_Subgroup_UMAP", feature_name_reactive = reactive("Subgroups"))
    output$download_atac_ref_umap_tab4 <- create_download_handler("atac_reference_umap_tab4", atac_reference_umap_tab4_plot, "ATAC_Reference_UMAP", feature_name_reactive = reactive(input$atac_reference_variable_tab4))

    output$download_motif_percent_barchart <- create_download_handler("motif_percent_positive_barchart", motif_percent_positive_barchart_plot, "Motif_Percent_Positive", feature_name_reactive = selected_viz_tf)
    output$download_motif_barchart <- create_download_handler("motif_activity_barchart", motif_activity_barchart_object, "Motif_Mean_Activity", feature_name_reactive = selected_viz_tf)
    output$download_motif_dist_plot <- create_download_handler("motif_activity_distribution_plot", motif_activity_distribution_plot_object, "Motif_Distribution", feature_name_reactive = selected_viz_tf)

    output$download_target_heatmap <- downloadHandler(
      filename = function() {
        paste0("Target_Gene_Heatmap_", input$target_analysis_tf_select, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        args <- target_heatmap_args()

        if (is.character(args)) {
          pdf(file, width = 8, height = 6)
          grid::grid.newpage()
          grid::grid.text(stringr::str_wrap(args, width = 30), x = 0.5, y = 0.5)
          dev.off()
          return()
        }

        if (!is.list(args) || is.null(args$mat)) {
          pdf(file, width = 8, height = 6)
          grid::grid.newpage()
          grid::grid.text("No data available for heatmap.", x = 0.5, y = 0.5)
          dev.off()
          return()
        }

        plot_size <- input$target_gene_heatmap_size
        width_in <- if (!is.null(plot_size$width)) as.numeric(plot_size$width) / 96 else 8
        height_in <- if (!is.null(plot_size$height)) as.numeric(plot_size$height) / 96 else 10

        # Build custom color list for pheatmap (Consistent with screen output)
        annot_colors <- list()
        if (!is.null(args$annotation_col) && is.data.frame(args$annotation_col)) {
          for (col_name in names(args$annotation_col)) {
            levels <- levels(factor(args$annotation_col[[col_name]]))
            default_palette <- setNames(generate_palette(length(levels)), levels)
            custom_palette <- rv$color_maps$heatmap_annotations[[col_name]]
            if (!is.null(custom_palette)) {
              final_palette <- unlist(utils::modifyList(as.list(default_palette), custom_palette))
            } else {
              final_palette <- default_palette
            }
            annot_colors[[col_name]] <- final_palette
          }
        }

        pheatmap::pheatmap(args$mat,
          filename = file, width = width_in, height = height_in,
          show_colnames = args$show_colnames, cluster_cols = args$cluster_cols,
          cluster_rows = args$cluster_rows, annotation_col = args$annotation_col,
          annotation_colors = annot_colors,
          main = args$main
        )
      }
    )

    output$download_tf_expression_violin <- create_download_handler(
      "tf_expression_violin_plot",
      tf_expression_violin_plot_object,
      "TF_RNA_Expression_Violin",
      feature_name_reactive = selected_viz_tf
    )
  }

  output$debug_motif_text <- renderPrint({
    cat("--- Debugging for TF Motif Activity Analysis ---\n\n")
    cat("1. Live Filter State:\n")
    if (!is.null(MOTIF_ASSAY_NAME)) {
      cat(paste(" Combinatorial Mode Enabled:", isTRUE(input$motif_combinatorial_filter), "\n"))
      cat(paste(" ATAC Cells in Selection:", scales::comma(length(motif_filtered_cells())), "\n"))
    } else {
      cat(" Analysis disabled.\n")
    }

    cat("\n2. Post-Analysis Object Status:\n")
    subset_obj <- motif_analysis_subset_object()
    if (is.null(subset_obj)) {
      cat(" motif_analysis_subset_object(): NULL\n")
    } else {
      cat(paste0(" motif_analysis_subset_object(): Seurat object with ", ncol(subset_obj), " cells and ", nrow(subset_obj), " features.\n"))
      cat(" Final Group Counts in Subset:\n")
      print(table(Idents(subset_obj)))
    }

    cat("\n3. Post-Analysis Results Status:\n")
    if (is.null(rv$motif_results_raw)) {
      cat(" motif_results_raw: NULL\n")
    } else {
      cat(paste0(" motif_results_raw: data.frame with ", nrow(rv$motif_results_raw), " rows.\n"))
      print(head(rv$motif_results_raw))
    }
  })

  output$debug_tf_expression_text <- renderPrint({
    cat("--- Debugging for TF Gene Expression Violin Plot ---\n\n")

    tf_name_upper <- tryCatch(selected_viz_tf(), error = function(e) "Not Selected")
    cat(paste("1. TF selected in UI:", tf_name_upper, "\n"))

    rna_subset <- tryCatch(upstream_plot_data_base(), error = function(e) NULL)
    if (is.null(rna_subset)) {
      cat("2. Corresponding RNA subset object (upstream_plot_data_base) is NULL.\n")
      return()
    }
    cat(paste("2. Corresponding RNA subset object is valid with", ncol(rna_subset), "cells.\n"))

    all_rna_genes <- rownames(rna_subset)
    cat(paste("3. Total genes in RNA subset:", length(all_rna_genes), "\n"))

    cat("   Example RNA gene names:", paste(head(all_rna_genes), collapse = ", "), "\n")

    matched_gene <- all_rna_genes[toupper(all_rna_genes) == toupper(tf_name_upper)]

    cat(paste("4. Case-insensitive search result for '", tf_name_upper, "':\n", sep = ""))
    if (length(matched_gene) == 0) {
      cat("   -> NO MATCH FOUND\n")
    } else {
      cat(paste("   -> MATCH FOUND:", paste(matched_gene, collapse = ", "), "\n"))
    }
  })

  # --- Custom Color Menu Logic ---

  # Dynamic UI for Tab 1 Color Pickers
  output$tab1_color_pickers_ui <- renderUI({
    req(input$tab1_color_target)

    if (input$tab1_color_target == "Population Groups") {
      plot_data <- selections_for_tab1()
      req(plot_data)

      groups <- plot_data %>%
        filter(plot_group != "Not Selected" & plot_group != "All Cells") %>%
        pull(plot_group) %>%
        unique() %>%
        sort()

      if (length(groups) == 0) {
        return(p("No population groups selected.", style = "color: #666; font-size: 0.9em;"))
      }

      default_palette <- generate_palette(length(groups))
      default_color_map <- setNames(default_palette, groups)

      tagList(lapply(groups, function(group_name) {
        sanitized_id <- gsub("[^a-zA-Z0-9_]", "_", group_name)
        input_id <- paste0("tab1_color_pop_", sanitized_id)
        current_color <- rv$color_maps$population[[group_name]] %||% default_color_map[[group_name]]

        shinyWidgets::colorPickr(
          inputId = input_id,
          label = group_name,
          selected = current_color
        )
      }))
    } else if (input$tab1_color_target == "UMAP Density") {
      palette_choices <- list(
        "Viridis" = c("viridis", "magma", "plasma", "inferno", "cividis"),
        "ColorBrewer" = c("RdYlBu", "Spectral", "Blues", "Reds", "Greens", "YlOrRd")
      )

      selectInput(
        inputId = "tab1_density_palette",
        label = "2. Select Palette:",
        choices = palette_choices,
        selected = rv$color_maps$density_palette
      )
    } else {
      p("This coloring feature is not available for this tab.", style = "color: #666; font-size: 0.9em;")
    }
  })

  # Observer to update central color state from Tab 1 inputs
  observe({
    plot_data <- selections_for_tab1()
    if (is.null(plot_data)) {
      return()
    }

    groups <- plot_data %>%
      filter(plot_group != "Not Selected" & plot_group != "All Cells") %>%
      pull(plot_group) %>%
      unique()

    for (group_name in groups) {
      sanitized_id <- gsub("[^a-zA-Z0-9_]", "_", group_name)
      input_id <- paste0("tab1_color_pop_", sanitized_id)

      if (!is.null(input[[input_id]])) {
        if (!isTRUE(all.equal(rv$color_maps$population[[group_name]], input[[input_id]]))) {
          rv$color_maps$population[[group_name]] <- input[[input_id]]
        }
      }
    }
  })

  observeEvent(input$tab1_density_palette, {
    req(input$tab1_density_palette)
    rv$color_maps$density_palette <- input$tab1_density_palette
  })

  # --- New Logic for Heatmap Color Menus ---

  # Dynamic UI for Tab 2 (DEG Heatmap)
  output$tab2_color_pickers_ui <- renderUI({
    req(input$tab2_color_target == "Heatmap Annotations")
    args <- tryCatch(deg_heatmap_args(), error = function(e) NULL)
    req(is.list(args), !is.null(args$annotation_col))
    annot_df <- args$annotation_col

    if (ncol(annot_df) == 0) {
      return(p("No annotations to color.", style = "color:#666;"))
    }

    lapply(names(annot_df), function(col_name) {
      levels <- levels(factor(annot_df[[col_name]]))
      default_palette <- setNames(generate_palette(length(levels)), levels)
      tagList(
        h5(str_to_title(gsub("_", " ", col_name))),
        lapply(levels, function(level_name) {
          sanitized_col <- gsub("[^a-zA-Z0-9_]", "_", col_name)
          sanitized_level <- gsub("[^a-zA-Z0-9_]", "_", level_name)
          input_id <- paste0("tab2_color_annot_", sanitized_col, "_", sanitized_level)
          current_color <- rv$color_maps$heatmap_annotations[[col_name]][[level_name]] %||% default_palette[[level_name]]
          shinyWidgets::colorPickr(inputId = input_id, label = level_name, selected = current_color)
        }),
        hr(style = "border-top:1px dotted #ccc;")
      )
    })
  })

  # Dynamic UI for Tab 4 (Target Gene Heatmap)
  output$tab4_color_pickers_ui <- renderUI({
    req(input$tab4_color_target == "Heatmap Annotations")
    args <- tryCatch(target_heatmap_args(), error = function(e) NULL)
    req(is.list(args), !is.null(args$annotation_col))
    annot_df <- args$annotation_col

    if (ncol(annot_df) == 0) {
      return(p("No annotations to color.", style = "color:#666;"))
    }

    lapply(names(annot_df), function(col_name) {
      levels <- levels(factor(annot_df[[col_name]]))
      default_palette <- setNames(generate_palette(length(levels)), levels)
      tagList(
        h5(str_to_title(gsub("_", " ", col_name))),
        lapply(levels, function(level_name) {
          sanitized_col <- gsub("[^a-zA-Z0-9_]", "_", col_name)
          sanitized_level <- gsub("[^a-zA-Z0-9_]", "_", level_name)
          input_id <- paste0("tab4_color_annot_", sanitized_col, "_", sanitized_level)
          current_color <- rv$color_maps$heatmap_annotations[[col_name]][[level_name]] %||% default_palette[[level_name]]
          shinyWidgets::colorPickr(inputId = input_id, label = level_name, selected = current_color)
        }),
        hr(style = "border-top:1px dotted #ccc;")
      )
    })
  })

  # Observer for Tab 2 Heatmap Colors
  observe({
    args <- tryCatch(deg_heatmap_args(), error = function(e) NULL)
    req(is.list(args), !is.null(args$annotation_col))
    annot_df <- args$annotation_col

    for (col_name in names(annot_df)) {
      for (level_name in levels(factor(annot_df[[col_name]]))) {
        sanitized_col <- gsub("[^a-zA-Z0-9_]", "_", col_name)
        sanitized_level <- gsub("[^a-zA-Z0-9_]", "_", level_name)
        input_id <- paste0("tab2_color_annot_", sanitized_col, "_", sanitized_level)

        input_val <- input[[input_id]]
        if (!is.null(input_val)) {
          if (is.null(rv$color_maps$heatmap_annotations[[col_name]])) rv$color_maps$heatmap_annotations[[col_name]] <- list()
          if (!isTRUE(all.equal(rv$color_maps$heatmap_annotations[[col_name]][[level_name]], input_val))) {
            rv$color_maps$heatmap_annotations[[col_name]][[level_name]] <- input_val
          }
        }
      }
    }
  })

  # Observer for Tab 4 Heatmap Colors
  observe({
    args <- tryCatch(target_heatmap_args(), error = function(e) NULL)
    req(is.list(args), !is.null(args$annotation_col))
    annot_df <- args$annotation_col

    for (col_name in names(annot_df)) {
      for (level_name in levels(factor(annot_df[[col_name]]))) {
        sanitized_col <- gsub("[^a-zA-Z0-9_]", "_", col_name)
        sanitized_level <- gsub("[^a-zA-Z0-9_]", "_", level_name)
        input_id <- paste0("tab4_color_annot_", sanitized_col, "_", sanitized_level)

        input_val <- input[[input_id]]
        if (!is.null(input_val)) {
          if (is.null(rv$color_maps$heatmap_annotations[[col_name]])) rv$color_maps$heatmap_annotations[[col_name]] <- list()
          if (!isTRUE(all.equal(rv$color_maps$heatmap_annotations[[col_name]][[level_name]], input_val))) {
            rv$color_maps$heatmap_annotations[[col_name]][[level_name]] <- input_val
          }
        }
      }
    }
  })
}
