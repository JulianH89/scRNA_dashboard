# R/04_umap_module.R
# Contains the UI and server logic for the UMAP Visualization tab.

# ======================================================================================
# --------------------------------- UI FUNCTION ----------------------------------------
# ======================================================================================
#' @param id A string, the namespace id for the module.
umap_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        h4("UMAP Plot Controls"),
        # Dropdown to select metadata for coloring
        selectInput(ns("color_by"), "Color Cells By:", choices = NULL),
        
        # Dropdown to select specific groups to highlight
        selectizeInput(ns("highlight_clusters"), "Highlight Groups:", choices = NULL, multiple = TRUE),
        
        hr(),
        h4("Top Marker Genes"),
        p("Select a single group in the 'Highlight Groups' dropdown to display its top 10 marker genes below."),
        width = 3
      ),
      mainPanel(
        h3("Interactive UMAP Embedding"),
        plotlyOutput(ns("umap_plot"), height = "600px"),
        hr(),
        h3("Top 10 Marker Genes for Selected Group"),
        DT::dataTableOutput(ns("marker_table"))
      )
    )
  )
}

# ======================================================================================
# -------------------------------- SERVER FUNCTION -------------------------------------
# ======================================================================================
#' @param id A string, the namespace id for the module.
#' @param seurat_data A reactive expression containing the loaded Seurat object.
umap_server <- function(id, seurat_data) {
  moduleServer(id, function(input, output, session) {
    
    # --- 1. Dynamic UI Updates ---
    # Update the 'color_by' dropdown based on the metadata of the loaded object
    observe({
      req(seurat_data())
      # Get metadata columns that are factors or have a reasonable number of unique values
      metadata <- seurat_data()@meta.data
      potential_cols <- sapply(metadata, function(col) is.factor(col) || is.character(col) || length(unique(col)) < 50)
      updateSelectInput(session, "color_by", choices = names(which(potential_cols)), selected = "seurat_clusters")
    })
    
    # Update the 'highlight_clusters' dropdown based on the selected 'color_by' column
    observeEvent(input$color_by, {
      req(seurat_data(), input$color_by)
      choices <- unique(seurat_data()@meta.data[[input$color_by]])
      updateSelectizeInput(session, "highlight_clusters", choices = sort(choices))
    })
    
    # --- 2. UMAP Plot Rendering ---
    output$umap_plot <- renderPlotly({
      req(seurat_data(), input$color_by)
      
      seurat_obj <- seurat_data()
      
      # FIX: Dynamically find the UMAP reduction name instead of hardcoding 'umap'
      reduc_name <- NULL
      if ("umap" %in% Reductions(seurat_obj)) {
        reduc_name <- "umap"
      } else {
        # Find any reduction that contains "umap" in its name as a fallback
        reduc_name <- grep("umap", Reductions(seurat_obj), value = TRUE)[1]
      }
      
      # Stop if no UMAP reduction is found at all
      if (is.na(reduc_name) || is.null(reduc_name)) {
        showModal(modalDialog(
          title = "UMAP Not Found",
          "Could not find a UMAP reduction in the Seurat object.",
          easyClose = TRUE
        ))
        req(FALSE) # Stop execution
      }
      
      # Extract data for plotting using the dynamically found reduction name
      umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction = reduc_name))
      metadata <- seurat_obj@meta.data
      plot_data <- cbind(umap_coords, metadata)
      
      # FIX: Get the actual column names (e.g., "UMAP_1" or "umap_1")
      umap_cols <- colnames(umap_coords)
      req(length(umap_cols) >= 2) # Ensure we have at least 2 dimensions
      
      # Logic for highlighting selected clusters
      if (!is.null(input$highlight_clusters) && length(input$highlight_clusters) > 0) {
        plot_data$highlight_alpha <- ifelse(plot_data[[input$color_by]] %in% input$highlight_clusters, 1, 0.1)
        plot_data$highlight_color <- ifelse(plot_data[[input$color_by]] %in% input$highlight_clusters, as.character(plot_data[[input$color_by]]), "unselected")
      } else {
        plot_data$highlight_alpha <- 1
        plot_data$highlight_color <- as.character(plot_data[[input$color_by]])
      }

      # FIX: Create the plot using the dynamic column names
      plot_ly(
        data = plot_data,
        x = ~get(umap_cols[1]), # Use get() to refer to the column by its string name
        y = ~get(umap_cols[2]),
        color = ~highlight_color,
        type = 'scatter', mode = 'markers',
        marker = list(size = 3, opacity = ~highlight_alpha),
        text = ~paste("Cell:", rownames(plot_data), "<br>Group:", plot_data[[input$color_by]]),
        hoverinfo = 'text'
      ) %>% layout(
          title = "UMAP Visualization", 
          legend = list(title = list(text = input$color_by)),
          xaxis = list(title = umap_cols[1]), # Set axis titles dynamically
          yaxis = list(title = umap_cols[2])
      )
    })

    # --- 3. Marker Gene Table ---
    marker_results <- reactive({
      req(seurat_data(), input$highlight_clusters)
      # Only run if exactly one cluster is selected
      if (length(input$highlight_clusters) == 1) {
        find_top_markers(
          seurat_obj = seurat_data(),
          group_by_var = input$color_by,
          ident_1 = input$highlight_clusters
        )
      } else {
        NULL
      }
    })
    
    output$marker_table <- DT::renderDataTable({
      req(marker_results())
      DT::datatable(
        marker_results(),
        options = list(pageLength = 10),
        rownames = FALSE
      )
    })
    
  })
}
