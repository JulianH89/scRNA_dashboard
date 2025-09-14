# R/02_dge_module.R
# Contains the UI and server logic for the Differential Gene Expression tab.

# ======================================================================================
# --------------------------------- UI FUNCTION ----------------------------------------
# ======================================================================================
#' @param id A string, the namespace id for the module.
dge_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        h4("DGE Parameters"),
        selectInput(ns("grouping_var"), "Group Cells By:", choices = NULL),
        selectInput(ns("group1"), "Group 1 (Test Group):", choices = NULL),
        selectInput(ns("group2"), "Group 2 (Control Group):", choices = NULL),
        helpText("Select 'All Other Cells' to compare Group 1 against every other cell."),
        actionButton(ns("run_dge"), "Run Differential Expression"),
        
        # Add a separator and controls for plot thresholds
        hr(),
        h4("Plotting Thresholds"),
        numericInput(ns("pval_threshold"), "Adjusted P-value Threshold:", 
                     value = 0.05, min = 0, max = 1, step = 0.01),
        sliderInput(ns("fc_threshold"), "Log2 Fold Change Threshold:", 
                    min = 0, max = 5, value = 1, step = 0.1),
        
        width = 3
      ),
      mainPanel(
        fluidRow(
          column(width = 12,
            h3("Volcano Plot"),
            p("This plot visualizes the relationship between fold change and statistical significance. Hover over points to see gene details."),
            plotlyOutput(ns("volcano_plot"), height = "500px")
          )
        ),
        hr(),
        fluidRow(
          column(width = 12,
            h3("Differential Expression Results"),
            p("This table shows the top markers differentiating the two selected groups. ",
              "Positive values in 'avg_log2FC' indicate the gene is more highly expressed in Group 1."),
            DT::dataTableOutput(ns("dge_table"))
          )
        )
      )
    )
  )
}

# ======================================================================================
# -------------------------------- SERVER FUNCTION -------------------------------------
# ======================================================================================
#' @param id A string, the namespace id for the module.
#' @param seurat_data A reactive expression containing the loaded Seurat object.
dge_server <- function(id, seurat_data) {
  moduleServer(id, function(input, output, session) {
    
    # --- 1. Dynamic UI Updates ---
    # ... (code is unchanged) ...
    observe({
      req(seurat_data())
      metadata <- seurat_data()@meta.data
      potential_groups <- sapply(metadata, function(x) !is.numeric(x) && length(unique(x)) > 1)
      updateSelectInput(session, "grouping_var", choices = names(metadata)[potential_groups])
    })
    observeEvent(input$grouping_var, {
      req(seurat_data(), input$grouping_var)
      groups <- unique(seurat_data()@meta.data[[input$grouping_var]])
      sorted_groups <- sort(as.character(groups))
      updateSelectInput(session, "group1", choices = sorted_groups)
      updateSelectInput(session, "group2", choices = c("All Other Cells", sorted_groups))
    })

    # --- 2. DGE Analysis ---
    # ... (code is unchanged) ...
    dge_results <- reactiveVal(NULL)
    observeEvent(input$run_dge, {
      req(seurat_data(), input$grouping_var, input$group1, input$group2)
      if (input$group1 == input$group2) {
        showModal(modalDialog(
          title = "Invalid Selection",
          "Group 1 and Group 2 cannot be the same. Please select different groups to compare.",
          easyClose = TRUE, footer = NULL
        ))
        return()
      }
      progress <- shiny::Progress$new()
      progress$set(message = "Running DGE Analysis", value = 0.3)
      on.exit(progress$close())
      results <- perform_dge(
        obj = seurat_data(),
        grouping_var = input$grouping_var,
        group1 = input$group1,
        group2 = input$group2,
        progress = progress
      )
      dge_results(results)
      progress$set(value = 1, detail = "Done!")
    })
    
    # --- 3. Render Outputs ---
    
    # Render the volcano plot, passing the reactive threshold inputs
    output$volcano_plot <- renderPlotly({
      req(dge_results())
      generate_volcano_plot(
        dge_results(),
        fc_threshold = input$fc_threshold,
        pval_threshold = input$pval_threshold
      )
    })
    
    # Render the output table
    output$dge_table <- DT::renderDataTable({
      req(dge_results())
      DT::datatable(
        dge_results(),
        options = list(scrollX = TRUE, pageLength = 10),
        rownames = FALSE
      )
    })
    
    # --- 4. Return reactive values for other modules to use ---
    return(
      list(
        # The DGE results data frame
        results = dge_results,
        # The thresholds from the UI inputs
        thresholds = reactive({
          list(
            pval = input$pval_threshold,
            fc = input$fc_threshold
          )
        })
      )
    )

  })
}