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
        # Dropdown to select the metadata column for grouping
        selectInput(ns("grouping_var"), "Group Cells By:", choices = NULL),
        
        # Dropdown to select the first group for comparison
        selectInput(ns("group1"), "Group 1 (Test Group):", choices = NULL),
        
        # Dropdown to select the second group for comparison
        selectInput(ns("group2"), "Group 2 (Control Group):", choices = NULL),
        helpText("Select 'All Other Cells' to compare Group 1 against every other cell."),
        
        # Action button to trigger the DGE analysis
        actionButton(ns("run_dge"), "Run Differential Expression"),
        width = 3
      ),
      mainPanel(
        h3("Differential Expression Results"),
        p("This table shows the top markers differentiating the two selected groups. ",
          "Positive values in 'avg_log2FC' indicate the gene is more highly expressed in Group 1."),
        # Output for the interactive results table
        DT::dataTableOutput(ns("dge_table"))
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
    
    # Update the grouping variable choices when a new Seurat object is loaded
    observe({
      req(seurat_data())
      # We only want metadata columns that define discrete groups
      metadata <- seurat_data()@meta.data
      potential_groups <- sapply(metadata, function(x) !is.numeric(x) && length(unique(x)) > 1)
      updateSelectInput(session, "grouping_var", choices = names(metadata)[potential_groups])
    })
    
    # Update the choices for Group 1 and Group 2 based on the selected grouping variable
    observeEvent(input$grouping_var, {
      req(seurat_data(), input$grouping_var)
      groups <- unique(seurat_data()@meta.data[[input$grouping_var]])
      sorted_groups <- sort(as.character(groups))
      
      updateSelectInput(session, "group1", choices = sorted_groups)
      # Add the special option for comparing against all other cells
      updateSelectInput(session, "group2", choices = c("All Other Cells", sorted_groups))
    })
    
    # --- 2. DGE Analysis ---
    
    # Create a reactive value to store the results data frame
    dge_results <- reactiveVal(NULL)
    
    # Run the analysis when the user clicks the action button
    observeEvent(input$run_dge, {
      req(seurat_data(), input$grouping_var, input$group1, input$group2)
      
      # Ensure Group 1 and Group 2 are not the same
      if (input$group1 == input$group2) {
        showModal(modalDialog(
          title = "Invalid Selection",
          "Group 1 and Group 2 cannot be the same. Please select different groups to compare.",
          easyClose = TRUE, footer = NULL
        ))
        return() # Stop execution
      }
      
      progress <- shiny::Progress$new()
      progress$set(message = "Running DGE Analysis", value = 0.3)
      on.exit(progress$close())
      
      # Call the helper function from R/utils.R to perform the analysis
      results <- perform_dge(
        obj = seurat_data(),
        grouping_var = input$grouping_var,
        group1 = input$group1,
        group2 = input$group2,
        progress = progress
      )
      
      dge_results(results) # Store the results
      progress$set(value = 1, detail = "Done!")
    })
    
    # --- 3. Render the output table ---
    output$dge_table <- DT::renderDataTable({
      req(dge_results())
      DT::datatable(
        dge_results(),
        options = list(scrollX = TRUE, pageLength = 15),
        rownames = FALSE
      )
    })
    
  })
}