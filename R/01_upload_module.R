# R/01_upload_module.R
# Contains the UI and server logic for the data upload and metadata display tab.

# ======================================================================================
# --------------------------------- UI FUNCTION ----------------------------------------
# ======================================================================================
#' @param id A string, the namespace id for the module.
upload_ui <- function(id) {
  # Create a namespace function using the provided id
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        h4("Data Input"),
        # File input for user to upload their own Seurat object
        fileInput(ns("upload_seurat"), "Upload Seurat Object (.rds)",
                  accept = c(".rds")),
        helpText("Upload a Seurat object with raw counts. Normalization and UMAP will be generated if missing."),
        hr(),
        # Button to load an example dataset
        actionButton(ns("load_example"), "Load Example PBMC Data"),
        p("Note: The example dataset must be in a 'data' subfolder.", style = "font-size: smaller;"),
        width = 3
      ),
      mainPanel(
        h3("Cell Metadata"),
        p("This table displays the metadata for each cell in the dataset. Use the filters at the top of each column to search and subset the data."),
        # Output for the interactive metadata table
        DT::dataTableOutput(ns("metadata_table"))
      )
    )
  )
}

# ======================================================================================
# -------------------------------- SERVER FUNCTION -------------------------------------
# ======================================================================================
#' @param id A string, the namespace id for the module.
#' @return A reactive expression containing the loaded Seurat object.
upload_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # Create a reactive value to store the loaded Seurat object.
    seurat_object <- reactiveVal(NULL)

    # --- Main Data Loading and Processing Logic ---
    handle_data_load <- function(data_path) {
      progress <- shiny::Progress$new()
      progress$set(message = "Loading & Processing Data", value = 0)
      on.exit(progress$close())
      
      tryCatch({
        # Step 1: Read the object
        progress$inc(0.2, detail = "Reading file...")
        loaded_obj <- readRDS(data_path)
        
        # Step 2: Validate the object (calls function from R/utils.R)
        progress$inc(0.2, detail = "Validating structure...")
        validation_result <- validate_seurat_object(loaded_obj)
        
        if (validation_result$is_valid) {
          # Step 3: Process if necessary (calls function from R/utils.R)
          processed_obj <- process_seurat_object(loaded_obj, progress)
          seurat_object(processed_obj) # Set the reactive value
          progress$set(value = 1, detail = "Ready!")
          showNotification("Seurat object loaded and processed successfully.", type = "message")
        } else {
          seurat_object(NULL) # Reset on failure
          showModal(modalDialog(
            title = "Invalid Seurat Object",
            p(validation_result$message),
            easyClose = TRUE, footer = NULL
          ))
        }
      }, error = function(e) {
        seurat_object(NULL)
        showModal(modalDialog(
          title = "Error Loading Data",
          p("An error occurred while reading or processing the file:"),
          p(e$message),
          easyClose = TRUE, footer = NULL
        ))
      })
    }

    # --- 1. Handle user file upload ---
    observeEvent(input$upload_seurat, {
      req(input$upload_seurat)
      handle_data_load(input$upload_seurat$datapath)
    })

    # --- 2. Handle example data loading ---
    observeEvent(input$load_example, {
      example_file <- "data/pbmc_tutorial.rds"
      
      if (file.exists(example_file)) {
        handle_data_load(example_file)
      } else {
        showModal(modalDialog(
          title = "Error: Example Data Not Found",
          paste0("The example file was not found at '", example_file, "'. ",
                 "Please create a 'data' folder and place a pre-processed Seurat object inside."),
          easyClose = TRUE, footer = NULL
        ))
      }
    })

    # --- 3. Render the interactive metadata table ---
    output$metadata_table <- DT::renderDataTable({
      req(seurat_object())
      DT::datatable(
        seurat_object()@meta.data,
        options = list(scrollX = TRUE, pageLength = 15, searching = TRUE),
        rownames = TRUE, filter = 'top'
      )
    })
    
    # --- 4. Return the reactive Seurat object ---
    return(seurat_object)
  })
}
