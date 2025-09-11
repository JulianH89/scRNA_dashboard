# app.R

# Set the maximum file upload size to 500MB
options(shiny.maxRequestSize = 500*1024^2)

# Load necessary libraries
library(shiny)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(DT)
library(plotly)

# Source module and utility files
source("R/utils.R")
source("R/01_upload_module.R")
source("R/02_dge_module.R")
# source("R/03_pathway_module.R")
# source("R/04_umap_module.R")

# Define the UI
ui <- navbarPage(
  title = "scRNA-Seq Analysis Dashboard",
  theme = shinytheme("cosmo"),
  
  # Sector 1: Data Upload and Metadata
  tabPanel(
    "1. Data Upload & Metadata",
    upload_ui("upload_module") 
  ),
  
  # Sector 2: Differential Gene Expression
  tabPanel(
    "2. Differential Expression",
    dge_ui("dge_module") 
  ),
  
  # Sector 3: Pathway Enrichment (Placeholder)
  tabPanel(
    "3. Pathway Enrichment",
    h3("Gene Pathway Enrichment Analysis (Coming Soon)")
    # pathway_ui("pathway_module")
  ),
  
  # Sector 4: UMAP Visualization (Placeholder)
  tabPanel(
    "4. UMAP Visualization",
    h3("Interactive UMAP Embedding (Coming Soon)")
    # umap_ui("umap_module")
  )
)

# Define the server logic
server <- function(input, output, session) {
  
  # Call the server module for the upload section.
  # This module returns the reactive Seurat object.
  seurat_data_reactive <- upload_server("upload_module")
  
  # Call the DGE server module, passing it the reactive Seurat object.
  dge_server("dge_module", seurat_data_reactive)
  
  # pathway_server("pathway_module", ...) 
  # umap_server("umap_module", seurat_data_reactive)
  
}

# Run the application
shinyApp(ui = ui, server = server)