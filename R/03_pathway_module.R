# R/03_pathway_module.R
# Contains the UI and server logic for the Pathway Enrichment tab.

# ======================================================================================
# --------------------------------- UI FUNCTION ----------------------------------------
# ======================================================================================
#' @param id A string, the namespace id for the module.
pathway_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    sidebarLayout(
      sidebarPanel(
        h4("Enrichment Parameters"),
        p("This analysis uses the significant genes from the DGE results."),
        
        # Select which genes to analyze
        selectInput(ns("gene_set"), "Select Gene Set:",
                    choices = c("Upregulated" = "up", 
                                "Downregulated" = "down", 
                                "All Significant" = "all")),
        
        # Select Gene Ontology
        selectInput(ns("ontology"), "Select Ontology:",
                    choices = c("Biological Process" = "BP", 
                                "Molecular Function" = "MF", 
                                "Cellular Component" = "CC")),
        
        actionButton(ns("run_enrichment"), "Run Enrichment Analysis"),
        
        hr(),
        p("Note: Analysis is performed on genes meeting the significance thresholds set on the DGE page."),
        width = 3
      ),
      mainPanel(
        h3("Gene Ontology Enrichment Results"),
        p("The plot below shows the top enriched pathways. The size of the dot corresponds to the number of genes found in that pathway, and the color indicates the statistical significance."),
        plotOutput(ns("enrichment_plot"), height = "600px"),
        hr(),
        h3("Enrichment Results Table"),
        DT::dataTableOutput(ns("enrichment_table"))
      )
    )
  )
}

# ======================================================================================
# -------------------------------- SERVER FUNCTION -------------------------------------
# ======================================================================================
#' @param id A string, the namespace id for the module.
#' @param dge_results A reactive expression containing the DGE results data frame.
#' @param dge_thresholds A reactive list containing p_val and fc thresholds from the DGE module.
pathway_server <- function(id, dge_results, dge_thresholds) {
  moduleServer(id, function(input, output, session) {
    
    enrichment_results <- eventReactive(input$run_enrichment, {
      req(dge_results())
      
      # First, call the reactive expression to get the list of thresholds
      thresholds <- dge_thresholds() 
      # Then, access the elements from the list
      fc_thresh <- thresholds$fc
      pval_thresh <- thresholds$pval

      # Filter the DGE results to get significant genes
      significant_genes <- dge_results() %>%
        filter(p_val_adj < pval_thresh, abs(avg_log2FC) > fc_thresh)

      # Select gene set based on user input
      gene_list <- switch(input$gene_set,
        "up"   = filter(significant_genes, avg_log2FC > 0) %>% pull(gene),
        "down" = filter(significant_genes, avg_log2FC < 0) %>% pull(gene),
        "all"  = pull(significant_genes, gene)
      )

      if (length(gene_list) < 5) {
        showModal(modalDialog(
          title = "Not Enough Genes",
          "Fewer than 5 significant genes were found with the current thresholds. Please adjust the DGE thresholds and re-run.",
          easyClose = TRUE, footer = nil
        ))
        return(NULL)
      }

      progress <- shiny::Progress$new()
      progress$set(message = "Running GO Analysis", value = 0)
      on.exit(progress$close())

      # Call the helper function
      perform_pathway_analysis(
        gene_list = gene_list,
        ontology = input$ontology,
        progress = progress
      )
    })

    output$enrichment_plot <- renderPlot({
      req(enrichment_results())
      # The result from the helper is an enrichResult object, which enrichplot can use
      dotplot(enrichment_results(), showCategory = 20)
    })
    
    output$enrichment_table <- DT::renderDataTable({
      req(enrichment_results())
      # Convert the enrichResult object to a data frame for display
      results_df <- as.data.frame(enrichment_results())
      DT::datatable(
        results_df,
        options = list(scrollX = TRUE, pageLength = 10),
        rownames = FALSE
      )
    })
    
  })
}