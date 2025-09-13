# R/utils.R
# This script contains helper functions for validating, processing, and analyzing the Seurat object.


#' Validate Seurat Object
#'
#' Checks if the uploaded object is a valid Seurat object with a raw count matrix
#' in the 'RNA' assay. This function is compatible with both Seurat v3/v4 (Assay)
#' and v5 (Assay5) objects.
#' @param obj The object to validate.
#' @return A list containing `is_valid` (boolean) and a `message` (string).
validate_seurat_object <- function(obj) {
  # Check 1: Is it a Seurat object?
  if (!is(obj, "Seurat")) {
    return(list(is_valid = FALSE, message = "Error: The uploaded file is not a Seurat object."))
  }
  
  # Check 2: Does the 'RNA' assay exist?
  if (!("RNA" %in% names(obj@assays))) {
      return(list(is_valid = FALSE, message = "Error: The Seurat object must contain an assay named 'RNA'."))
  }
  
  # Check 3: Does the 'RNA' assay contain a raw count matrix?
  # This logic handles both Seurat v5 (Assay5) and older (Assay) objects.
  rna_assay <- obj@assays$RNA
  counts_exist <- FALSE
  
  if (inherits(rna_assay, "Assay5")) {
    # For Seurat v5, check for the 'counts' layer
    if ("counts" %in% names(rna_assay@layers) && nrow(rna_assay@layers$counts) > 0) {
      counts_exist <- TRUE
    }
  } else {
    # For older Seurat versions, check the @counts slot
    if (nrow(rna_assay@counts) > 0) {
      counts_exist <- TRUE
    }
  }
  
  if (!counts_exist) {
    return(list(
      is_valid = FALSE,
      message = "Error: The 'RNA' assay does not contain a raw count matrix. Please ensure it has count data."
    ))
  }
  
  # If all checks pass, return TRUE
  return(list(is_valid = TRUE, message = "Seurat object is valid."))
}


#' Process Seurat Object
#'
#' Runs the standard Seurat workflow. This function is compatible with both
#' Seurat v3/v4 and v5 objects.
#' @param obj The Seurat object to process.
#' @param progress A shiny::Progress object to report progress to the user.
#' @return A fully processed Seurat object.
process_seurat_object <- function(obj, progress) {
  
  # Ensure all subsequent operations are performed on the RNA assay
  DefaultAssay(obj) <- "RNA"
  rna_assay <- obj@assays$RNA
  
  # --- Step 1: Normalization ---
  # Check if the data has been normalized, handling both v5 and older objects.
  is_normalized <- FALSE
  if (inherits(rna_assay, "Assay5")) {
    # For Seurat v5, check for the 'data' layer
    if ("data" %in% names(rna_assay@layers) && nrow(rna_assay@layers$data) > 0) {
      is_normalized <- TRUE
    }
  } else {
    # For older Seurat versions, check the @data slot
    if (nrow(rna_assay@data) > 0) {
      is_normalized <- TRUE
    }
  }

  if (!is_normalized) {
    progress$set(value = 0.5, detail = "Normalizing data...")
    showNotification("Raw counts detected. Normalizing data now.", duration = 5)
    obj <- NormalizeData(obj)
  }
  
  # --- Step 2: Dimensionality Reduction ---
  if (!("umap" %in% names(obj@reductions))) {
    showNotification("UMAP embedding not found. Generating it now. This may take a few minutes.", duration = 10)
    
    progress$set(value = 0.6, detail = "Finding variable features...")
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    
    progress$set(value = 0.7, detail = "Scaling data...")
    all.genes <- rownames(obj)
    obj <- ScaleData(obj, features = all.genes)
    
    progress$set(value = 0.8, detail = "Running PCA...")
    obj <- RunPCA(obj, features = VariableFeatures(object = obj))
    
    progress$set(value = 0.9, detail = "Running UMAP...")
    obj <- RunUMAP(obj, dims = 1:15)
  }
  
  return(obj)
}


#' Perform Differential Gene Expression
#'
#' Runs Seurat's FindMarkers and formats the output.
#' @param obj The Seurat object.
# ... (function is unchanged) ...
perform_dge <- function(obj, grouping_var, group1, group2, progress) {
  
  progress$set(value = 0.5, detail = "Running FindMarkers...")
  
  # Set the identity of the cells based on the user's selection
  Idents(obj) <- obj@meta.data[[grouping_var]]
  
  # Define the second group (NULL means all other cells)
  group2_ident <- if (group2 == "All Other Cells") NULL else group2
  
  # Run FindMarkers
  markers <- FindMarkers(obj, ident.1 = group1, ident.2 = group2_ident, verbose = FALSE)
  
  progress$set(value = 0.9, detail = "Formatting results...")
  
  # Format the results for display
  markers$gene <- rownames(markers)
  markers <- markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
  # Ensure p_val_adj is not zero for log transformation, set to a very small number if it is
  markers$p_val_adj[markers$p_val_adj == 0] <- .Machine$double.xmin
  markers[, 2:6] <- round(markers[, 2:6], 4) # Round numeric columns
  
  return(markers)
}


#' Generate an interactive volcano plot
#'
#' Creates a volcano plot using plotly from DGE results.
#' @param dge_results A data frame from perform_dge.
#' @param fc_threshold Log2 fold-change threshold for significance.
#' @param pval_threshold Adjusted p-value threshold for significance.
#' @return A plotly object.
generate_volcano_plot <- function(dge_results, fc_threshold = 1, pval_threshold = 0.05) {
  
  req(dge_results)
  
  # Create a column to determine significance for coloring based on user-defined thresholds
  dge_results <- dge_results %>%
    mutate(significance = case_when(
      abs(avg_log2FC) > fc_threshold & p_val_adj < pval_threshold ~ "Significant",
      TRUE                                                        ~ "Not Significant"
    ))
  
  # Define colors for the plot: red for significant, grey for not
  colors <- c("Significant" = "#E41A1C", "Not Significant" = "grey")
  
  plot_ly(
    data = dge_results,
    x = ~avg_log2FC,
    y = ~-log10(p_val_adj),
    text = ~paste("Gene:", gene, "<br>log2FC:", round(avg_log2FC, 3), "<br>p.adj:", format.pval(p_val_adj, digits = 3)),
    hoverinfo = "text",
    type = 'scatter',
    mode = 'markers',
    # Use top-level color/colors arguments for a more robust mapping
    color = ~factor(significance, levels = c("Significant", "Not Significant")),
    colors = colors,
    marker = list(
      size = 5,
      opacity = 0.7
    )
  ) %>%
  layout(
    title = "Volcano Plot",
    xaxis = list(title = "Average Log2 Fold Change"),
    yaxis = list(title = "-log10(Adjusted P-value)")
  )
}