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
#' @param grouping_var The metadata column to group cells by.
#' @param group1 The first group for comparison (test group).
#' @param group2 The second group for comparison (control group). Can be "All Other Cells".
#' @param progress A shiny::Progress object to report progress.
#' @return A formatted data frame with DGE results.
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
  markers[, -1] <- round(markers[, -1], 4) # Round numeric columns
  
  return(markers)
}