# R/utils.R
# This script contains helper functions for validating and processing the Seurat object.
# By keeping them here, we can reuse them in other modules and keep our server logic clean.

#' Validate Seurat Object
#'
#' Checks if the uploaded object is a valid Seurat object with raw count data.
#' @param obj The object to validate.
#' @return A list containing `is_valid` (boolean) and a `message` (string).
validate_seurat_object <- function(obj) {
  # Check 1: Is it a Seurat object?
  if (!is(obj, "Seurat")) {
    return(list(is_valid = FALSE, message = "Error: The uploaded file is not a Seurat object."))
  }
  
  # Check 2: Does it have a default assay?
  if (length(obj@assays) == 0) {
      return(list(is_valid = FALSE, message = "Error: The Seurat object does not contain any assays."))
  }
  
  # Check 3: Does it have a raw count matrix?
  if (nrow(obj@assays[[DefaultAssay(obj)]]@counts) == 0) {
    return(list(
      is_valid = FALSE,
      message = "Error: The Seurat object's default assay does not contain a raw count matrix (`counts` slot is empty)."
    ))
  }
  
  # If all checks pass, return TRUE
  return(list(is_valid = TRUE, message = "Seurat object is valid."))
}


#' Process Seurat Object
#'
#' Runs the standard Seurat workflow (Normalize, Scale, PCA, UMAP) if the steps
#' have not already been completed.
#' @param obj The Seurat object to process.
#' @param progress A shiny::Progress object to report progress to the user.
#' @return A fully processed Seurat object.
process_seurat_object <- function(obj, progress) {
  
  # --- Step 1: Normalization ---
  # Check if the data slot is empty, which implies it hasn't been normalized.
  if (nrow(obj@assays[[DefaultAssay(obj)]]@data) == 0) {
    progress$set(value = 0.5, detail = "Normalizing data...")
    showNotification("Raw counts detected. Normalizing data now.", duration = 5)
    obj <- NormalizeData(obj)
  }
  
  # --- Step 2: Dimensionality Reduction ---
  # Check if UMAP has been run. If not, run the full pipeline.
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