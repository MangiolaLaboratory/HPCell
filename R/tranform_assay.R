# Define the generic function
#' @export
tranform_assay <- function(input_hpc, fx = input_hpc$initialisation$input_hpc |> map(~identity), target_input = "data_object", target_output = "sce_transformed", ...) {
  UseMethod("tranform_assay")
}

#' @importFrom purrr map
#' 
#' @export
tranform_assay.HPCell = function(
    input_hpc,
    
    # This might be carrying the environment
    fx = input_hpc$initialisation$input_hpc |> map(~"identity"), 
    target_input = "data_object", 
    target_output = "sce_transformed", 
    ...
) {
  
  fx |> saveRDS("temp_fx.rds")
  
  input_hpc |> 
    
    # Track the file
    hpc_single("transform_file", "temp_fx.rds", format = "file") |> 
    
    hpc_iterate(
      target_output = "transform", 
      user_function = readRDS |> quote() ,
      file = "transform_file" |> is_target() 
      # ,
      # iteration = "list", 
      # deployment = "main"
    ) |> 
    
    hpc_iterate(
      target_output = target_output, 
      user_function = transform_utility |> quote() , 
      input_read_RNA_assay = "target_input" |> is_target(), 
      transform_fx = "transform" |> is_target()  ,
      external_path = glue("{input_hpc$initialisation$store}/external") |> as.character()
    )
  
}

#' Harmonize Counts Data in a SingleCellExperiment Object
#'
#' This function harmonizes the counts data in a \code{SingleCellExperiment} object by adjusting negative values,
#' scaling counts to avoid downstream failures, applying a transformation method, and removing cells with zero counts.
#'
#' @param data A \code{SingleCellExperiment} object containing assays with counts data.
#' @param transform_method A function or the name of a function (as a character string) to transform the counts data (e.g., \code{"log1p"}, \code{"exp"}).
#'
#' @return The modified \code{SingleCellExperiment} object with harmonized counts.
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' # Using a function object
#' transformed_sce <- census_harmonise_anndata_counts(sce, log1p)
#' # Using a function name as a character string
#' transformed_sce <- census_harmonise_anndata_counts(sce, "exp")
#' }
#'
#' @export
census_harmonise_anndata_counts <- function(data, transform_method) {
  # Convert transform_method to a function if it is a character string
  transform_function <- match.fun(transform_method)
  
  # Get the name of the first assay in the data object
  assay_name <- names(assays(data))[1]
  
  # Extract the counts matrix from the assay
  counts <- assay(data, assay_name)
  
  # Compute the density estimate of the counts
  density_est <- density(as.matrix(counts))
  
  # Find the mode (peak) value of the counts
  mode_value <- density_est$x[which.max(density_est$y)]
  
  # If the mode value is negative, shift counts to make the mode zero
  if (mode_value < 0) {
    counts <- counts + abs(mode_value)
  }
  
  # Scale counts to a maximum of 20 to avoid downstream failures
  # Check if the transformation method is not 'identity' and counts exceed 20
  if (!identical(transform_function, identity) && (max(counts) > 20)) {
    scale_factor <- 20 / max(counts)
    counts <- counts * scale_factor
  }
  
  # Apply the transformation method to counts
  counts <- transform_function(counts)
  
  # Round counts to avoid potential subtraction errors due to floating-point precision
  counts <- round(counts, 5)
  
  # Find the most frequent count value (mode) in the counts
  majority_gene_counts <- as.numeric(names(which.max(table(as.vector(counts)))))
  
  # Subtract the mode value from counts if it is not zero
  if (majority_gene_counts != 0) {
    counts <- counts - majority_gene_counts
  }
  
  # Replace negative counts with zero to avoid downstream failures
  if (min(counts[, seq_len(min(10000, ncol(counts)))]) < 0) {
    counts[counts < 0] <- 0
  }
  
  # Assign the modified counts back to the data object
  assay(data, assay_name) <- counts
  
  # Calculate the column sums (total counts per cell)
  col_sums <- colSums(counts)
  
  # Remove cells with zero total counts
  data <- data[, col_sums > 0]
  
  # Remove row data to avoid downstream binding errors
  rowData(data) <- NULL
  
  # Return the modified data object
  data
}

#' Apply a transformation to an assay and save as HDF5
#'
#' This function applies a specified transformation to the assay of a 
#' SummarizedExperiment object and saves the transformed object in HDF5 format.
#'
#' @param input_read_RNA_assay A SummarizedExperiment object to be transformed.
#' @param transform A function to apply to the assay of the SummarizedExperiment object.
#' @param external_path A character string specifying the directory path to save the transformed object.
#'
#' @return The function does not return an object. It saves the transformed SummarizedExperiment object to the specified path.
#'
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom glue glue
#' @importFrom tools digest
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#'
#' @export
transform_utility  = function(input_read_RNA_assay, transform_fx, external_path, data_container_type) {
  
  #input_read_RNA_assay = input_read_RNA_assay |> read_data_container(container_type = data_container_type) 
  
  dir.create(external_path, showWarnings = FALSE, recursive = TRUE)
  
  file_name = glue("{external_path}/{digest(input_read_RNA_assay)}")
  
  #assay(input_read_RNA_assay) = assay(input_read_RNA_assay) |> transform_fx()
  
  input_read_RNA_assay = input_read_RNA_assay |> census_harmonise_anndata_counts(transform_fx)
  
  input_read_RNA_assay |> 
    
    save_experiment_data(dir = file_name,
                         
                         container_type = data_container_type )
  
  extension <- switch(data_container_type,
                      
                      "sce_rds" = ".rds",
                      
                      "seurat_rds" = ".rds",
                      
                      "seurat_h5" = ".h5Seurat",
                      
                      "anndata" = ".h5ad",
                      
                      "sce_hdf5" = "")
  
  file_name = paste0(file_name, extension)
  
  # Return data as target instead of file_name pointer
  
  input_read_RNA_assay
  
}

