#' Apply a Custom Assay Transformation
#'
#' @description
#' Adds a custom per-cell assay transformation step to the HPCell pipeline.
#' The transformation function (e.g. log1p, identity) is applied to each
#' sample's count matrix, optionally capping values at a maximum scale.
#'
#' @param input_hpc An `HPCell` object.
#' @param fx A list of transformation functions (one per sample) to apply to
#'   the assay matrix. Defaults to `identity` for all samples.
#' @param target_input Name of the targets target providing the data object.
#' @param target_output Name of the targets target to write the transformed
#'   data to.
#' @param scale_max A list of numeric upper-bound values (one per sample)
#'   passed to `limit_max_to_scale()`. Default: `10` per sample.
#' @param ... Additional arguments (unused; for method dispatch).
#' @return The updated `HPCell` object with the transformation step appended.
#' @export
transform_assay <- function(input_hpc, fx = input_hpc$initialisation$input_hpc |> map(~identity), 
                            target_input = "data_object", target_output = "sce_transformed", 
                            # make every element to be 20 as a placeholder unless further indicate
                            scale_max = input_hpc$initialisation$input_hpc |> map(~10), ...) {
  UseMethod("transform_assay")
}

#' @rdname transform_assay
#' @importFrom purrr map
#' @export
transform_assay.HPCell = function(
    input_hpc,
    
    # This might be carrying the environment
    fx = input_hpc$initialisation$input_hpc |> map(~"identity"), 
    target_input = "data_object", 
    target_output = "sce_transformed", 
    scale_max = input_hpc$initialisation$input_hpc |> map(~10),
    ...
) {
  
  fx |> saveRDS("temp_fx.rds")
  scale_max |> saveRDS("count_upper_bound.rds")
  
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
    
    # Track the file
    hpc_single("count_upper_bound_file", "count_upper_bound.rds", format = "file") |> 
    hpc_iterate(
      target_output = "count_upper_bound", 
      user_function = readRDS |> quote() ,
      file = "count_upper_bound_file" |> is_target() 
    ) |> 
    
    hpc_iterate(
      target_output = target_output, 
      user_function = transform_utility |> quote() , 
      input_read_RNA_assay = "data_object" |> is_target(), 
      transform_fx = "transform" |> is_target()  ,
      external_path = glue("{input_hpc$initialisation$store}/external") |> as.character(),
      container_type = "data_container_type" |> is_target(),
      scale_max = "count_upper_bound" |> is_target()
      
    )
  
}

#' Limit Maximum Value of Counts by Scaling
#'
#' Scales a numeric vector down so that its maximum value does not exceed
#' \code{scale_max}. If the maximum of \code{counts} is already within the
#' limit, the vector is returned unchanged.
#'
#' @param counts A numeric vector of count values to be scaled.
#' @param scale_max A numeric scalar specifying the upper bound for the maximum
#'   value of \code{counts}.
#'
#' @return A numeric vector of the same length as \code{counts}, with all
#'   values scaled so that the maximum does not exceed \code{scale_max}.
#' @export
limit_max_to_scale <- function(counts, scale_max) {
  if (max(counts) > scale_max) {
    scale_factor <- scale_max / max(counts)
    counts <- counts * scale_factor
  }
  counts
}

#' Safe Exponential Minus One Transformation
#'
#' Applies a two-stage \code{\link[base]{expm1}} transformation to a numeric
#' vector with scaling applied before each stage to prevent overflow. Prior to
#' each call to \code{expm1()}, \code{\link{limit_max_to_scale}} is used to
#' ensure the maximum value does not exceed \code{scale_max}, making the
#' transformation numerically stable for large count values.
#'
#' The transformation pipeline is:
#' \deqn{counts \rightarrow \text{scale} \rightarrow expm1 \rightarrow \text{scale} \rightarrow expm1}
#'
#' @param counts A numeric vector of count values to be transformed.
#' @param scale_max A numeric scalar passed to \code{\link{limit_max_to_scale}}
#'   specifying the upper bound applied before each \code{expm1} step.
#'
#' @return A numeric vector of the same length as \code{counts} with the
#'   two-stage \code{expm1} transformation applied.
#' @seealso \code{\link{expm1}}, \code{\link{limit_max_to_scale}}
#' @export
safe_expm1 <- function(counts, scale_max) {
  counts <- counts |> limit_max_to_scale(scale_max) |> expm1() |> limit_max_to_scale(scale_max) |> expm1()
  counts
}

#' Identity Transformation with Maximum Value Scaling
#'
#' Scales a numeric vector so that its maximum does not exceed \code{scale_max}
#' before returning it unchanged via \code{\link[base]{identity}}. This is
#' useful when no mathematical transformation is desired but the counts still
#' need to be bounded for downstream stability.
#'
#' @param counts A numeric vector of count values to be scaled.
#' @param scale_max A numeric scalar specifying the upper bound for the maximum
#'   value of \code{counts}. Passed directly to \code{\link{limit_max_to_scale}}.
#'
#' @return A numeric vector of the same length as \code{counts}, with all
#'   values scaled so that the maximum does not exceed \code{scale_max}.
#' @seealso \code{\link{limit_max_to_scale}}, \code{\link{safe_expm1}}
#' @export
identity_with_max_limit <- function(counts, scale_max) {
  counts <- counts |> limit_max_to_scale(scale_max) |> identity()
  counts
}

#' Apply a transformation to an assay and save as HDF5
#'
#' This function applies a specified transformation to the assay of a 
#' SummarizedExperiment object and saves the transformed object in HDF5 format.
#'
#' @param input_read_RNA_assay A SummarizedExperiment object to be transformed.
#' @param transform_fx A function to apply to the assay of the SummarizedExperiment object.
#' @param external_path A character string specifying the directory path to save the transformed object.
#' @param container_type A character vector specifying the output file type. Ideally it should match to the input file type.
#' @param scale_max An integer specifying the allowed max count value. It is used to define 
#'     scaling factor in the process of transformation.  
#' @return The function does not return an object. It saves the transformed SummarizedExperiment object to the specified path.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assay<-
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SingleCellExperiment reducedDim<-
#' @importFrom dplyr select
#' @importFrom glue glue
#' @importFrom digest digest
#' @importFrom stats density
#'
#' @export
transform_utility  = function(input_read_RNA_assay, transform_fx, 
                              external_path, container_type, scale_max) {
  
  numer_of_cells_to_sample = 5e3
  
  if(ncol(input_read_RNA_assay) == 0) return(NULL)
  
  # Rename assay names to for consistency
  if (length(names(assays(input_read_RNA_assay))) == 1 && 
      names(assays(input_read_RNA_assay)) != "X") names(assays(input_read_RNA_assay)) <- "X"
  
  # strip metadata that we don't need
  input_read_RNA_assay =
    input_read_RNA_assay |>
    select(.cell, observation_joinid, observation_originalid, donor_id, dataset_id, sample_id, cell_type)

  # Remove reduced dimensions
  reducedDim(input_read_RNA_assay) = NULL
  
  # Remove row data to avoid downstream binding errors
  rowData(input_read_RNA_assay) <- NULL
  
  # Clear memory  
  gc()
  
  dir.create(external_path, showWarnings = FALSE, recursive = TRUE)
  
  # Get the name of the first assay in the data object
  assay_name <- names(assays(input_read_RNA_assay))[1]
  
  # Extract the counts matrix from the assay
  counts <- assay(input_read_RNA_assay, assay_name)
  
  # Convert transform_method to a function if it is a character string
  transform_function <- match.fun(transform_fx)
  
  # Scale counts to a maximum of 20 to avoid downstream failures.  
  # This Check needs ~13Gb to run for 5000+ cell datasets
  # Check if the transformation method is not 'identity' and counts exceed counts upper bound
  # Scale for other transform function is handled internally in `transform_function`
  if (identical(transform_function, identity) || identical(transform_function, expm1)) {
    counts <- limit_max_to_scale(counts, scale_max)
    counts <- transform_function(counts)
  } else {
    # identity_with_max_limit and safe_expm1 handle scale_max internally
    counts <- transform_function(counts, scale_max)
  }
  
  # Clear memory
  gc()
  
  # This is to avoid memory explosion
  set.seed(42)
  counts_light_for_checks = counts[,sample(seq_len(ncol(counts)), size = min(numer_of_cells_to_sample, ncol(counts))),drop=FALSE]
  
  # Compute the density estimate of the counts. This needs ~13Gb to run for 5000+ cell datasets
  density_est <- counts_light_for_checks |> as.matrix() |> density()
  
  # Clear memory
  gc()
  
  # Find the mode (peak) value of the counts
  mode_value <- density_est$x[which.max(density_est$y)]
  
  # If the mode value is negative, shift counts to make the mode zero
  if (mode_value < 0) {
    counts <- counts + abs(mode_value)
    counts_light_for_checks = counts_light_for_checks + abs(mode_value)
  }
  
  # Round counts to avoid potential subtraction errors due to floating-point precision
  counts <- round(counts, 5)
  counts_light_for_checks = round(counts_light_for_checks, 5)
  
  # Find the most frequent count value (mode) in the counts
  majority_gene_counts <- compute_mode_delayedarray(counts_light_for_checks)$mode
  
  # Subtract the mode value from counts if it is not zero
  if (majority_gene_counts != 0) {
    counts <- counts - majority_gene_counts
    counts_light_for_checks <- counts_light_for_checks - majority_gene_counts
  }
  
  # Replace negative counts with zero to avoid downstream failures
  if (min(counts_light_for_checks) < 0) {
    counts[counts < 0] <- 0
  }
  
  # Clear memory
  gc()
  
  # Assign the modified counts back to the data object
  assay(input_read_RNA_assay, assay_name) <- counts
  
  # Remove cells with zero total counts
  # !!! MAYBE WE SHOULD LKEEP THESE CELLS AND LEAVE THEM TO THE FILTERING STEP
  input_read_RNA_assay <- input_read_RNA_assay[, colSums(counts) > 0]
  
  if (ncol(input_read_RNA_assay) == 0) return(NULL)
  
  # Rebuild the SCE to stay light, and to set the assay with the right name
  input_read_RNA_assay = SingleCellExperiment(
    assays = list(X = input_read_RNA_assay |> assay() ), 
    colData = colData(input_read_RNA_assay)
  )
  
  # Return the modified data object
  input_read_RNA_assay |> 
    
    save_experiment_data(
      dir = glue("{external_path}/{digest(input_read_RNA_assay)}"), 
      container_type = container_type
    )
  
  # extension <- switch(container_type,
  #                     
  #                     "sce_rds" = ".rds",
  #                     "seurat_rds" = ".rds",
  #                     
  #                     "seurat_h5" = ".h5Seurat",
  #                     
  #                     "anndata" = ".h5ad",
  #                     
  #                     "sce_hdf5" = "")
  
  # file_name = paste0(file_name, extension)
  # 
  # # Return data as target instead of file_name pointer
  # 
  # input_read_RNA_assay
  # 
  # 
  # extension <- switch(container_type,
  #                     "sce_rds" = ".rds",
  #                     "seurat_rds" = ".rds",
  #                     "seurat_h5" = ".h5Seurat",
  #                     "anndata" = ".h5ad",
  #                     "sce_hdf5" = "")
  # file_name = paste0(file_name, extension)
  
  # Return data as target instead of file_name pointer
  
}
