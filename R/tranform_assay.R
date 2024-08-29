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
    fx = input_hpc$initialisation$input_hpc |> map(~identity), 
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
    
    # # Load data container type
    # hpc_single("data_container_type_file", "data_container_type.rds", format = "file") |>
    # 
    # hpc_single(
    #   target_output = "data_type",
    #   user_function = readRDS |> quote(),
    #   file = "data_container_type_file" |> is_target()
    # ) |>
    
    hpc_iterate(
      target_output = target_output, 
      user_function = transform_utility |> quote() , 
      input_read_RNA_assay = as.name(target_input), 
      transform_fx = transform |> quote() ,
      external_path = glue("{input_hpc$initialisation$store}/external"),
      data_container_type = input_hpc$initialisation$data_container_type
    )
  
}

#' Apply a transformation to an assay and save as HDF5
#'
#' This function applies a specified transformation to the assay of a 
#' SummarizedExperiment object and saves the transformed object in HDF5 format.
#'
#' @param input_read_RNA_assay A SummarizedExperiment object to be transformed.
#' @param transform_fx A function to apply to the assay of the SummarizedExperiment object.
#' @param external_path A character string specifying the directory path to save the transformed object.
#' @param data_container_type A character vector specifying the output file type. Ideally it should match to the input file type.
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
  input_read_RNA_assay = input_read_RNA_assay |> transform_fx()
  
  input_read_RNA_assay |> 
    save_experiment_data(dir = file_name,
                         container_type = data_container_type )
    # saveHDF5SummarizedExperiment(
    #   dir = file_name,
    #   replace=TRUE,
    #   as.sparse=TRUE
    # )
  
  extension <- switch(data_container_type,
                      "sce_rds" = ".rds",
                      "seurat_rds" = ".rds",
                      "seurat_h5" = ".h5Seurat",
                      "anndata" = ".h5ad",
                      "sce_hdf5" = "")
  file_name = paste0(file_name, extension)
  file_name
}
  
  