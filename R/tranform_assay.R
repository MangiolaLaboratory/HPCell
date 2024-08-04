# Define the generic function
#' @export
tranform_assay <- function(input_hpc, fx = "identity", target_input = "read_file", target_output = "sce_transformed", ...) {
  UseMethod("tranform_assay")
}

#' @export
tranform_assay.HPCell = function(input_hpc, fx = identity, target_input = "read_file", target_output = "sce_transformed", ...) {
  
  # Capture all arguments including defaults
  args_list <- as.list(environment())[-1]
  
  # Optionally, you can evaluate the arguments if they are expressions
  args_list <- lapply(args_list, eval, envir = parent.frame())
  
  args_list$factory = function(tiers, target_input, target_output, external_path){
    
    if(length(readRDS("temp_fx.rds"))==1) 
      {
      target_fx = tar_target_raw("transform", readRDS("temp_fx.rds") |> quote(), deployment = "main")
      arguments_to_tier = NULL
    }
    else 
      {
        target_fx = tar_target_raw("transform", readRDS("temp_fx.rds") |> quote(), iteration = "list", deployment = "main")
        arguments_to_tier ="transform"
       }
    
    
    
    list(
      target_fx,
      
      factory_split(
        target_output, 
        i |> 
          read_data_container(container_type = data_container_type) |> 
          transform_utility(transform, e) |> 
          substitute(env = list(i=as.symbol(target_input), e = external_path)),
        tiers, 
        arguments_to_tier = arguments_to_tier,
        other_arguments_to_tier = target_input,
        other_arguments_to_map = target_input,
        iteration = "list", 
        packages = "HPCell"
      )  
    )
    
  }
  
  # We don't want recursive when we call factory
  if(input_hpc |> length() > 0) {

    fx |> saveRDS("temp_fx.rds")
    
    tar_tier_append(
      fx = quote(dummy_hpc |> tranform_assay() %$% tranform_assay %$% factory),
      tiers = input_hpc$initialisation$tier |> get_positions() ,
      target_input = target_input,
      target_output = target_output,
      external_path = glue("{input_hpc$initialisation$store}/external"),
      script = glue("{input_hpc$initialisation$store}.R")
    )
    
  }
  
  
  # Add pipeline step
  input_hpc |>
    c(list(tranform_assay = args_list)) |>
    add_class("HPCell")
  
}

#' Apply a transformation to an assay and save as HDF5
#'
#' This function applies a specified transformation to the assay of a 
#' SummarizedExperiment object and saves the transformed object in HDF5 format.
#'
#' @param i A SummarizedExperiment object to be transformed.
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
transform_utility  = function(i, transform, external_path) {
  #i = i |> read_data_container(container_type = data_container_type) 
  
  dir.create(external_path, showWarnings = FALSE, recursive = TRUE)
  file_name = glue("{external_path}/{digest(i)}")
  
  assay(i) = assay(i) |> transform()
  
  i |> 
    saveHDF5SummarizedExperiment(
      dir = file_name, 
      replace=TRUE, 
      as.sparse=TRUE
    )
  
  file_name
}
