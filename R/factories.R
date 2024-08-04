#' Parse a Function Call String
#'
#' This function takes a string representing a function call and parses it into
#' its constituent parts: the function name and the argument names.
#'
#' @param input_string A character string representing a function call, e.g., "report(empty_droplets_tbl, arg1)".
#' 
#' @return A list containing two elements:
#'   \item{function_name}{A character string representing the name of the function.}
#'   \item{arguments}{A character vector containing the names of the arguments.}
#'
#' @importFrom rlang parse_expr
#'
#' @examples
#' # Example usage:
#' command <- quote(report(empty_droplets_tbl, arg1))
#' output <- parse_function_call(command)
#' print(output)
#'
#' @noRd
parse_function_call <- function(command) {
  # Parse the input string to an expression
  expr <- command |> deparse() |> parse_expr()
  
  # Extract the function name
  function_name <- as.character(expr[[1]])
  
  # Extract the arguments
  args <- as.list(expr[-1])
  
  # Convert arguments to character vector
  args <- sapply(args, as.character)
  
  # Create the output list
  result <- list(function_name = function_name, arguments = args)
  
  return(result)
}

#' Parse a Function Call String and Expand Tiered Arguments
#'
#' This function takes a string representing a function call, parses it into
#' its constituent parts (function name and arguments), and expands the specified
#' tiered arguments based on the provided tier labels.
#'
#' @param command A character string representing a function call, e.g., "report(empty_droplets_tbl, arg1)".
#' @param tiers A character vector indicating the tier labels for the specified tiered arguments, e.g., c("_1", "_2", "_3").
#' @param tiered_args A character vector specifying which arguments should be tiered, e.g., c("empty_droplets_tbl").
#' 
#' @return A character string representing the modified function call with tiered arguments expanded.
#'
#' @importFrom rlang parse_expr
#'
#' @examples
#' # Example usage:
#' input_string <- "report(empty_droplets_tbl, arg1)"
#' tiers <- c("_1", "_2", "_3")
#' tiered_args <- c("empty_droplets_tbl", "another_arg")
#' output <- expand_tiered_arguments(input_string, tiers, tiered_args)
#' print(output)
#'
#' @export
expand_tiered_arguments <- function(command, tiers, tiered_args) {
  # Parse the input command to get function name and arguments
  command_character = command |> deparse() 
  
  for(t in tiered_args){
    command_character = command_character |> str_replace(t, paste0(t, "_", tiers) |> paste(collapse = ", "))
  } 
  
  command_character |> rlang::parse_expr()

}



#' @importFrom stringr str_extract
#' @export
factory_split = function(
    name_output, command, tiers, arguments_to_tier = c(), other_arguments_to_tier = c(), 
    other_arguments_to_map = c(), packages = targets::tar_option_get("packages") , ...
  ){
  
  if(command |> deparse() |> str_detect("%>%") |> any()) 
    stop("HPCell says: no \"%>%\" allowed in the command, please use \"|>\" ")
  
  #input = command |> deparse() |> paste(collapse = "") |> str_extract("[a-zA-Z0-9_]+\\(([a-zA-Z0-9_]+),.*", group=1) 
  
  # Filter out arguments to be tiered from the input command
  if(arguments_to_tier |> length() > 0)
    other_arguments_to_tier <- other_arguments_to_tier |> str_subset(paste(arguments_to_tier, collapse = "|"), negate = TRUE)
  
  map2(tiers, names(tiers), ~ {
    
    my_index = .x
    
    # Pattern
    pattern = NULL
    if(
      arguments_to_tier |> length() > 0 |
      other_arguments_to_map |> length() > 0
    ){
      
      pattern = as.name("map")
      
      if(arguments_to_tier |> length() > 0)
        pattern = pattern |> c(
          arguments_to_tier |>
            map(~ substitute(slice(input, index  = arg ), list(input = as.symbol(.x), arg=my_index)) )
        )
      
      if(other_arguments_to_map |> length() > 0)
        pattern = pattern |> c(glue("{other_arguments_to_map}_{.y}") |> lapply(as.name))
      
      pattern = as.call(pattern)
      
    }

    
    
    # Resources
    if(length(tiers) == 1)
      resources = targets::tar_option_get("resources")
    else 
      resources = tar_resources(crew = tar_resources_crew(.y)) 
    
    # # Process additional arguments in ... 
    # additional_args <- list(...)
    # 
    tar_target_raw(
      name = glue("{name_output}_{.y}") |> as.character(), 
      command = command |>  add_tier_inputs(other_arguments_to_tier, .y),
      pattern = pattern,
      iteration = "list",
      resources = resources,
      packages = packages
    )
  })
}

factory_collapse = function(name_output, command, tiered_input, tiers, ...){
  
  command = command |> expand_tiered_arguments(names(tiers), tiered_input)
  
  tar_target_raw(name_output, command, ...) 
}

# factory_tiering = function(preparation, tiering, collapsing, tiers){
#   
#   t1 = tar_target_raw(preparation[[1]], preparation[[2]])
#   t2 = factory_split(
#     tiering[[1]], 
#     tiering[[2]], 
#     tiers,
#     tiering[[3]]
#   )
#   t3 = factory_collapse(
#     collapsing[[1]],
#     collapsing[[2]],
#     collapsing[[3]],
#     tiers
#     
#   )
#   
#   list(t1, t2, t3)
# }


factory_merge_pseudobulk = function(se_list_input, output_se, tiers, external_path){
  
  dir.create(external_path, showWarnings = FALSE, recursive = TRUE)
  
  list(
    factory_split(
      name_output = "pseudobulk_group", 
      command = i |> pseudobulk_merge() |> 
        substitute(env = list(i = as.name(se_list_input))), 
      tiers = tiers, 
      arguments_to_tier = c(), 
      other_arguments_to_tier = c(se_list_input), 
      packages = c("tidySummarizedExperiment")
    ) ,
    
    factory_collapse(
      name_output = output_se, 
      command =
        cbind(pseudobulk_group) |>
        
        # Conver to monolytic H5 cebause otherwise the memory blows up
        saveHDF5SummarizedExperiment(dir = f, replace=TRUE, as.sparse=TRUE) |> 
        substitute(env = list(f = glue("{external_path}/output_se"))), 
      tiers = tiers, 
      tiered_input = "pseudobulk_group",
      packages = c("HDF5Array")
    )
  )
}



#' @importFrom tidybulk test_differential_abundance
#' 
#' 
#' @export
factory_de_fix_effect = function(se_list_input, output_se, formula, method, tiers, factor_of_interest = NULL, .abundance = NULL){
  list(
    
    
    tar_target_raw("chunk_tbl", 
                   pseudobulk_se |> 
                     rownames() |> 
                     feature_chunks() |> 
                     quote()
                   
    ),
    
    tar_target_raw(
      "pseudobulk_group_list",
      pseudobulk_se |> 
        group_split(!!sym(pseudobulk_group_by)) |>  
        quote(),
      iteration = "list",
      packages = c("tidySummarizedExperiment", "S4Vectors", "targets")
    ),
    
    
    # Analyse
    tar_target_raw(
      output_se,
      pseudobulk_group_list |> 
        keep_abundant(factor_of_interest = i, .abundance = a) |> 
        test_differential_abundance(
          f,
          .abundance = !!sym(a),
          method = m
        ) |> 
        pivot_transcript() |> 
        substitute(
          env = list(
            f=as.formula(formula), a = .abundance, m=method, 
            i = factor_of_interest
          )),
      pattern = map(pseudobulk_group_list) |> quote(), 
      packages="tidybulk"
    )
    
    
    # factory_collapse(
    #   "pseudobulk_joined", 
    #   command = bind_rows(create_pseudobulk_sample) |> left_join(chunk_tbl) |> quote(),
    #   tiered_input = "create_pseudobulk_sample", 
    #   tiers = tiers, 
    #   packages = c("dplyr") #, pattern = map(create_pseudobulk_sample) |> quote()
    # )
  )
}

#' @importFrom rlang sym
#' @importFrom dplyr left_join
#' @importFrom dplyr nest
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr unnest
#' @importFrom purrr map
#' @importFrom S4Vectors split
#' @importFrom purrr compact
#' @importFrom purrr map
#' 
#' @export
factory_de_random_effect = function(se_list_input, output_se, formula, tiers, factor_of_interest = NULL, .abundance = NULL){
  
  
  list(
    
    tar_target_raw(
      "pseudobulk_group_list",
      pseudobulk_se |> 
        group_split(!!sym(pseudobulk_group_by)) |>  
        quote(),
      iteration = "list",
      packages = c("tidySummarizedExperiment", "S4Vectors", "targets", "rlang")
    ),
    
    tar_target_raw(
      "pseudobulk_table_dispersion_gene",
      pseudobulk_group_list  |> 
        keep_abundant(factor_of_interest = i, .abundance = a, minimum_proportion = 0.2, minimum_counts = 1,) |> 
        scale_abundance() |> 
        se_add_dispersion(f, a ) |> 
        split_summarized_experiment(chunk_size = 10) |> 
        substitute(env = list(
          f=as.formula(formula), 
          a = glue("{.abundance}_scaled"),
          i = factor_of_interest 
          )), 
      pattern = map(pseudobulk_group_list) |> quote(), 
      iteration = "list",
      packages = c("tidySummarizedExperiment", "S4Vectors", "purrr", "dplyr" ,"tidybulk", "HPCell")
    ), 
    
   tar_target_raw(
     "pseudobulk_table_dispersion_gene_unlist", 
     pseudobulk_table_dispersion_gene |> unlist() |> unlist() |> quote(),
     iteration = "list"
     #, 
     #pattern = map(pseudobulk_table_dispersion_gene) |> quote()
    ),
    
    # Analyse
    tar_target_raw(
      "de_split",
      pseudobulk_table_dispersion_gene_unlist |> 
        map_de( f, a, 
          "glmmseq_lme4", 
          max_rows_for_matrix_multiplication = 10000, 
          cores = 1,
          .scaling_factor = !!sym(s),
          .group_by = !!sym(pseudobulk_group_by)
        ) |> 

        substitute(env = list(f=as.formula(formula), a = .abundance, s = "TMM")),
      pattern = map(pseudobulk_table_dispersion_gene_unlist) |> quote(),
      iteration = "list",
      packages = c("rlang")
    ),
   
   tar_target_raw(output_se, 
                  de_split |>
                    compact() |> 
                    bind_rows() |>
                    nest(summary = -!!sym(pseudobulk_group_by)) |>
                    quote(), 
                  packages = c("tidyr", "dplyr", "purrr")
                )
    
  )
}
