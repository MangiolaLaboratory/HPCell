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





#' @importFrom stringr str_extract
#' @export
factory_split = function(
    name_output, command, tiers, arguments_to_tier = c(), arguments_already_tiered = c(), 
    other_arguments_to_map = c(),
    packages = targets::tar_option_get("packages") , 
    deployment = targets::tar_option_get("deployment"), 
    ...
){
  
  if(command |> deparse() |> str_detect("%>%") |> any()) 
    stop("HPCell says: no \"%>%\" allowed in the command, please use \"|>\" ")
  
  #input = command |> deparse() |> paste(collapse = "") |> str_extract("[a-zA-Z0-9_]+\\(([a-zA-Z0-9_]+),.*", group=1) 
  
  # Filter out arguments to be tiered from the input command
  if(arguments_to_tier |> length() > 0)
    arguments_already_tiered <- arguments_already_tiered |> str_subset(paste(arguments_to_tier, collapse = "|"), negate = TRUE)
  
  map2(tiers, names(tiers), ~ {
    
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
      command = command |>  add_tier_inputs(arguments_already_tiered, .y),
      pattern = build_pattern(
        other_arguments_to_map = glue("{other_arguments_to_map}_{.y}"), 
        arguments_to_tier = arguments_to_tier, 
        index = .x
      ) ,
      iteration = "list",
      resources = resources,
      packages = packages, 
      deployment = deployment
    )
  })
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

#' @export
hpc_add_internal = function(
    tiers = NULL, 
    target_output, 
    user_function,
    arguments_to_tier = c(), 
    arguments_already_tiered = c(), 
    other_arguments_to_map = c(), 
    packages = targets::tar_option_get("packages") , 
    deployment = targets::tar_option_get("deployment"),
    ...
){
  
  args <- list(...)  # Capture the ... arguments as a list
  
  # Construct the full call expression with the pipeline substituted into the function
  fx_call <- as.call(c(user_function, args))
  
  if(tiers |> is.null()){
    
      tar_target_raw(
        name = target_output, 
        command = fx_call,
        packages = packages,
        deployment = deployment, 
        iteration = "list", 
        
        # This is in case I am not tiering (e.g. DE analyses) but I need to map
        pattern = build_pattern(other_arguments_to_map = other_arguments_to_map)
      )
      
  }

  
  else 
    list(
      
      factory_split(
        target_output, 
        fx_call,
        tiers, 
        arguments_to_tier = arguments_to_tier,
        arguments_already_tiered = arguments_already_tiered ,
        other_arguments_to_map = other_arguments_to_map,
        packages = packages,
        deployment = deployment,
        ...
      ) 
      # ,
      # 
      # factory_collapse(
      #   "my_report",
      #   bind_rows(o) |> substitute(env = list(o = as.symbol(target_output))),
      #   target_output,
      #   tiers, packages = c("dplyr")
      # )
    )
  
}


#' Add HPC step to pipeline
#'
#' This function adds a new step to the HPC pipeline by appending the appropriate
#' targets to the target script. It allows the user to specify the input and output
#' targets, as well as a custom user function to be applied.
#'
#' @param input_hpc The input HPC object.
#' @param target_output The output target name (default: NULL).
#' @param user_function A custom function provided by the user (default: NULL).
#' @param ... Additional arguments to pass to the internal functions.
#'
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom purrr set_names
#' @export
hpc_iterate = 
  function(input_hpc, target_output = NULL, user_function = NULL, ...) {
    
    # Check for argument consistency
    check_for_name_value_conflicts(...)
    
    # Target script
    target_script = glue("{input_hpc$initialisation$store}.R")
    
    # Delete line with target in case the user execute the command, without calling initialise_hpc
    target_output |>  delete_lines_with_word(target_script)
    
    # please, because sometime we set up list target that do not depend on any other ones
    # if tiers is set to NULL, then the target will not acquire the _<tier> suffix
    # I HAVE TO MAKE THIS MORE ELEGANT, AND NOT RELY ON tiers ARGUMENT
    if(
      list(...) |> arguments_to_action(input_hpc, "tiered") |> length() == 0 &
      list(...) |> arguments_to_action(input_hpc, "tier") |> length() == 0
    ){
      iterate_value = "tier"
      tiers_value = NULL
    } else {
      iterate_value = "tiered"
      tiers_value = input_hpc$initialisation$tier |> get_positions()
    }

    tar_append(
      fx = hpc_add_internal |> quote(),
      tiers = tiers_value ,
      target_output = target_output,
      script = target_script,
      user_function = user_function,
      
      # I HAVE TO IMPROVE the fact that I have to convert to character 
      # because arguments_to_action is also used in expand_tiered_arguments, which needs a named vector
      arguments_to_tier = list(...) |> arguments_to_action(input_hpc, "tier") |> as.character()  , # This "tier" value is decided for each new target below. Usually just at the beginning of the piepline
      arguments_already_tiered = list(...) |> arguments_to_action(input_hpc, "tiered") |> as.character(), # This "tiered" value is decided for each new target below. Ususally every other list targets.
      other_arguments_to_map = list(...) |> arguments_to_action(input_hpc, c("tiered", "map")) |> as.character(), # This "tiered" value is decided for each new target below. Ususally every other list targets.
      ...
    )
  
      
    # Add pipeline step
    input_hpc |>
      c(
        as.list(environment())[-1] |> 
          c(list(iterate = iterate_value)) |> 
          list() |> 
          set_names(target_output) 
      ) |>
      add_class("HPCell")
    
    
  }

#' Add HPC step to pipeline
#'
#' This function adds a new step to the HPC pipeline by appending the appropriate
#' targets to the target script. It allows the user to specify the input and output
#' targets, as well as a custom user function to be applied.
#'
#' @param input_hpc The input HPC object.
#' @param target_output The output target name (default: NULL).
#' @param user_function A custom function provided by the user (default: NULL).
#' @param ... Additional arguments to pass to the internal functions.
#'
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom purrr set_names
#' @export
hpc_single = 
  function(input_hpc, target_output = NULL, user_function = NULL, iterate = "none", ...) {
    
    # tar_script_append2(script = glue("{store}.R"), append = FALSE)
    # 
    # tar_target_raw(target_output, readRDS("input_reference.rds") |> quote()),
    # 
    # # Check for argument consistency
    # check_for_name_value_conflicts(...)
    
    # Target script
    target_script = glue("{input_hpc$initialisation$store}.R")
    
    # Delete line with target in case the user execute the command, without calling initialise_hpc
    target_output |>  delete_lines_with_word(target_script)
    
    tar_append(
      fx = hpc_add_internal |> quote(),
      target_output = target_output,
      script = target_script,
      user_function = user_function,
      ...
    )
    
    
    # Add pipeline step
    input_hpc |>
      c(
        as.list(environment())[-1] |> 
          c(list(iterate = iterate)) |> 
          list() |> 
          set_names(target_output)
      ) |>
      add_class("HPCell")
    
    
  }

#' Add HPC step to pipeline
#'
#' This function adds a new step to the HPC pipeline by appending the appropriate
#' targets to the target script. It allows the user to specify the input and output
#' targets, as well as a custom user function to be applied.
#'
#' @param input_hpc The input HPC object.
#' @param target_output The output target name (default: NULL).
#' @param user_function A custom function provided by the user (default: NULL).
#' @param ... Additional arguments to pass to the internal functions.
#'
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom purrr set_names
#' @export
hpc_merge = 
  function(input_hpc, target_output = NULL, user_function = NULL, ...) {
    
    # Check for argument consistency
    check_for_name_value_conflicts(...)
    
    # Target script
    target_script = glue("{input_hpc$initialisation$store}.R")
    
    # Delete line with target in case the user execute the command, without calling initialise_hpc
    target_output |>  delete_lines_with_word(target_script)
    
    # name_target_intermediate = glue("{target_output}_merge_within_tier")
    
    # tar_append(
    #   fx = hpc_add_internal |> quote(),
    #   tiers = input_hpc$initialisation$tier |> get_positions() ,
    #   target_output = name_target_intermediate,
    #   script = target_script,
    #   user_function = user_function,
    #   arguments_already_tiered = list(...) |> arguments_to_action(input_hpc, "tiered") , # This "tiered" value is decided for each new target below. Ususally every other list targets.
    #   ...
    # )
    
    modified_args = 
      list(...) |> 
      expand_tiered_arguments(
        tiers = input_hpc$initialisation$tier |> get_positions() |> names(), 
        argument_to_replace = list(...) |> arguments_to_action(input_hpc, "tiered") |> names(),
        tiered_args = list(...) |> arguments_to_action(input_hpc, "tiered") |> names()
      )
    
    # this is needed because I cannot use ellipse (...) anymore, I have to use do.call.
    do.call(tar_append, c(
      list(
        fx = hpc_add_internal |> quote() |> quote(),
        #tiers = input_hpc$initialisation$tier |> get_positions(),
        target_output = t |> substitute(env = list(t = target_output)) ,
        script = target_script,
        user_function = u |> quote() |> substitute(env = list(u = user_function))
      ),
      modified_args
    ))
    
    
    
    # Add pipeline step
    input_hpc |>
      c(
        as.list(environment())[-1] |> 
          c(list(iterate = "single")) |> 
          list() |> 
          set_names(target_output) 
      ) |>
      add_class("HPCell")
    
    
  }