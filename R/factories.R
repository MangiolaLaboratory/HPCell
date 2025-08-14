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


#' @export
hpc_internal = function(
    tiers = NULL, 
    target_output, 
    user_function,
    arguments_to_tier = c(), 
    arguments_already_tiered = c(), 
    other_arguments_to_map = c(), 
    packages = targets::tar_option_get("packages") , 
    deployment = targets::tar_option_get("deployment"),
    format = targets::tar_option_get("format"),
    ...
){
  
  args <- list(...)  # Capture the ... arguments as a list
  
  
  # If format is file just pass the argument
  if(format != "file")
    
    # Construct the full call expression with the pipeline substituted into the function
    user_function <- as.call(c(user_function, args))
  
  if(tiers |> is.null() || tiers |> length() < 2){
    
      tar_target_raw(
        name = target_output |> as.character(), 
        command = user_function,
        
        # This is in case I am not tiering (e.g. DE analyses) but I need to map
        pattern = build_pattern(other_arguments_to_map = other_arguments_to_map),
        
        iteration = "list", 
        packages = packages,
        deployment = deployment,
        format = format


      )
      
  }

  
  else {
    
    if(user_function |> deparse() |> str_detect("%>%") |> any()) 
      stop("HPCell says: no \"%>%\" allowed in the command, please use \"|>\" ")
    
    # Filter out arguments to be tiered from the input command
    if(arguments_to_tier |> length() > 0)
      arguments_already_tiered <- arguments_already_tiered |> str_subset(paste(arguments_to_tier, collapse = "|"), negate = TRUE)
    
    map2(tiers, names(tiers), ~ {
    
      tar_target_raw(
        name = 
          glue("{target_output}_{.y}") |> 
          
          # This is needed because using glue
          as.character() , 
        command = user_function |>  add_tier_inputs(arguments_already_tiered, .y),
        pattern = 
          build_pattern(
            other_arguments_to_map = glue("{other_arguments_to_map}_{.y}"), 
            arguments_to_tier = arguments_to_tier, 
            index = .x
          ) ,
        iteration = "list",
        packages = packages, 
        deployment = deployment,
        resources = tar_resources(crew = tar_resources_crew(.y)) ,
        format = format
      )
    })
    
  }
   
}

#' @export
hpc_internal_report = function(
    tiers = NULL, 
    target_output, 
    rmd_path,
    render_arguments = quote(list()),
    output_file = NULL,
    arguments_to_tier = c(), 
    arguments_already_tiered = c(), 
    other_arguments_to_map = c(), 
    packages = targets::tar_option_get("packages") , 
    deployment = targets::tar_option_get("deployment"),
    ...
){
  

  if(tiers |> is.null() || tiers |> length() < 2){
    
    tar_quarto_raw(
      name = target_output |> as.character(), 
      path = rmd_path,
      output_file = output_file,
      execute_params = render_arguments,
      # This is in case I am not tiering (e.g. DE analyses) but I need to map
      # pattern = build_pattern(other_arguments_to_map = other_arguments_to_map),
      
      # iteration = "list", 
      packages = packages,
      deployment = deployment
      
      
    )
    
  }
  
  
  else {
    
    # Filter out arguments to be tiered from the input command
    if(arguments_to_tier |> length() > 0)
      arguments_already_tiered <- arguments_already_tiered |> str_subset(paste(arguments_to_tier, collapse = "|"), negate = TRUE)
    
    map2(tiers, names(tiers), ~ {
      
      tar_quarto_raw(
        name = 
          glue("{target_output}_{.y}") |> 
          
          # This is needed because using glue
          as.character() , 
        pattern = 
          build_pattern(
            other_arguments_to_map = glue("{other_arguments_to_map}_{.y}"), 
            arguments_to_tier = arguments_to_tier, 
            index = .x
          ) ,
        # iteration = "list",
        packages = packages, 
        deployment = deployment,
        resources = tar_resources(crew = tar_resources_crew(.y)) 
      )
    })
    
  }
  
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
  function(
    input_hpc, 
    target_output = NULL, 
    user_function = NULL, 
    user_function_source_path = NULL,
    ...
  ) {
    
    # Check for argument consistency
    check_for_name_value_conflicts(...)
    
    # Target script
    target_script = glue("{input_hpc$initialisation$store}.R")
    
    # Delete line with target in case the user execute the command, without calling initialise_hpc
    target_output |>  delete_lines_with_word(target_script)
    
    # Append source if any
    write_source(user_function_source_path, target_script)
      
    # please, because sometime we set up list target that do not depend on any other ones
    # if tiers is set to NULL, then the target will not acquire the _<tier> suffix
    # I HAVE TO MAKE THIS MORE ELEGANT, AND NOT RELY ON tiers ARGUMENT
    if(input_hpc$initialisation$tier |> get_positions() |> length() < 2){
      iterate_value = "map"
      tiers_value = NULL
    }
      
    else if(
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
      fx = hpc_internal |> quote(),
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
  function(
    input_hpc, 
    target_output = NULL, 
    user_function = NULL, 
    user_function_source_path = NULL,
    iterate = "none", 
    ...) {
    
    # Target script
    target_script = glue("{input_hpc$initialisation$store}.R")
    
    # Delete line with target in case the user execute the command, without calling initialise_hpc
    target_output |>  delete_lines_with_word(target_script)
    
    # Append source if any
    write_source(user_function_source_path, target_script)
    
    
    tar_append(
      fx = hpc_internal |> quote(),
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
  function(
    input_hpc, 
    target_output = NULL, 
    user_function = NULL, 
    user_function_source_path = NULL,
    ...
  ) {
    
    # Check for argument consistency
    check_for_name_value_conflicts(...)
    
    # Target script
    target_script = glue("{input_hpc$initialisation$store}.R")
    
    # Delete line with target in case the user execute the command, without calling initialise_hpc
    target_output |>  delete_lines_with_word(target_script)
    
    # Append source if any
    write_source(user_function_source_path, target_script)
    
    
    # If no tiers
    if(input_hpc$initialisation$tier |> get_positions() |> length() < 2)
      tar_append(
          fx = hpc_internal |> quote(),
          target_output = target_output,
          script = target_script,
          user_function = user_function,
          ...
      )
      
    else{
      
      args = 
        list(...)  |> 
        expand_tiered_arguments(
          tiers = input_hpc$initialisation$tier |> get_positions() |> names(), 
          argument_to_replace = list(...) |> arguments_to_action(input_hpc, "tiered") |> names(),
          tiered_args = list(...) |> arguments_to_action(input_hpc, "tiered") |> names()
        )
      
      # this is needed because I cannot use ellipse (...) anymore, I have to use do.call.
      do.call(tar_append, c(
        list(
          fx = hpc_internal |> quote() |> quote(),
          #tiers = input_hpc$initialisation$tier |> get_positions(),
          target_output = t |> substitute(env = list(t = target_output)) ,
          script = target_script,
          user_function = u |> quote() |> substitute(env = list(u = user_function))
        ),
        args
      ))
    }

    
    
    
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
#' @importFrom here here
#' 
#' 
#' @export
hpc_report = function(input_hpc, target_output = NULL, rmd_path = NULL, ...) {
    
    # # Check for argument consistency
    # check_for_name_value_conflicts(...)
    # 
    # Target script
    target_script = glue("{input_hpc$initialisation$store}.R")
    
    # Delete line with target in case the user execute the command, without calling initialise_hpc
    target_output |>  delete_lines_with_word(target_script)
    
    external_dir = glue("{input_hpc$initialisation$store}/external") |> as.character() |> here()
    dir.create(external_dir, showWarnings = FALSE, recursive = TRUE)
    
    # If no tiers
    if(input_hpc$initialisation$tier |> get_positions() |> length() < 2)
      tar_append(
        fx = hpc_internal_report |> quote(),
        target_output = target_output,
        script = target_script,
        rmd_path = rmd_path,
        output_file = glue("{external_dir}/{target_output}") |> as.character(),
        render_arguments = substitute(quote(expr), list(expr = list(params = list(...)))) # Add quotation
      )
    
    else{
      
      args = 
        list(...)  |> 
        expand_tiered_arguments(
          tiers = input_hpc$initialisation$tier |> get_positions() |> names(), 
          argument_to_replace = list(...) |> arguments_to_action(input_hpc, "tiered") |> names(),
          tiered_args = list(...) |> arguments_to_action(input_hpc, "tiered") |> names()
        )
      
      # this is needed because I cannot use ellipse (...) anymore, I have to use do.call.
      do.call(tar_append, c(
        list(
          fx = hpc_internal |> quote() |> quote(),
          #tiers = input_hpc$initialisation$tier |> get_positions(),
          target_output = t |> substitute(env = list(t = target_output)) ,
          script = target_script,
          user_function = u |> quote() |> substitute(env = list(u = user_function))
        ),
        args
      ))
    }
    
    
    
    
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