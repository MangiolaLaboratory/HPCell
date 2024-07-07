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
  parsed_call <- parse_function_call(command)
  
  function_name <- parsed_call$function_name
  arguments <- parsed_call$arguments
  
  # Expand tiered arguments
  expanded_arguments <- unlist(lapply(arguments, function(arg) {
    if (arg %in% tiered_args) {
      # Generate tiered arguments using provided tier labels
      return(paste0(arg, "_", tiers))
    } else {
      # Return the argument as is
      return(arg)
    }
  }))
  
  # Construct the new function call string
  expanded_call <- paste0(function_name, "(", paste(expanded_arguments, collapse = ", "), ")")
  
  expanded_call |> rlang::parse_expr()
  
}



#' @importFrom stringr str_extract
factory_split = function(name_output, command, tiers){
  
  if(command |> deparse() |> str_detect("%>%")) 
    stop("HPCell says: no \"%>%\" allowed in the command, please use \"|>\" ")
  
  input = command |> deparse() |> str_extract("[a-zA-Z0-9_]+\\((.+),.*", group=1) |> as.symbol() 
  
  tiers |> imap(~ {
    
    pattern = substitute(slice(input, index  = arg ), list(input = input, arg=.x))
    if(length(tiers) == 1)
      resources = targets::tar_option_get("resources")
    else 
      resources = substitute(tar_resources(crew = tar_resources_crew(arg)) , list(arg = .x))
    
    
    tar_target_raw(
      glue("{name_output}_{.y}"),
      command,
      pattern = pattern,
      iteration = "list",
      resources = resources
    )
  })
}

factory_collapse = function(name_output, command, tiered_input, tiers){
  
  command = command |> expand_tiered_arguments(names(tiers), tiered_input)
  
  tar_target_raw(name_output, command) 
}

factory_tiering = function(preparation, tiering, collapsing, tiers){
  
  t1 = tar_target_raw(preparation[[1]], preparation[[2]])
  t2 = factory_split(
    tiering[[1]], 
    tiering[[2]], 
    tiers
  )
  t3 = factory_collapse(
    collapsing[[1]],
    collapsing[[2]],
    collapsing[[3]],
    tiers
    
  )
  
  list(t1, t2, t3)
}
