# Helper function to add class to an object
add_class <- function(obj, class_name) {
  class(obj) <- c(class_name, class(obj))
  return(obj)
}

# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

# Raise to the power
pow = function(a,b){	a^b }

# Equals
eq = function(a,b){	a==b }

#' Read various types of single-cell data
#' @param file A character vector of length one specifies the file path, or directory path.
#'   For data format anndata, rds and seurat_h5, use file path.
#'   For data format hdf5, use directory path.
#' @param container_type A character vector of length one specifies the input data type.
#' @return A `[Seurat::Seurat-class]` object
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
#' @export
read_data_container <- function(file,
                                container_type = "anndata"){
  
  if (container_type == "seurat_h5") {
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      stop("HPCell says: You need to install the SeuratDisk package.")
    }
  }
  
  if (container_type == "anndata") {
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
      stop("HPCell says: You need to install the zellkonverter package.")
    }
  }
  
  switch(container_type,
         "anndata" = zellkonverter::readH5AD(file, reader = "R", use_hdf5 = TRUE, 
                                             obs = FALSE, raw = FALSE, layers = FALSE),
         "sce_rds" = readRDS(file),
         "seurat_rds" = readRDS(file),
         "sce_hdf5" = loadHDF5SummarizedExperiment(file),
         "seurat_h5" = SeuratDisk::LoadH5Seurat(file)
  )
}

#' Save various types of single-cell data
#' @param data A data object to save.
#' @param dir A character vector of length one specifies the file path, or directory path.
#' @param container_type A character vector of length one specifies the input data type.
#' @return An object stored in the defined path.
#' @importFrom HDF5Array loadHDF5SummarizedExperiment saveHDF5SummarizedExperiment
#' @importFrom SummarizedExperiment assay
#' @export
save_experiment_data <- function(data,
                                 dir,
                                 container_type = "anndata"){
  
  if (container_type == "seurat_h5") {
    if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
      stop("HPCell says: You need to install the SeuratDisk package.")
    }
  }
  
  if (container_type == "anndata") {
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
      stop("HPCell says: You need to install the zellkonverter package.")
    }
  }
  
  switch(container_type,
         "anndata" = {
           if (ncol(assay(data)) == 1) data = data |> duplicate_single_column_assay()
           zellkonverter::writeH5AD(data, paste0(dir, ".h5ad"), compression = "gzip") 
           read_data_container(paste0(dir, ".h5ad"), "anndata")
         },
         "sce_rds" = data,
         "seurat_rds" = data,
         "sce_hdf5" = saveHDF5SummarizedExperiment(data,
                                                   dir,
                                                   replace = TRUE,
                                                   as.sparse = TRUE),
         
         "seurat_h5" = SeuratDisk::SaveH5Seurat(data,
                                                paste0(dir, ".h5Seurat"),
                                                overwrite = TRUE)
  )
}

#' Duplicate Single-Column Assay in a SingleCellExperiment or Seurat Object
#'
#' This function handles a `SingleCellExperiment` or `Seurat` object where a specified assay 
#' contains only one column. It duplicates the single-column assay to avoid potential 
#' errors during saving or downstream analysis that require at least two columns. 
#' The duplicated column is marked with a prefix `DUMMY___` to distinguish it. 
#' Corresponding entries in the column metadata (`colData`) are also duplicated.
#'
#' @param data A `SingleCellExperiment` or `Seurat` object.
#' @importFrom SummarizedExperiment assay assays colData
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Seurat GetAssayData CreateSeuratObject
#' @importFrom rlang set_names
#' @return A modified `SingleCellExperiment` or `Seurat` object with the single-column assay 
#' duplicated if applicable.
duplicate_single_column_assay <- function(data) {
  
  assay_name = data@assays|> names() |> magrittr::extract2(1)
  
  if (inherits(data, "SingleCellExperiment") && ncol(data) == 1){
    
    # Duplicate the assay to prevent saving errors due to single-column matrices
    my_assay = cbind(assay(data), assay(data))
    # Rename the second column to distinguish it
    colnames(my_assay)[2] = paste0("DUMMY", "___", colnames(my_assay)[2])
    
    cd = colData(data)
    cd = cd |> rbind(cd)
    rownames(cd)[2] = paste0("DUMMY", "___", rownames(cd)[2])
    
    data =  SingleCellExperiment(assay = list(my_assay) |> set_names(assay_name), colData = cd)
    data
  } 
  
  if (inherits(data, "Seurat") && ncol(data) == 1) {
    
    my_assay <- GetAssayData(data, layer = "counts", assay = assay_name)
    my_assay <- cbind(my_assay, my_assay)
    colnames(my_assay)[2] <- paste0("DUMMY___", colnames(my_assay)[2])
    
    cd <- data[[]]
    cd <- rbind(cd, cd)
    rownames(cd)[2] <- paste0("DUMMY___", rownames(cd)[2])
    
    data <- CreateSeuratObject(counts = my_assay, meta.data = cd, assay = assay_name)
    data
  }
  data
}

#' Gene name conversion using ensembl database
#' 
#' @param id Character vector of gene names
#' @param current_nomenclature Character vector of input gene nomenclature
#' @return A data frame of gene name before and after conversion
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importMethodsFrom ensembldb genes   
#' @export
convert_gene_names <- function(id,
                               current_nomenclature) {
  if (current_nomenclature == "symbol"){
    edb <- EnsDb.Hsapiens.v86
    edb_df <- genes(edb,
                    columns = c("gene_name", "entrezid", "gene_biotype"),
                    filter = AnnotationFilter::GeneNameFilter(id),
                    return.type = "data.frame")   
  } 
  else if (current_nomenclature == "ensembl") {
    edb <- EnsDb.Hsapiens.v86
    edb_df <- genes(edb,
                    columns = c("gene_name", "entrezid", "gene_biotype"),
                    filter = AnnotationFilter::GeneIdFilter(id),
                    return.type = "data.frame")   
  }
  edb_df
}

#' Transform counts to continous data
#' @param counts A SummarizedExperiment object
#' @importFrom tidyr pivot_longer
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom magrittr extract2
get_count_per_gene_df <- function(counts) {
  #assay_name = data@assays |> names() |> extract2(1)
  #counts <- assay(data, assay_name) |> as.data.frame() |> rownames_to_column(var = "features")
  counts_tidy <- counts |> as.data.frame() |> tibble::rownames_to_column(var = "features") |>
    as_tibble() |> pivot_longer(!features, names_to = "cells",
                                values_to = "counts")
  counts_tidy
}







#' Reference Label Fine Identification
#'
#' @description
#' Adds a reference label depending on the provided tissue type of the input data.
#'
#' @param tissue Type of tissue.
#'
#' @return Appropriate reference label for fine categorization.
#' @export
reference_label_fine_id <- function(tissue) {
  return(
    ifelse(tissue == "pbmc", "monaco_first.labels.fine",
           ifelse(tissue == "solid", "blueprint_first.labels.fine",
                  ifelse(tissue == "atypical", "none",
                         ifelse(tissue == "none", "monaco_first.labels.fine", NA)))))
}


#' Reference Label Coarse Identification
#'
#' @description
#' Adds a reference label depending on the provided tissue type of the input data.
#'
#' @param tissue Type of tissue.
#'
#' @return Appropriate reference label for coarse categorization.
#' @export
reference_label_coarse_id <- function(tissue) {
  return(
    ifelse(tissue == "pbmc", "monaco_first.labels.coarse",
           ifelse(tissue == "solid", "blueprint_first.labels.coarse",
                  ifelse(tissue == "atypical", "none",
                         ifelse(tissue == "none", "monaco_first.labels.coarse", NA)))))
}

#' Change Default Assay to RNA
#'
#' @importFrom Seurat DefaultAssay
#' @description
#' `add_RNA_assay` changes the default assay in a Seurat object to RNA.
#'
#' @param input_read Seurat object.
#' @param RNA_assay_name Name of the RNA assay to set as default.
#'
#' @return Modified Seurat object with the default assay set to RNA.
#'
#' @importFrom Seurat DefaultAssay
#' @noRd
add_RNA_assay <- function(input_read, RNA_assay_name){
  
  
  if (RNA_assay_name != "RNA"){
    input_read[["RNA"]] = input_read[[RNA_assay_name]]
    Seurat::DefaultAssay(object = input_read) <- "RNA"
    input_read[[RNA_assay_name]] = NULL
  }
  
  # names(input_read@assays)<- names(input_read@assays) |> sapply(function(x) if(x == RNA_assay_name) "RNA" else x)
  input_read
}


#' Identify Variable Features by Cell Type in Seurat
#'
#' @description
#' Extracts variable features for each cell type within a Seurat object.
#'
#' @param counts Seurat object.
#' @param assay Name of the assay to use.
#' @param .cell_group Optional grouping variable for cells.
#' @param features_number_per_cell_group Number of features to identify per cell group.
#'
#' @return Tibble of variable features by cell type.
#'
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat VariableFeatures
#' @importFrom dplyr mutate
#' @importFrom data.table :=
#' @import dplyr
#' @importFrom glue glue
#' @importFrom purrr map_int
#' @noRd
seurat_to_variable_features_by_cell_type = function(counts, assay, .cell_group = NULL, features_number_per_cell_group = 300){
  
  #Fix GitChecks 
  feature = NULL 
  group = NULL 
  
  .cell_group = enquo(.cell_group)
  
  # Nest
  counts =
    counts |>
    nest(data = -!!.cell_group) |>
    
    # Filter more than 10 cells
    filter(map_int(data, ncol) > 100)
  
  # If I have enough information per cell type
  if(counts |> nrow() |> gt(1)){
    
    # Per cell type
    counts |>
      
      # Get feature within each cluster/cell-type
      mutate(feature = map(
        data,
        ~ .x |>
          FindVariableFeatures(nfeatures = features_number_per_cell_group, assay=assay) |>
          VariableFeatures(assay=assay)
      )) |>
      select(-data) |>
      unnest(feature) |>
      
      # Rename
      rename(group := !!.cell_group)
    
  }
  else {
    
    warning(glue("HPCell says: you have only one distinct `{quo_name(.cell_group)}`, the per-cell-group variable gene detection will be skipped as it would olverlap with the global detection."))
    
    tibble()
  }
}


#' Identify Overall Variable Features in Seurat
#'
#' @description
#' Extracts a set number of overall variable features from a Seurat object.
#'
#' @param counts Seurat object.
#' @param assay Name of the assay to use.
#' @param features_number Number of features to identify.
#'
#' @return Tibble of overall variable features.
#'
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat VariableFeatures
#' @noRd
seurat_to_variable_features_overall = function(counts, assay, features_number = 300){
  
  counts |>
    FindVariableFeatures(nfeatures = features_number, assay=assay) |>
    VariableFeatures(assay=assay) |>
    as_tibble() |>
    rename("feature" = "value") |>
    mutate(group= "variable_overall")
  
}

#' Extract Variable Features from Seurat Object
#'
#' @description
#' Combines the extraction of overall variable features and per cell type variable features.
#'
#' @param counts Seurat object.
#' @param assay Name of the assay to use.
#' @param .sample Sample identifier.
#' @param .cell_group Optional grouping variable for cells.
#' @param features_number_independent_of_cell_groups Number of overall features to identify.
#' @param features_number_per_cell_group Number of features to identify per cell group.
#'
#' @return Combined tibble of overall and per cell type variable features.
#'
#' @importFrom rlang enquo
#' @importFrom rlang is_symbolic
#' @import tidyseurat
#' @importFrom Seurat NormalizeData
#' @importFrom stringr str_subset
#' @importFrom purrr map_int
#' @noRd
seurat_to_variable_features = function(
    counts,
    assay,
    .sample,
    .cell_group = NULL,
    features_number_independent_of_cell_groups = 300,
    features_number_per_cell_group = 300
){
  
  #Fix GitChecks 
  upper_quantile = NULL 
  number_features_overall = NULL 
  
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  
  # If more than one sample balance the size
  if(counts |> distinct(!!.sample) |> nrow() |> gt(1))
    counts =
    counts |>
    
    # Sample up to a plateau to avoid extreme cell_type bias
    nest(data = -!!.sample) %>%
    mutate(n = map_int(data, ~ ncol(.x))) %>%
    mutate(upper_quantile = quantile(n, 0.75) %>% as.integer()) %>%
    mutate(data = map2(
      data, upper_quantile,
      ~ sample_n(.x, min(ncol(.x), .y), replace = FALSE)
    )) %>%
    filter(n>1) %>%
    unnest(data)
  
  # Normalise before - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#normalizing-the-data-1
  counts  = counts |> NormalizeData(assay=assay)
  
  # Drop TCR, MT and RPL, RPS
  counts = counts[rownames(counts) |> str_subset("^MT|^RPL|^RPS|TRAV|TRBV|TRDV|TRGV", negate = TRUE), ]
  
  # Variable overall
  variable_df_overall = seurat_to_variable_features_overall(
    counts,
    assay,
    features_number = features_number_independent_of_cell_groups
  )
  
  # If cell_type_column_for_subsetting == null calculate clusters\
  if(!is_symbolic(.cell_group)){
    
    counts =
      counts |>
      FindVariableFeatures(nfeatures = number_features_overall, assay=assay)  |>
      NormalizeData(assay=assay) |>
      ScaleData(assay=assay) |>
      RunPCA(assay=assay) |>
      FindNeighbors(dims = 1:20) |>
      FindClusters(resolution = 0.5)
    
    .cell_group = as.symbol("seurat_clusters")
  }
  
  variable_df_by_cell_type = seurat_to_variable_features_by_cell_type(
    counts,
    assay,
    .cell_group = !!.cell_group,
    features_number_per_cell_group = features_number_per_cell_group
  )
  
  variable_df_overall  |>
    bind_rows(variable_df_by_cell_type)
  
}

#' Subset Top Rank Variable Genes Across Batches
#'
#' @description
#' Selects top rank variable genes across different batches for a given cell group.
#'
#' @param table_across_cell_groups Tibble with features across cell groups.
#' @param table_within_cell_groups Tibble with features within cell groups.
#' @param .cell_group Cell group variable.
#' @param .batch Batch variable.
#' @param features_number_independent_of_cell_groups Number of overall features to identify.
#' @param features_number_per_cell_group Number of features to identify per cell group.
#'
#' @return Vector of unique top rank variable genes.
#'
#' @importFrom dplyr count
#' @importFrom dplyr with_groups
#' @importFrom dplyr pull
#' @importFrom purrr map_int
#' @noRd
subset_top_rank_variable_genes_across_batches = function(
    table_across_cell_groups,
    table_within_cell_groups,
    .cell_group,
    .batch,
    features_number_independent_of_cell_groups = 2000,
    features_number_per_cell_group = 300
){
  #Fix GitChecks 
  feature = NULL
  
  .cell_group = enquo(.cell_group)
  .batch = enquo(.batch)
  
  batches_with_more_than_10_cell_types =
    table_within_cell_groups |>
    nest(data = -!!.batch) |>
    filter(map_int(data, ~ .x |> distinct(!!.cell_group) |> nrow()) > 10) |>
    unnest(data) |>
    pull(!!.batch)
  
  # Across cell types
  variable_across_cell_types =
    table_across_cell_groups |>
    
    # Filter files that have more than 10 cell types
    # because the genes will be overlapped with the cell type specific
    filter(!!.batch %in% batches_with_more_than_10_cell_types)  |>
    
    count(feature, !!.cell_group) |>
    with_groups(!!.cell_group, ~ .x |> arrange(desc(n)) |> slice(1:features_number_independent_of_cell_groups))
  
  
  # Within cell type
  variable_within_cell_types =
    table_within_cell_groups |>
    count(feature, !!.cell_group) |>
    with_groups(!!.cell_group, ~ .x |> arrange(desc(n)) |> slice(1:features_number_per_cell_group))
  
  # Unique
  bind_rows(
    
    # Cell type specific
    variable_within_cell_types,
    
    # Overall
    variable_across_cell_types
    
  ) |>
    pull(feature) |>
    unique()
}



splitColData <- function(x, f) {
  # This is by @jma1991
  # at https://github.com/drisso/SingleCellExperiment/issues/55
  
  i <- split(seq_along(f), f)
  
  v <- vector(mode = "list", length = length(i))
  
  names(v) <- names(i)
  
  for (n in names(i)) { v[[n]] <- x[, i[[n]]] }
  
  return(v)
  
}

splitRowData <- function(x, f) {
  
  i <- split(seq_along(f), f)
  
  v <- vector(mode = "list", length = length(i))
  
  names(v) <- names(i)
  
  for (n in names(i)) { v[[n]] <- x[i[[n]], ] }
  
  return(v)
  
}
#' Append Code to a Targets Script
#'
#' @description
#' Appends given code to a 'targets' package script.
#'
#' @param code Code to append.
#' @param script Path to the script file.
#'
#' @importFrom readr write_lines
#' @importFrom targets tar_config_get
#' @noRd
tar_script_append = function(code, script = targets::tar_config_get("script")){
  substitute(code) |>
    deparse() |>
    head(-1) |>
    tail(-1) |>
    write_lines(script, append = TRUE)
}

tar_append = function(fx, tiers = NULL, script = targets::tar_config_get("script"), ...){
  
  # Deal with additional argument
  additional_args <- 
    list(...) |> 
    
    # I need this because otherwise the quotation of for example the function names 
    # and the target names will be lost, so those object will be evaluated and 
    # triggered because they do not exist in the environment
    quote_name_classes()
  
  arguments_to_pass  = c(fx)
  
  if(tiers |> is.null() |> not())
    arguments_to_pass = arguments_to_pass |> c(list(tiers = tiers))
  
  if (length(additional_args) > 0)
    arguments_to_pass = arguments_to_pass |> c(additional_args)
  
  # Construct the call with substitute
  # if (length(additional_args) > 0) {
  call_expr = 
    as.call(arguments_to_pass) |> 
    deparse()
  
  # } else {
  #   call_expr <- substitute(fx(x), env = list(fx = fx, x = tiers)) |> 
  #     deparse() 
  # }
  
  # Add prefix
  "target_list |> target_append(" |> 
    c(call_expr ) |> 
    c(")") |> 
    
    paste(collapse = " ") |> 
    
    # Write
    write_lines(script, append = TRUE)
  
}

#' Append Code to a Targets Script
#'
#' @description
#' Appends given code to a 'targets' package script.
#'
#' @param code Code to append.
#' @param script Path to the script file.
#'
#' @importFrom readr write_lines
#' @importFrom targets tar_config_get
#' @noRd
tar_script_append2 = function(code, script = targets::tar_config_get("script"), append = TRUE){
  code |>
    deparse() |>
    head(-1) |>
    tail(-1) |>
    write_lines(script, append = append)
}

#' Append Code to a Targets Script
#'
#' @description
#' Appends given code to a 'targets' package script.
#'
#' @param code Code to append.
#' @param script Path to the script file.
#'
#' @importFrom readr write_lines
#' @importFrom targets tar_config_get
#' @noRd
tar_script_append3 = function(code, script = targets::tar_config_get("script")){
  code |>
    head(-1) |>
    tail(-1) |>
    write_lines(script, append = TRUE)
}

#' Append Code to a Targets Script
#'
#' @description
#' Appends given code to a 'targets' package script.
#'
#' @param code Code to append.
#' @param script Path to the script file.
#'
#' @importFrom readr write_lines
#' @importFrom targets tar_config_get
#' @noRd
append_chunk_fix = function(chunk, script = targets::tar_config_get("script")){
  
  # Add prefix
  "target_list = c(target_list, list(" |> 
    c(
      chunk |>  # cannot start with pipe
        deparse() |> 
        head(-1) |>
        tail(-1) 
    ) |> 
    
    # Add suffix
    c("))") |> 
    
    paste(collapse = " ") |> 
    
    
    # Write
    write_lines(script, append = TRUE)
}

#' Append Code to a Targets Script
#'
#' @description
#' Appends given code to a 'targets' package script.
#'
#' @param code Code to append.
#' @param script Path to the script file.
#'
#' @importFrom readr write_lines
#' @importFrom targets tar_config_get
#' @noRd
append_chunk_tiers = function(chunk, tiers, script = targets::tar_config_get("script")){
  
  # This does not work with purrr:::imap
  # As chunk does not like passed to a function
  tiers = tiers |> get_positions()
  .y = 1
  for(.x in tiers |> names()  ){
    
    
    "target_list = c(target_list, list(" |> 
      c(
        substitute(chunk) |>  # cannot start with pipe
          deparse() |> 
          head(-1) |>
          tail(-1) |> 
          str_replace_all("TIER_PLACEHOLDER", as.character(.y)) |> 
          str_replace_all("SLICE_PLACEHOLDER", tiers[[.y]] |> as.numeric() |> vector_to_code()) |> 
          str_replace("RESOURCE_PLACEHOLDER",  
                      
                      # If not tiering ignore resource naming
                      if_else(
                        length(tiers) == 1,
                        "targets::tar_option_get(\"resources\")",
                        glue("tar_resources(crew = tar_resources_crew(\"{.x}\"))" )
                      )
          )
      ) |> 
      
      # Add suffix
      c("))") |> 
      
      paste(collapse = " ") |> 
      
      # Write
      write_lines(script, append = TRUE)
    
    .y = .y + 1
  }
  
  
}

#' Simple Addition Function
#'
#' @description
#' Performs addition of two numbers.
#'
#' @param a First number.
#' @param b Second number.
#'
#' @return Sum of a and b.
#' @noRd
addition = function(a, b){
  c<-a+b
  c
}

#' Calculate UMAP Embeddings for Seurat Object
#' 
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat RunUMAP
#' @importFrom dplyr as_tibble
#' 
#' @param input_seurat Input data 
#' @return  A tibble containing UMAP coordinates and cluster assignments, to be used for plotting and further analysis 
#'
#' @description
#' Identify variable features, scale the data, run Principal Component Analysis (PCA), find neighbors, identify clusters,
#' and compute UMAP (Uniform Manifold Approximation and Projection) embeddings for 
#' visualization of cell clusters. 
#'
#' @noRd
#' 
#' @importFrom Seurat RunUMAP
#' @export
calc_UMAP <- function(input_seurat) {
  assay_name <- input_seurat@assays |> names() |> extract2(1)
  
  # Check if variable features are already present, if not calculate them
  if (length(VariableFeatures(input_seurat)) == 0) {
    input_seurat <- FindVariableFeatures(input_seurat)
  }
  
  # Extract variable features using VariableFeatures() for Seurat v5
  var_genes <- VariableFeatures(input_seurat)
  
  # Ensure that there are variable features before proceeding
  if (length(var_genes) > 0) {
    # Scale data and run PCA on variable genes
    x <- ScaleData(input_seurat) |>
      RunPCA(features = var_genes) |>
      FindNeighbors(dims = 1:30) |>
      FindClusters(resolution = 0.5) |>
      RunUMAP(dims = 1:30, spread = 0.5, min.dist = 0.01, n.neighbors = 10L) |>
      as_tibble()
  } else {
    stop("No variable features available for UMAP calculation.")
  }
  
  return(x)
}

#' Subsetting input dataset into a list of SingleCellExperiment or Seurat objects by pre-specified sample column tissue 
#' 
#' @importFrom dplyr quo_name pull
#' @importFrom SummarizedExperiment colData
#' 
#' @param sce_obj A `SingleCellExperiment` or `Seurat` object containing input single-cell data
#' @param sample_column The column name specifying sample information
#' 
#' @return The unique sample types in the sample column
#' 
#'  @description
#' Function to subset `SingleCellExperiment` or `Seurat` object by tissue
#' @export
get_unique_tissues <- function(sce_obj, sample_column) {
  sample_column <- quo_name(sample_column)
  
  if (inherits(sce_obj, "Seurat")) {
    unique_sample <-
      sce_obj@meta.data |> pull(sample_column) |> unique()
  } else if (inherits(sce_obj, "SingleCellExperiment")) {
    unique_sample <-
      colData(sce_obj) |> as.data.frame() |> pull(sample_column) |> unique()
  }
  
  unique_sample
}

#' Check for Strong Evidence
#'
#' This function checks for strong evidence in cell annotations.
#'
#' @importFrom rlang enquo
#' @importFrom dplyr case_when
#'
#' @param single_cell_data A data frame containing single-cell data.
#' @param cell_annotation_azimuth_l2 A column representing Azimuth L2 cell annotation.
#' @param cell_annotation_blueprint_singler A column representing Blueprint Singler cell annotation.
#'
#' @return A data frame with a column indicating strong evidence.
#'
# @examples
# single_cell_data <- data.frame(
#   cell_annotation_azimuth_l2 = c("cd14 mono", "b naive"),
#   cell_annotation_blueprint_singler = c("monocytes", "naive b")
#' 
# strong_evidence_data <- is_strong_evidence(single_cell_data, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler)
#'
is_strong_evidence = function(single_cell_data, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler){
  
  cell_annotation_azimuth_l2 = enquo(cell_annotation_azimuth_l2)
  cell_annotation_blueprint_singler = enquo(cell_annotation_blueprint_singler)
  
  single_cell_data |>  
    mutate(strong_evidence = !!cell_annotation_azimuth_l2 == cell_annotation_blueprint_singler) |>
    mutate(strong_evidence = case_when(
      strong_evidence == TRUE ~ strong_evidence,
      !!cell_annotation_azimuth_l2 == "cd14 mono" & cell_annotation_blueprint_singler == "monocytes" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "b naive" & cell_annotation_blueprint_singler == "naive b" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "cd16 mono" & cell_annotation_blueprint_singler == "monocytes" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "plasma" & cell_annotation_blueprint_singler == "plasma" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "b memory" & cell_annotation_blueprint_singler == "class-switched memory b" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "nk_cd56bright" & cell_annotation_blueprint_singler == "nk" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "treg" & cell_annotation_blueprint_singler == "tregs" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "b memory" & cell_annotation_blueprint_singler == "memory b" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "nk proliferating" & cell_annotation_blueprint_singler == "nk" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "cd8 tem" & cell_annotation_blueprint_singler == "cd8 t" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "cdc2" & cell_annotation_blueprint_singler == "dc" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "cd4 naive" & cell_annotation_blueprint_singler == "cd4 t" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "cd4 tcm" & cell_annotation_blueprint_singler == "cd4 t" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "cd8 tcm" & cell_annotation_blueprint_singler == "cd8 t" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "cd4 tem" & cell_annotation_blueprint_singler == "cd4 t" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "sdc" & cell_annotation_blueprint_singler == "dc" ~ TRUE,
      !!cell_annotation_azimuth_l2 == "mait" & cell_annotation_blueprint_singler == "cd8 t" ~ TRUE,
      TRUE ~ FALSE
    ))
}

#' reference_annotation_to_consensus
#'
#' This function takes cell type annotations from multiple datasets (Azimuth, Monaco, Blueprint) and harmonizes them into a consensus annotation. The function utilizes predefined mappings between cell type labels in these datasets to generate standardized cell types across references.
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#' @importFrom dplyr left_join
#' @importFrom dplyr tribble
#' @importFrom tidyr expand_grid
#' @importFrom stringr str_detect
#' @importFrom tibble deframe
#' 
#' @param azimuth_input A vector of cell type annotations from the Azimuth dataset.
#' @param monaco_input A vector of cell type annotations from the Monaco dataset.
#' @param blueprint_input A vector of cell type annotations from the Blueprint dataset.
#'
#' @return A vector of consensus cell type annotations, merging inputs from the three datasets.
#'
#' @examples
#' # Example usage:
#' tibble::tibble(
#'   azimuth_predicted.celltype.l2 = c("CD8 TEM", "NK", "CD4 Naive"),
#'   monaco_first.labels.fine = c("Effector memory CD8 T cells", "Natural killer cells", "Naive CD4 T cells"),
#'   blueprint_first.labels.fine = c("CD8+ Tem", "NK cells", "Naive B-cells")
#' ) |>
#'   dplyr::mutate(consensus = reference_annotation_to_consensus(
#'     azimuth_predicted.celltype.l2, monaco_first.labels.fine, blueprint_first.labels.fine))
#' 
#' @note This function is designed to harmonize specific cell types, especially T cells, B cells, monocytic cells, and innate lymphoid cells (ILCs), across reference datasets.
#'
#' @seealso \code{\link[dplyr]{mutate}}, \code{\link[stringr]{str_detect}}, \code{\link[tidyr]{expand_grid}}
#'
#' @export
reference_annotation_to_consensus = function(azimuth_input, monaco_input, blueprint_input){
  
  # azimuth_pbmc = enquo(azimuth_pbmc)
  # monaco_fine = enquo(monaco_fine)
  # blueprint_fine = enquo(blueprint_fine)
  
  monaco = 
    tribble(
      ~Query,                              ~Reference,
      "Naive CD8 T cells",                 "cd8 naive",
      "Central memory CD8 T cells",        "cd8 tcm",
      "Effector memory CD8 T cells",       "cd8 tem",
      "Terminal effector CD8 T cells",     "cd8 tem",         # Adjusting for the closest match
      "MAIT cells",                        "mait",
      "Vd2 gd T cells",                    "tgd",
      "Non-Vd2 gd T cells",                "tgd",             # No direct match, leaving as NA
      "Follicular helper T cells",         "cd4 fh",
      "T regulatory cells",                "treg",
      "Th1 cells",                         "cd4 th1",
      "Th1/Th17 cells",                    "cd4 th1/th17",
      "Th17 cells",                        "cd4 th17",
      "Th2 cells",                         "cd4 th2",
      "Naive CD4 T cells",                 "cd4 naive",
      "Progenitor cells",                  "progenitor",
      "Naive B cells",                     "b naive",
      "Naive B",                           "b naive",
      "Non-switched memory B cells",       "b memory",        # No direct match, leaving as NA
      "Nonswitched memory B",              "b memory",        # No direct match, leaving as NA
      "Exhausted B cells",                 "plasma",          # Removed " cell"
      "Switched memory B cells",           "b memory",
      "Switched memory B",                 "b memory",
      "Plasmablasts",                      "plasma",          # Removed " cell"
      "Classical monocytes",               "cd14 mono",
      "Intermediate monocytes",            "cd14 mono",       # Mapping to a closely related term
      "Non classical monocytes",           "cd16 mono",
      "Natural killer cells",              "nk",
      "Natural killer",                    "nk",
      "Plasmacytoid dendritic cells",      "pdc",
      "Myeloid dendritic cells",           "cdc",
      "Myeloid dendritic",                 "cdc",
      "Low-density neutrophils",           "granulocyte",
      "Lowdensity neutrophils",            "granulocyte",
      "Low-density basophils",             "granulocyte",     # No direct match, leaving as NA
      "Lowdensity basophils",              "granulocyte",     # No direct match, leaving as NA
      "Terminal effector CD4 T cells",     "cd4 tem",
      "progenitor",                        "progenitor"       # Removed " cell"
    )
  
  azimuth = 
    tribble(
      ~Query,              ~Reference,
      "NK",                "nk",
      "CD8 TEM",           "cd8 tem",
      "CD4 CTL",           "cd4 helper",       # CD4 cytotoxic T lymphocytes often relate to Th1 cells
      "dnT",               "dnt",
      "CD8 Naive",         "cd8 naive",
      "CD4 Naive",         "cd4 naive",
      "CD4 TCM",           "cd4 tcm",          # Central memory cells often relate to Th1 or Th17
      "gdT",               "tgd",
      "CD8 TCM",           "cd8 tcm",
      "MAIT",              "mait",
      "CD4 TEM",           "cd4 tem",          # Effector memory cells can relate to terminal effector cells
      "ILC",               "ilc",
      "CD14 Mono",         "cd14 mono",
      "cDC1",              "cdc",              # Conventional dendritic cell 1 is commonly referred to as CDC
      "pDC",               "pdc",
      "cDC2",              "cdc",              # No specific reference for cDC2, but using CDC as a general category
      "B naive",           "b naive",
      "B intermediate",    "b memory",         # No direct match, leaving as NA
      "B memory",          "b memory",
      "Platelet",          "platelet",
      "Eryth",             "erythrocyte",
      "CD16 Mono",         "cd16 mono",
      "HSPC",              "progenitor",
      "Treg",              "treg",
      "NK_CD56bright",     "nk",               # CD56bright NK cells are a subset of NK cells
      "Plasmablast",       "plasma",
      "NK Proliferating",  "nk",               # NK cells can be proliferative
      "ASDC",              "cdc",              # No direct match, leaving as NA
      "CD8 Proliferating", "cd8 tem",
      "CD4 Proliferating", "cd4 tem",
      "doublet",           "non immune",
      NA,                  NA
    )
  
  blueprint = tribble(
    ~Query,                             ~Reference,
    "Neutrophils",                      "granulocyte",
    "Monocytes",                        "monocyte",
    "MEP",                              "progenitor",      # MEP typically refers to megakaryocyte-erythroid progenitor
    "CD4+ T-cells",                     "cd4 t",
    "Tregs",                            "treg",
    "CD4+ Tcm",                         "cd4 tcm",
    "CD4+ Tem",                         "cd4 tem",
    "CD8+ Tcm",                         "cd8 tcm",
    "CD8+ Tem",                         "cd8 tem",
    "NK cells",                         "nk",
    "naive B-cells",                    "b naive",
    "Memory B-cells",                   "b memory",
    "Class-switched memory B-cells",    "b memory",           # No direct match, leaving as NA
    "HSC",                              "progenitor",      # HSC typically refers to hematopoietic stem cell
    "MPP",                              "progenitor",      # MPP typically refers to multipotent progenitor
    "CLP",                              "progenitor",      # CLP typically refers to common lymphoid progenitor
    "GMP",                              "progenitor",      # GMP typically refers to granulocyte-macrophage progenitor
    "Macrophages",                      "macrophage",
    "CD8+ T-cells",                     "cd8",
    "CD8 T",                            "cd8",
    "Erythrocytes",                     "erythrocyte",
    "Megakaryocytes",                   "megakaryocytes",
    "CMP",                              "progenitor",      # CMP typically refers to common myeloid progenitor
    "Macrophages M1",                   "macrophage m1",      # Specific polarization states (M1, M2) not explicitly listed
    "Macrophages M2",                   "macrophage m2",
    "Endothelial cells",                "endothelial",        # Removed " cell"
    "DC",                               "cdc",                # Assuming DC refers to dendritic cells
    "Eosinophils",                      "granulocyte",        # No direct match, leaving as NA
    "Plasma cells",                     "plasma",             # Removed " cell"
    "Chondrocytes",                     "chondrocyte",
    "Fibroblasts",                      "fibroblast",
    "Smooth muscle",                    "smooth muscle",      # Removed " cell" from "smooth muscle cell"
    "Epithelial cells",                 "epithelial",         # Removed " cell"
    "Melanocytes",                      "melanocyte",
    "Skeletal muscle",                  "muscle",             # "muscle cell" becomes "muscle"
    "Keratinocytes",                    "keratinocyte",
    "mv Endothelial cells",             "endothelial",        # Removed " cell"
    "Myocytes",                         "myocyte",
    "Adipocytes",                       "fat",                # "fat cell" becomes "fat" after removing " cell"
    "Neurons",                          "neuron",
    "Pericytes",                        "pericyte",           # "pericyte cell" becomes "pericyte"
    "Preadipocytes",                    "adipocyte",          # No direct match, leaving as NA
    "Astrocytes",                       "astrocyte",
    "Mesangial cells",                  "mesangial"           # Removed " cell"
  )
  
  non_immune_cells <- c(
    "megakaryocytes",
    "endothelial",
    "chondrocyte",
    "fibroblast",
    "smooth muscle",
    "epithelial",
    "melanocyte",
    "muscle",
    "keratinocyte",
    "endothelial",  # Appears again in the original vector
    "myocyte",
    "fat",
    "neuron",
    "pericyte",
    "adipocyte",
    "astrocyte",
    "mesangial"
  )
  
  t_cells <- c(
    "cd8 naive",
    "cd8 tcm",
    "cd8 tem",
    "cd4 tem",
    "cd4 tcm",
    "cd4 effector",
    "treg",
    "cd4 th1/th17",
    "cd4 th1",
    "cd4 th17",
    "cd4 th2",
    "cd4 t",
    "t_nk",
    "cd8 effector",
    "dnt",
    "cd4 naive",
    "cd4 th2",
    "tgd",
    "cd4 fh",
    "mait"
  )
  
  b_cells <- c(
    "b naive",
    "b memory",
    "plasma"
  )
  
  myeloid_cells <- c(
    "cd14 mono",
    "monocyte",
    "cd16 mono",
    "macrophage",
    "macrophage m1",
    "macrophage m2",
    "macrophages",
    #"pdc",  # Plasmacytoid dendritic cells
    "cdc",  # Conventional dendritic cells
    "kupffer"
  )
  
  ilcs <- c(
    "ilc",
    "nk"
  )
  
  
  all_combinations = 
    expand_grid(
      blueprint_fine = blueprint |> pull(Reference) |> unique(), 
      monaco_fine = monaco |> pull(Reference) |> unique(),
      azimuth_pbmc = azimuth |> pull(Reference) |> unique()
    ) |> 
    
    # Find consensus manually
    mutate(consensus =
             case_when(
               
               # Non immune
               blueprint_fine %in% non_immune_cells ~ "non immune",
               
               # Full consensus
               blueprint_fine  ==  monaco_fine  &
                 blueprint_fine  ==  azimuth_pbmc  ~  blueprint_fine ,
               
               # This goes before partial exact consensus because is a special case
               monaco_fine %in% c("cd4 fh","cd4 th1","cd4 th1/th17", "cd4 th17", "cd4 th2") & blueprint_fine %in% c("cd4 tcm") & azimuth_pbmc %in% c("cd4 tcm") ~ glue("{monaco_fine} cm") , # Because most Th cells are central and effector memory CD4 T cells (CM and EM), PMID: 30726743
               monaco_fine %in% c("cd4 fh","cd4 th1","cd4 th1/th17", "cd4 th17", "cd4 th2") & blueprint_fine %in% c("cd4 tem") & azimuth_pbmc %in% c("cd4 tem") ~ glue("{monaco_fine} em") , # Because most Th cells are central and effector memory CD4 T cells (CM and EM), PMID: 30726743
               blueprint_fine |> str_detect("macrophage") & monaco_fine |> str_detect(" mono") & azimuth_pbmc |> str_detect(" mono") ~ blueprint_fine,
               
               # Partial consensus
               blueprint_fine  ==  monaco_fine  ~  blueprint_fine ,
               blueprint_fine  ==  azimuth_pbmc  ~  blueprint_fine ,
               monaco_fine  ==  azimuth_pbmc  ~  monaco_fine ,
               
               ##################
               # More difficoult combination if nothing above matched
               ##################
               
               # T cells
               str_detect( blueprint_fine , "cd8") & str_detect( monaco_fine , "cd8") & str_detect( azimuth_pbmc , "cd8") ~ "t cd8",
               
               
               str_detect( blueprint_fine , "cd4|treg") & str_detect( monaco_fine , "cd4|treg") & str_detect( azimuth_pbmc , "cd4|treg") ~ "t cd4",
               blueprint_fine  %in% t_cells &  monaco_fine  %in% t_cells &  azimuth_pbmc  %in% t_cells ~ "t",
               
               # B cells
               blueprint_fine  %in% b_cells &  monaco_fine  %in% b_cells &  azimuth_pbmc  %in% b_cells ~ "b",
               
               # Monocytic cells
               blueprint_fine  %in% myeloid_cells &  monaco_fine  %in% myeloid_cells &  azimuth_pbmc  %in% myeloid_cells ~ "monocytic",
               
               # ILCs
               blueprint_fine  %in% ilcs &  monaco_fine  %in% ilcs &  azimuth_pbmc  %in% ilcs ~ "ilc",
               
               # cytotoxic
               (  blueprint_fine  %in% ilcs | str_detect( blueprint_fine , "cd8") ) & 
                 (  monaco_fine  %in% ilcs | str_detect( monaco_fine , "cd8") ) &
                 (  azimuth_pbmc  %in% ilcs | str_detect( azimuth_pbmc , "cd8") ) ~ "cytotoxic",
               
               ##################
               # Partial consensus broad cell types
               ##################
               
               # T cells
               str_detect( blueprint_fine , "cd8") & str_detect( monaco_fine , "cd8") ~ "t cd8",
               str_detect( blueprint_fine , "cd8") & str_detect( azimuth_pbmc , "cd8") ~ "t cd8",
               str_detect( monaco_fine , "cd8") & str_detect( azimuth_pbmc , "cd8") ~ "t cd8",
               
               monaco_fine %in% c("cd4 fh","cd4 th1","cd4 th1/th17", "cd4 th17", "cd4 th2") & blueprint_fine %in% c("cd4 tcm") ~ glue("{monaco_fine} cm") , # Because most Th cells are central and effector memory CD4 T cells (CM and EM), PMID: 30726743
               monaco_fine %in% c("cd4 fh","cd4 th1","cd4 th1/th17", "cd4 th17", "cd4 th2") & azimuth_pbmc %in% c("cd4 tcm") ~ glue("{monaco_fine} cm") , # Because most Th cells are central and effector memory CD4 T cells (CM and EM), PMID: 30726743
               
               
               monaco_fine %in% c("cd4 fh","cd4 th1","cd4 th1/th17", "cd4 th17", "cd4 th2") & blueprint_fine %in% c("cd4 tem") ~ glue("{monaco_fine} em") , # Because most Th cells are central and effector memory CD4 T cells (CM and EM), PMID: 30726743
               monaco_fine %in% c("cd4 fh","cd4 th1","cd4 th1/th17", "cd4 th17", "cd4 th2") & azimuth_pbmc %in% c("cd4 tem") ~ glue("{monaco_fine} em") , # Because most Th cells are central and effector memory CD4 T cells (CM and EM), PMID: 30726743
               
               
               str_detect( blueprint_fine , "cd4|treg") & str_detect( monaco_fine , "cd4|treg") ~ "t cd4",
               str_detect( blueprint_fine , "cd4|treg") & str_detect( azimuth_pbmc , "cd4|treg") ~ "t cd4",
               str_detect( monaco_fine , "cd4|treg") & str_detect( azimuth_pbmc , "cd4|treg") ~ "t cd4",
               
               blueprint_fine  %in% t_cells &  monaco_fine  %in% t_cells  ~ "t",
               blueprint_fine  %in% t_cells  &  azimuth_pbmc  %in% t_cells ~ "t",
               monaco_fine  %in% t_cells &  azimuth_pbmc  %in% t_cells ~ "t",
               
               # B cells
               blueprint_fine  %in% b_cells &  monaco_fine  %in% b_cells  ~ "b",
               blueprint_fine  %in% b_cells  &  azimuth_pbmc  %in% b_cells ~ "b",
               monaco_fine  %in% b_cells &  azimuth_pbmc  %in% b_cells ~ "b",
               
               # Monocytic cells
               blueprint_fine |> str_detect("monocyte") & monaco_fine |> str_detect(" mono") ~ monaco_fine,  # This is because blueprint does not have CDC16 or CD14
               blueprint_fine |> str_detect("monocyte") & azimuth_pbmc |> str_detect(" mono") ~ azimuth_pbmc,  # This is because blueprint does not have CDC16 or CD14
               
               blueprint_fine |> str_detect("macrophage") & monaco_fine |> str_detect(" mono") ~ blueprint_fine,  # This is because only blueprint has mac M1 M2
               blueprint_fine |> str_detect("macrophage") & azimuth_pbmc |> str_detect(" mono") ~ blueprint_fine,  # This is because only blueprint has mac M1 M2
               
               blueprint_fine  %in% myeloid_cells &  monaco_fine  %in% myeloid_cells  ~ "monocytic",
               blueprint_fine  %in% myeloid_cells  &  azimuth_pbmc  %in% myeloid_cells ~ "monocytic",
               monaco_fine  %in% myeloid_cells &  azimuth_pbmc  %in% myeloid_cells ~ "monocytic",
               
               # ILCs
               blueprint_fine  %in% ilcs &  monaco_fine  %in% ilcs  ~ "ilc",
               blueprint_fine  %in% ilcs  &  azimuth_pbmc  %in% ilcs ~ "ilc",
               monaco_fine  %in% ilcs &  azimuth_pbmc  %in% ilcs ~ "ilc",
               
               # cytotoxic
               (  blueprint_fine  %in% ilcs | str_detect( blueprint_fine , "cd8") ) & 
                 (  monaco_fine  %in% ilcs | str_detect( monaco_fine , "cd8") )  ~ "cytotoxic",
               
               (  blueprint_fine  %in% ilcs | str_detect( blueprint_fine , "cd8") ) & 
                 (  azimuth_pbmc  %in% ilcs | str_detect( azimuth_pbmc , "cd8") ) ~ "cytotoxic",
               
               (  monaco_fine  %in% ilcs | str_detect( monaco_fine , "cd8") ) &
                 (  azimuth_pbmc  %in% ilcs | str_detect( azimuth_pbmc , "cd8") ) ~ "cytotoxic",
               
               
               TRUE ~ NA_character_
             )) |> 
    
    # simplify Thelper cm to tcm
    mutate(consensus = if_else(consensus |> str_detect("cd4 .* cm"), "cd4 tcm", consensus  ))
  
  # |> 
  #   rowid_to_column("combination_id") |> 
  #   pivot_longer(-combination_id, names_to = "Database", values_to = "Reference") |> 
  #   left_join(conversion_table, relationship = "many-to-many") |> 
  #   select(-Reference) |> 
  #   pivot_wider(names_from = Database, values_from = Query, values_fn = function(x) paste(unique(x), collapse = ","))
  
  # parse names, chenge to lower case for all
  tibble(
    blueprint_fine = blueprint |> mutate(across(everything(), tolower)) |> deframe() |> _[!!tolower(blueprint_input)],
    monaco_fine = monaco |> mutate(across(everything(), tolower)) |> deframe() |> _[!!tolower(monaco_input)],
    azimuth_pbmc = azimuth |> mutate(across(everything(), tolower)) |> deframe() |> _[!!tolower(azimuth_input)],
  ) |> 
    left_join(
      all_combinations,
      by = join_by(
        blueprint_fine == blueprint_fine,
        monaco_fine == monaco_fine,
        azimuth_pbmc == azimuth_pbmc
      )
    ) |> 
    pull(consensus)
  
  
}

#' Clean and Standardize Cell Type Names
#'
#' Cleans and standardizes a vector of cell type names by applying a series of string transformations to improve consistency.
#' This function is particularly useful for preprocessing cell type labels in biological datasets where consistent naming conventions are important.
#'
#' @param x A character vector of cell type names to be cleaned and standardized.
#'
#' @return A character vector of cleaned and standardized cell type names.
#'
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_remove
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_trim
#'
#' @examples
#' cell_types <- c("CD4+ T-cells", "NK cells", "Blast-cells", "Terminally differentiated macrophage")
#' cleaned_cell_types <- clean_cellxgene_cell_types(cell_types)
#' print(cleaned_cell_types)
#'
#' # Output:
#' # [1] "cd4 t"     "nk"        ""          "macrophage"
#' 
#' @export
clean_cellxgene_cell_types = function(x){
  
  x |>
    # Annotate
    tolower() |>
    str_remove_all(",") |>
    str_remove("alphabeta") |>
    str_remove_all("positive") |>
    str_replace("cd4  t", "cd4") |>
    str_replace("regulatory t", "treg") |>
    str_remove("thymusderived") |>
    str_remove("human") |>
    str_remove("igg ") |>
    str_remove("igm ") |>
    str_remove("iga ") |>
    str_remove("group [0-9]") |>
    str_remove("common") |>
    str_remove("cd45ro") |>
    str_remove("type i") |>
    str_remove("germinal center") |>
    str_remove("iggnegative") |>
    str_remove("terminally differentiated") |>
    
    str_replace(".*macrophage.*", "macrophage") |>
    str_replace("^mononuclear phagocyte$", "macrophage") |>
    str_replace(".* treg.*", "treg") |>
    str_replace(".* dendritic.*", "dendritic") |>
    str_replace(".* thelper.*", "thelper") |>
    str_replace(".*thelper .*", "thelper") |>
    str_replace(".*gammadelta.*", "tgd") |>
    str_replace(".*natural killer.*", "nk") |> 
    
    str_replace_all("  ", " ") |>
    
    str_replace("myeloid leukocyte", "myeloid") |>
    str_replace("effector memory", "tem") |>
    str_replace("effector", "tem") |>
    str_replace_all("cd8 t", "cd8") |>
    str_replace("central memory", "tcm") |>
    str_replace("gammadelta t", "gdt") |>
    str_replace("nonclassical monocyte", "cd16 monocyte") |>
    str_replace("classical monocyte", "cd14 monocyte") |>
    str_replace("follicular b", "b") |>
    str_replace("unswitched memory", "memory") |>
    
    str_trim() |>
    
    str_remove_all("\\+") |>
    str_remove_all("cells") |>
    str_remove_all("cell") |>
    str_remove_all("blast") |>
    str_remove_all("-") |>
    str_trim() |> 
    
    str_remove("^_+|_+$") |>  # Removes leading and trailing underscores
    
    # clean NON IMMUNE
    str_replace("(?i)\\bepithelial\\b", "epithelial_cell") |>
    str_replace("(?i)\\bfibroblast\\b", "fibroblast") |>
    str_replace("(?i)\\bendothelial\\b", "endothelial_cell") |>
    str_replace("(?i)^(Mueller cell|Muller cell)$", "Muller_cell") |>
    str_replace("(?i)\\bneuron\\b", "neuron") |>
    str_replace("(?i)amplifying cell", "amplifying_cell") |>
    str_replace("(?i)stem cell", "stem_cell") |>
    str_replace("(?i)progenitor cell", "progenitor_cell") |>
    str_replace("(?i)acinar cell", "acinar_cell") |>
    str_replace("(?i)goblet cell", "goblet_cell") |>
    str_replace("(?i)thymocyte", "thymocyte") |>
    str_replace("(?i)urothelial", "urothelial_cell") |>
    str_replace("(?i)\\bfat\\b", "fat_cell") |>
    str_replace("(?i)pneumocyte", "pneumocyte") |>
    str_replace("(?i)mesothelial", "mesothelial_cell") |>
    str_replace("(?i)enteroendocrine", "enteroendocrine_cell") |>
    str_replace("(?i)enterocyte", "enterocyte") |>
    str_replace("(?i)\\bbasal\\b", "basal_cell") |>
    str_replace("(?i)stromal", "stromal_cell") |>
    str_replace("(?i)retina", "retinal_cell") |>
    str_replace("(?i)ciliated", "ciliated_cell") |>
    str_replace("(?i)pericyte", "pericyte_cell") |>
    str_replace("(?i)trophoblast", "trophoblast") |>
    str_replace("(?i)brush", "brush_cell") |>
    str_replace("(?i)serous", "serous_cell") |>
    str_replace("(?i)hepatocyte", "hepatocyte") |>
    str_replace("(?i)melanocyte", "melanocyte") |>
    str_replace("(?i)myocyte", "myocyte") |>
    str_replace("(?i)promyelocyte", "promyelocyte") |>
    str_replace("(?i)cholangiocyte", "cholangiocyte") |>
    str_replace("(?i)myoblast", "myoblast") |>
    str_replace("(?i)satellite", "satellite_cell") |>
    str_replace("(?i)muscle", "muscle_cell") |>
    str_replace("(?i)progenitor", "progenitor_cell") |>
    str_replace("(?i)erythrocyte", "erythrocyte") |>
    str_replace("(?i)myoepithelial", "myoepithelial_cell") |>
    str_replace("(?i)myofibroblast", "myofibroblast_cell") |>
    str_replace("(?i)pancreatic", "pancreatic_cell") |>
    str_replace("(?i)renal", "renal_cell") |>
    str_replace("(?i)epidermal", "epidermal_cell") |>
    str_replace("(?i)cortical", "cortical_cell") |>
    str_replace("(?i)interstitial", "interstitial_cell") |>
    str_replace("(?i)neuroendocrine", "neuroendocrine_cell") |>
    str_replace("(?i)granular", "granular_cell") |>
    str_replace("(?i)kidney", "kidney_cell") |>
    str_replace("(?i)paneth", "paneth_cell") |>
    str_replace("(?i)bipolar", "bipolar_cell") |>
    str_replace_all(" ", "_")
}


#' Harmonize Non-Immune Cell Type Names
#'
#' This function harmonizes non-immune cell type names in the metadata.
#'
#'
#' @param metadata A data frame containing cell type information.
#'
#' @return The metadata with harmonized cell type names.
#'
#' @examples
#' metadata <- data.frame(cell_type = c("Myofibroblast", "Fibroblast", "Other Fibroblast"))
# harmonized_metadata <- harmonise_names_non_immune(metadata)
#'
harmonise_names_non_immune = function(metadata){
  
  # grep("\\bfibroblast\\b", c("myofibroblast", "fibroblast", "other fibroblast"))
  
  metadata$cell_type_harmonised =  metadata$cell_type
  
  metadata$cell_type_harmonised <- ifelse(grepl("\\bepithelial\\b", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "epithelial_cell",  ## does not include myoepithelial.
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("\\bfibroblast\\b", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "fibroblast",  ## does not include myofibroblast.
                                          metadata$cell_type_harmonised)
  
  
  metadata$cell_type_harmonised <- ifelse(grepl("\\bendothelial\\b", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "endothelial_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(metadata$cell_type_harmonised=="Mueller cell" | metadata$cell_type_harmonised=="Muller cell",
                                          "Muller_cell",  ## Merge Mueller cells and Muller cells.
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("\\bneuron\\b", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "neuron",  ## does not include interneuron and neuronal receptor cell.
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("amplifying cell", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "amplifying_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("stem cell", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "stem_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("progenitor cell", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "progenitor_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("acinar cell", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "acinar_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("goblet cell", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "goblet_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("thymocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "thymocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("urothelial", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "urothelial_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("fat", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "fat_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("pneumocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "pneumocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("pneumocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "pneumocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("mesothelial", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "mesothelial_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("enteroendocrine", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "enteroendocrine_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("enterocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "enterocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("basal", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "basal_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("stromal", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "stromal_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("retina", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "retinal_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("ciliated", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "ciliated_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("pericyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "pericyte_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("trophoblast", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "trophoblast",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("brush", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "brush_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("serous", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "serous cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("hepatocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "hepatocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("melanocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "melanocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("myocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "myocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("promyelocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "promyelocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("cholangiocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "cholangiocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("melanocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "melanocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("myoblast", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "myoblast", ## Discussed with Stefano on Teams on 16/12/2022.
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("satellite", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "satellite_cell", ## Discussed with Stefano on Teams on 16/12/2022.
                                          metadata$cell_type_harmonised)
  
  
  metadata$cell_type_harmonised <- ifelse(grepl("muscle", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "muscle_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("progenitor", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "progenitor_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("progenitor", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "progenitor_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("erythrocyte", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "erythrocyte",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("myoepithelial", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "myoepithelial_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("myofibroblast", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "myofibroblast_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("pancreatic", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "pancreatic_cell",  ## Discussed with Stefano on Teams on 19/12/2022.
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("renal", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "renal_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("epidermal", metadata$cell_type_harmonised, ignore.case=TRUE),  ## Includes epidermal cell and epidermal Langerhans cell.
                                          "epidermal_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("cortical", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "cortical_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("interstitial", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "interstitial_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("neuroendocrine", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "neuroendocrine_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("granular", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "granular_cell", ## Discussed with Stefano on Teams on 19/12/2022.
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("kidney", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "kidney_cell",  ## Discussed with Stefano on Teams on 19/12/2022.
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("paneth", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "paneth_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- ifelse(grepl("bipolar", metadata$cell_type_harmonised, ignore.case=TRUE),
                                          "bipolar_cell",
                                          metadata$cell_type_harmonised)
  
  metadata$cell_type_harmonised <- gsub(" " , "_", metadata$cell_type_harmonised)
  
  
  
  table(metadata$cell_type_harmonised[grepl("glial", metadata$cell_type_harmonised, ignore.case=TRUE)])
  ##  glial cell, microglial cell, radial glial cell
  ## https://www.simplypsychology.org/glial-cells.html#:~:text=Glial%20cells%20are%20a%20general,that%20keep%20the%20brain%20functioning.
  
  table(metadata$cell_type_harmonised[grepl("hematopoietic", metadata$cell_type_harmonised, ignore.case=TRUE)])
  ## hematopoietic cell, hematopoietic precursor cell
  ## precursor cell = stem_cell?
  
  table(metadata$cell_type_harmonised[grepl("papillary", metadata$cell_type_harmonised, ignore.case=TRUE)])
  ## papillary tips cell = renal_cell ?
  
  
  metadata
}

# get_manually_curated_immune_cell_types = function(){
#   
#   # library(zellkonverter)
#   # library(Seurat)
#   # library(SingleCellExperiment) # load early to avoid masking dplyr::count()
#   # library(tidySingleCellExperiment)
#   # library(dplyr)
#   # library(cellxgenedp)
#   # library(tidyverse)
#   #library(tidySingleCellExperiment)
#   # library(stringr)
#   # library(scMerge)
#   # library(glue)
#   # library(tidyseurat)
#   # library(celldex)
#   # library(SingleR)
#   # library(glmGamPoi)
#   # library(stringr)
#   # library(purrr)
#   
#   
#   #Fix GCHECKS 
#   metadata_file = NULL 
#   .cell = NULL 
#   cell_type = NULL
#   file_id = NULL 
#   .sample = NULL 
#   azhimut_confirmed = NULL 
#   blueprint_confirmed <- NULL
#   arrange <- NULL # This one is actually a function from dplyr, so you should use it with dplyr::arrange or import it
#   cell_type_clean <- NULL
#   blueprint_singler <- NULL
#   predicted.celltype.l2 <- NULL
#   strong_evidence <- NULL
#   cell_type_harmonised <- NULL
#   confidence_class <- NULL
#   lineage_1 <- NULL
#   monaco_singler <- NULL
#   cell_annotation_monaco_singler <- NULL
#   cell_annotation_azimuth_l2 <- NULL
#   cell_annotation_blueprint_singler <- NULL
#   confidence_class_manually_curated <- NULL
#   cell_type_harmonised_manually_curated <- NULL
#   file_curated_annotation_merged <- NULL
#   .sample <- NULL
#   cell_type_harmonised_non_immune <- NULL
# 
#   # library(zellkonverter)
#   # library(Seurat)
#   # library(SingleCellExperiment) # load early to avoid masking dplyr::count()
#   # library(tidySingleCellExperiment)
#   # library(dplyr)
#   # library(cellxgenedp)
#   # library(tidyverse)
#   # #library(tidySingleCellExperiment)
#   # library(stringr)
#   # library(scMerge)
#   # library(glue)
#   # library(DelayedArray)
#   # library(HDF5Array)
#   # library(tidyseurat)
#   # library(celldex)
#   # library(SingleR)
#   # library(glmGamPoi)
#   # library(stringr)
#   # library(purrr)
#   
#   # # source("utility.R")
#   # 
#   # metadata_file = "/vast/projects/cellxgene_curated//metadata_0.2.rds"
#   # file_curated_annotation_merged = "~/PostDoc/CuratedAtlasQueryR/dev/cell_type_curated_annotation_0.2.3.rds"
#   # file_metadata_annotated = "/vast/projects/cellxgene_curated/metadata_annotated_0.2.3.rds"
#   # annotation_directory = "/vast/projects/cellxgene_curated//annotated_data_0.2/"
#   # 
#   # # metadata_file = "/vast/projects/cellxgene_curated//metadata.rds"
#   # # file_curated_annotation_merged = "~/PostDoc/CuratedAtlasQueryR/dev/cell_type_curated_annotation.rds"
#   # # file_metadata_annotated = "/vast/projects/cellxgene_curated//metadata_annotated.rds"
#   # # annotation_directory = "/vast/projects/cellxgene_curated//annotated_data_0.1/"
#   # 
#   # 
#   # annotation_harmonised =
#   #   dir(annotation_directory, full.names = TRUE) |>
#   #   enframe(value="file") |>
#   #   tidyr::extract(  file,".sample", "/([a-z0-9]+)\\.rds", remove = F) |>
#   #   mutate(data = map(file, ~ .x |> readRDS() |> select(-contains("score")) )) |>
#   #   unnest(data) |>
#   #   
#   #   # Format
#   #   mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler, monaco_singler, ),	tolower	)) |>
#   #   mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler, monaco_singler, ),	clean_cell_types	)) |>
#   #   
#   #   # Format
#   # is_strong_evidence(predicted.celltype.l2, blueprint_singler) |> 
#   #   
#   # 
#   # 
#   # 
#   # job::job({
#   #   annotation_harmonised |>  saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/annotated_data_0.2_temp_table.rds")
#   # })
#   # 
# 
#   annotation_harmonised = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/annotated_data_0.2_temp_table.rds")
#   
#   # library(CuratedAtlasQueryR)
#   metadata_df = readRDS(metadata_file)
#   
#   # Integrate with metadata
#   
#   annotation =
#     metadata_df |>
#     select(.cell, cell_type, file_id, .sample) |>
#     as_tibble() |>
#     left_join(read_csv("~/PostDoc/CuratedAtlasQueryR/dev/metadata_cell_type.csv"),  by = "cell_type") |>
#     left_join(annotation_harmonised, by = c(".cell", ".sample")) |>
#     
#     # Clen cell types
#     mutate(cell_type_clean = cell_type |> clean_cell_types())
#   
#   # annotation |>
#   # 	filter(lineage_1=="immune") |>
#   # 	count(cell_type, predicted.celltype.l2, blueprint_singler, strong_evidence) |>
#   # 	arrange(!strong_evidence, desc(n)) |>
#   # 	write_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm.csv")
#   
#   
#   annotation_crated_confirmed =
#     read_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm_manually_curated.csv") |>
#     
#     # TEMPORARY
#     rename(cell_type_clean = cell_type) |>
#     
#     filter(!is.na(azhimut_confirmed) | !is.na(blueprint_confirmed)) |>
#     filter(azhimut_confirmed + blueprint_confirmed > 0) |>
#     
#     # Format
#     mutate(cell_type_harmonised = case_when(
#       azhimut_confirmed ~ predicted.celltype.l2,
#       blueprint_confirmed ~ blueprint_singler
#     )) |>
#     
#     mutate(confidence_class = 1)
#   
#   
#   
#   # To avoid immune cell annotation if very contrasting evidence
#   blueprint_definitely_non_immune = c(   "astrocytes" , "chondrocytes"  , "endothelial"  ,  "epithelial" ,  "fibros"  ,  "keratinocytes" ,    "melanocytes"  , "mesangial"  ,  "mv endothelial",   "myocytes" ,  "neurons"  ,  "pericytes" ,  "preadipocytes" , "skeletal muscle"  ,  "smooth muscle"      )
#   
#   
#   
#   annotation_crated_UNconfirmed =
#     
#     # Read
#     read_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm_manually_curated.csv") |>
#     
#     # TEMPORARY
#     rename(cell_type_clean = cell_type) |>
#     
#     filter(is.na(azhimut_confirmed) | (azhimut_confirmed + blueprint_confirmed) == 0) |>
#     
#     clean_cell_types_deeper() |> 
#     
#     mutate(cell_type_harmonised = "") |>
#     
#     # Classify strong evidence
#     mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("cd8 cytokine secreting tem t") & blueprint_singler == "nk", T, blueprint_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("cd8 cytotoxic t") & blueprint_singler == "nk",  T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("cd8alphaalpha intraepithelial t") & predicted.celltype.l2 == "cd8 tem" & blueprint_singler == "cd8 tem", T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("mature t") & strong_evidence & predicted.celltype.l2  |> str_detect("tem|tcm"), T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("myeloid") & strong_evidence & predicted.celltype.l2  == "cd16 mono", T, azhimut_confirmed) ) |>
#     
#     # Classify weak evidence
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("b", "B") & predicted.celltype.l2   == "b memory" & blueprint_singler == "classswitched memory b", T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("b", "B") & predicted.celltype.l2   %in% c("b memory", "b intermediate", "b naive", "plasma") & !blueprint_singler %in% c("classswitched memory b", "memory b", "naive b"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("b", "B") & !predicted.celltype.l2   %in% c("b memory", "b intermediate", "b naive") & blueprint_singler %in% c("classswitched memory b", "memory b", "naive b", "plasma"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "activated cd4" & predicted.celltype.l2  %in% c("cd4 tcm", "cd4 tem", "tregs"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "activated cd4" & blueprint_singler  %in% c("cd4 tcm", "cd4 tem", "tregs"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "activated cd8" & predicted.celltype.l2  %in% c("cd8 tcm", "cd8 tem"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "activated cd8" & blueprint_singler  %in% c("cd8 tcm", "cd8 tem"), T, blueprint_confirmed) ) |>
#     
#     # Monocyte macrophage
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 cd16 monocyte" & predicted.celltype.l2  %in% c("cd14 mono", "cd16 mono"), T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 cd16negative classical monocyte" & predicted.celltype.l2  %in% c("cd14 mono"), T, azhimut_confirmed) ) |>
#     mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14 cd16negative classical monocyte" & blueprint_singler  %in% c("monocytes"), "cd14 mono", cell_type_harmonised) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 monocyte" & predicted.celltype.l2  %in% c("cd14 mono"), T, azhimut_confirmed) ) |>
#     mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14 monocyte" & blueprint_singler  %in% c("monocytes"), "cd14 mono", cell_type_harmonised) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14low cd16 monocyte" & predicted.celltype.l2  %in% c("cd16 mono"), T, azhimut_confirmed) ) |>
#     mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14low cd16 monocyte" & blueprint_singler  %in% c("monocytes"), "cd16 mono", cell_type_harmonised) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd16 monocyte" & predicted.celltype.l2  %in% c("cd16 mono"), T, azhimut_confirmed) ) |>
#     mutate(cell_type_harmonised = if_else(cell_type_clean == "cd16 monocyte" & blueprint_singler  %in% c("monocytes"), "cd16 mono", cell_type_harmonised) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "monocyte" & blueprint_singler  |> str_detect("monocyte|macrophage") & !predicted.celltype.l2 |> str_detect(" mono"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "monocyte" & predicted.celltype.l2 |> str_detect(" mono"), T, azhimut_confirmed) ) |>
#     
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd4" & predicted.celltype.l2 |> str_detect("cd4|treg") & !blueprint_singler  |> str_detect("cd4"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "cd4" & !predicted.celltype.l2 |> str_detect("cd4") & blueprint_singler  |> str_detect("cd4|treg"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>
#     
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "memory t" & predicted.celltype.l2 |> str_detect("tem|tcm") & !blueprint_singler  |> str_detect("tem|tcm"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "memory t" & !predicted.celltype.l2 |> str_detect("tem|tcm") & blueprint_singler  |> str_detect("tem|tcm"), T, blueprint_confirmed) ) |>
#     
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8alphaalpha intraepithelial t" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8alphaalpha intraepithelial t" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8hymocyte" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8hymocyte" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>
#     
#     # B cells
#     mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("memory b") & predicted.celltype.l2 =="b memory", T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("memory b") & blueprint_singler |> str_detect("memory b"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "immature b" & predicted.celltype.l2 =="b naive", T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "immature b" & blueprint_singler |> str_detect("naive b"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "mature b" & predicted.celltype.l2 %in% c("b memory", "b intermediate"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "mature b" & blueprint_singler |> str_detect("memory b"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "naive b" & predicted.celltype.l2 %in% c("b naive"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "naive b" & blueprint_singler |> str_detect("naive b"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "transitional stage b" & predicted.celltype.l2 %in% c("b intermediate"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "transitional stage b" & blueprint_singler |> str_detect("naive b") & !predicted.celltype.l2 %in% c("b intermediate"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "memory b" & predicted.celltype.l2 %in% c("b intermediate"), T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & predicted.celltype.l2 %in% c("b naive") & !blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & blueprint_singler |> str_detect("naive b") & predicted.celltype.l2 %in% c("hspc"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & predicted.celltype.l2 %in% c("hspc"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
#     
#     # Plasma cells
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "plasma") & predicted.celltype.l2 == "plasma" , T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "plasma") & predicted.celltype.l2 == "plasma" , T, blueprint_confirmed) ) |>
#     
#     mutate(azhimut_confirmed = case_when(
#       cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & predicted.celltype.l2 == "cd4 ctl" & blueprint_singler != "cd4 tcm" ~ T,
#       cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & predicted.celltype.l2 == "cd4 tem" & blueprint_singler != "cd4 tcm" ~ T,
#       TRUE ~ azhimut_confirmed
#     ) ) |>
#     mutate(blueprint_confirmed = case_when(
#       cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & blueprint_singler == "cd4 tem" & predicted.celltype.l2 != "cd4 tcm" ~ T,
#       cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & blueprint_singler == "cd4 t" & predicted.celltype.l2 != "cd4 tcm" ~ T,
#       TRUE ~ blueprint_confirmed
#     ) ) |>
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "cd4hymocyte" & predicted.celltype.l2 |> str_detect("cd4|treg") & !blueprint_singler  |> str_detect("cd4"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "cd4hymocyte" & !predicted.celltype.l2 |> str_detect("cd4") & blueprint_singler  |> str_detect("cd4|treg"), T, blueprint_confirmed) ) |>
#     
#     mutate(azhimut_confirmed = case_when(
#       cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 == "cd8 tem" & blueprint_singler != "cd8 tcm" ~ T,
#       cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 == "cd8 tcm" & blueprint_singler != "cd8 tem" ~ T,
#       TRUE ~ azhimut_confirmed
#     ) ) |>
#     mutate(blueprint_confirmed = case_when(
#       cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 != "cd8 tem" & blueprint_singler == "cd8 tcm" ~ T,
#       cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 != "cd8 tcm" & blueprint_singler == "cd8 tem" ~ T,
#       TRUE ~ blueprint_confirmed
#     ) ) |>
#     
#     mutate(azhimut_confirmed = case_when(
#       cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 == "cd4 tem" & blueprint_singler != "cd8 tcm" ~ T,
#       cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 == "cd4 tcm" & blueprint_singler != "cd8 tem" ~ T,
#       TRUE ~ azhimut_confirmed
#     ) ) |>
#     mutate(blueprint_confirmed = case_when(
#       cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 != "cd4 tem" & blueprint_singler == "cd4 tcm" ~ T,
#       cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 != "cd4 tcm" & blueprint_singler == "cd4 tem" ~ T,
#       TRUE ~ blueprint_confirmed
#     ) ) |>
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "t") & blueprint_singler =="cd8 t" & predicted.celltype.l2 |> str_detect("cd8"), T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "t") & blueprint_singler =="cd4 t" & predicted.celltype.l2 |> str_detect("cd4|treg"), T, azhimut_confirmed) ) |>
#     
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "treg") & blueprint_singler %in% c("tregs"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "treg") & predicted.celltype.l2 == "treg", T, azhimut_confirmed) ) |>
#     
#     
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tcm cd4") & blueprint_singler %in% c("cd4 tcm"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tcm cd4") & predicted.celltype.l2 == "cd4 tcm", T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tcm cd8") & blueprint_singler %in% c("cd8 tcm"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tcm cd8") & predicted.celltype.l2 == "cd8 tcm", T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tem cd4") & blueprint_singler %in% c("cd4 tem"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tem cd4") & predicted.celltype.l2 == "cd4 tem", T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tem cd8") & blueprint_singler %in% c("cd8 tem"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tem cd8") & predicted.celltype.l2 == "cd8 tem", T, azhimut_confirmed) ) |>
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tgd") & predicted.celltype.l2 == "gdt", T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "activated cd4") & predicted.celltype.l2 == "cd4 proliferating", T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "activated cd8") & predicted.celltype.l2 == "cd8 proliferating", T, azhimut_confirmed) ) |>
#     
#     
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("naive cd4", "naive t") & predicted.celltype.l2 %in% c("cd4 naive"), T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("naive cd8", "naive t") & predicted.celltype.l2 %in% c("cd8 naive"), T, azhimut_confirmed) ) |>
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("cd4 naive") & !blueprint_singler |> str_detect("clp|hcs|mpp|cd8"), T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("cd8 naive") & !blueprint_singler |> str_detect("clp|hcs|mpp|cd4"), T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("hspc"), T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "prot") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
#     
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "dendritic" & predicted.celltype.l2 %in% c("asdc", "cdc2", "cdc1", "pdc"), T, azhimut_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean == "double negative t regulatory" & predicted.celltype.l2 == "dnt", T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "early t lineage precursor", "immature innate lymphoid") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "early t lineage precursor", "immature innate lymphoid") & predicted.celltype.l2 == "hspc" & blueprint_singler != "clp", T, azhimut_confirmed) ) |>
#     
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("ilc1", "ilc2", "innate lymphoid") & blueprint_singler == "nk", T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("ilc1", "ilc2", "innate lymphoid") & predicted.celltype.l2 %in% c( "nk", "ilc", "nk proliferating"), T, azhimut_confirmed) ) |>
#     
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "immature t") & blueprint_singler %in% c("naive t"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "immature t") & predicted.celltype.l2 == "t naive", T, azhimut_confirmed) ) |>
#     
#     mutate(cell_type_harmonised = if_else(cell_type_clean == "fraction a prepro b", "naive b", cell_type_harmonised))  |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == "granulocyte" & blueprint_singler %in% c("eosinophils", "neutrophils"), T, blueprint_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("immature neutrophil", "neutrophil") & blueprint_singler %in% c( "neutrophils"), T, blueprint_confirmed) ) |>
#     
#     mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("megakaryocyte") & blueprint_singler |> str_detect("megakaryocyte"), T, blueprint_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("macrophage") & blueprint_singler |> str_detect("macrophage"), T, blueprint_confirmed) ) |>
#     
#     mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "nk") & blueprint_singler %in% c("nk"), T, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "nk") & predicted.celltype.l2 %in% c("nk", "nk proliferating", "nk_cd56bright", "ilc"), T, azhimut_confirmed) ) |>
#     
#     
#     # If identical force
#     mutate(azhimut_confirmed = if_else(cell_type_clean == predicted.celltype.l2 , T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean == blueprint_singler , T, blueprint_confirmed) ) |>
#     
#     # Perogenitor
#     mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("progenitor|hematopoietic|precursor") & predicted.celltype.l2  == "hspc", T, azhimut_confirmed) ) |>
#     mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("progenitor|hematopoietic|precursor") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
#     
#     # Generic original annotation and stem for new annotations
#     mutate(azhimut_confirmed = if_else(
#       cell_type_clean  %in% c("T cell", "myeloid cell", "leukocyte", "myeloid leukocyte", "B cell") &
#         predicted.celltype.l2  == "hspc" &
#         blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, azhimut_confirmed) ) |>
#     
#     # Omit mature for stem
#     mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("mature") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), F, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("mature") & predicted.celltype.l2  == "hspc", F, azhimut_confirmed) ) |>
#     
#     # Omit megacariocyte for stem
#     mutate(blueprint_confirmed = if_else(cell_type_clean  == "megakaryocyte" & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), F, blueprint_confirmed) ) |>
#     mutate(azhimut_confirmed = if_else(cell_type_clean  == "megakaryocyte" & predicted.celltype.l2  == "hspc", F, azhimut_confirmed) ) |>
#     
#     # Mast cells
#     mutate(cell_type_harmonised = if_else(cell_type_clean == "mast", "mast", cell_type_harmonised))  |>
#     
#     
#     # Visualise
#     #distinct(cell_type_clean, predicted.celltype.l2, blueprint_singler, strong_evidence, azhimut_confirmed, blueprint_confirmed) |>
#     arrange(!strong_evidence, cell_type_clean) |>
#     
#     # set cell names
#     mutate(cell_type_harmonised = case_when(
#       cell_type_harmonised == "" & azhimut_confirmed ~ predicted.celltype.l2,
#       cell_type_harmonised == "" & blueprint_confirmed ~ blueprint_singler,
#       TRUE ~ cell_type_harmonised
#     )) |>
#     
#     # Add NA
#     mutate(cell_type_harmonised = case_when(cell_type_harmonised != "" ~ cell_type_harmonised)) |>
#     
#     # Add unannotated cells because datasets were too small
#     mutate(cell_type_harmonised = case_when(
#       is.na(cell_type_harmonised) & cell_type_clean  |> str_detect("progenitor|hematopoietic|stem|precursor") ~ "stem",
#       
#       is.na(cell_type_harmonised) & cell_type_clean == "cd14 monocyte" ~ "cd14 mono",
#       is.na(cell_type_harmonised) & cell_type_clean == "cd16 monocyte" ~ "cd16 mono",
#       is.na(cell_type_harmonised) & cell_type_clean %in% c("cd4 cytotoxic t", "tem cd4") ~ "cd4 tem",
#       is.na(cell_type_harmonised) & cell_type_clean %in% c("cd8 cytotoxic t", "tem cd8") ~ "cd8 tem",
#       is.na(cell_type_harmonised) & cell_type_clean |> str_detect("macrophage") ~ "macrophage",
#       is.na(cell_type_harmonised) & cell_type_clean %in% c("mature b", "memory b", "transitional stage b") ~ "b memory",
#       is.na(cell_type_harmonised) & cell_type_clean == "mucosal invariant t" ~ "mait",
#       is.na(cell_type_harmonised) & cell_type_clean == "naive b" ~ "b naive",
#       is.na(cell_type_harmonised) & cell_type_clean == "nk" ~ "nk",
#       is.na(cell_type_harmonised) & cell_type_clean == "naive cd4" ~"cd4 naive",
#       is.na(cell_type_harmonised) & cell_type_clean == "naive cd8" ~"cd8 naive",
#       is.na(cell_type_harmonised) & cell_type_clean == "treg" ~ "treg",
#       is.na(cell_type_harmonised) & cell_type_clean == "tgd" ~ "tgd",
#       TRUE ~ cell_type_harmonised
#     )) |>
#     
#     mutate(confidence_class = case_when(
#       !is.na(cell_type_harmonised) & strong_evidence ~ 2,
#       !is.na(cell_type_harmonised) & !strong_evidence ~ 3
#     )) |>
#     
#     # Lowest grade annotation UNreliable
#     mutate(cell_type_harmonised = case_when(
#       
#       # Get origincal annotation
#       is.na(cell_type_harmonised) & cell_type_clean %in% c("neutrophil", "granulocyte") ~ cell_type_clean,
#       is.na(cell_type_harmonised) & cell_type_clean %in% c("conventional dendritic", "dendritic") ~ "cdc",
#       is.na(cell_type_harmonised) & cell_type_clean %in% c("classical monocyte") ~ "cd14 mono",
#       
#       # Get Seurat annotation
#       is.na(cell_type_harmonised) & predicted.celltype.l2 != "eryth" & !is.na(predicted.celltype.l2) ~ predicted.celltype.l2,
#       is.na(cell_type_harmonised) & !blueprint_singler %in% c(
#         "astrocytes", "smooth muscle", "preadipocytes", "mesangial", "myocytes",
#         "doublet", "melanocytes", "chondrocytes", "mv endothelial", "fibros",
#         "neurons", "keratinocytes", "endothelial", "epithelial", "skeletal muscle", "pericytes", "erythrocytes", "adipocytes"
#       ) & !is.na(blueprint_singler) ~ blueprint_singler,
#       TRUE ~ cell_type_harmonised
#       
#     )) |>
#     
#     # Lowest grade annotation UNreliable
#     mutate(cell_type_harmonised = case_when(
#       
#       # Get origincal annotation
#       !cell_type_harmonised %in% c("doublet", "platelet") ~ cell_type_harmonised
#       
#     )) |>
#     
#     mutate(confidence_class = case_when(
#       is.na(confidence_class) & !is.na(cell_type_harmonised) ~ 4,
#       TRUE ~ confidence_class
#     ))
#   
#   # Another passage
#   
#   # annotated_samples = annotation_crated_UNconfirmed |> filter(!is.na(cell_type_harmonised)) |>  distinct( cell_type, .sample, file_id)
#   #
#   # annotation_crated_UNconfirmed |>
#   # 	filter(is.na(cell_type_harmonised))  |>
#   # 	count(cell_type ,    cell_type_harmonised ,predicted.celltype.l2 ,blueprint_singler) |>
#   # 	arrange(desc(n)) |>
#   # 	print(n=99)
#   
#   
#   annotation_all =
#     annotation_crated_confirmed |>
#     clean_cell_types_deeper() |> 
#     bind_rows(
#       annotation_crated_UNconfirmed
#     ) |>
#     
#     # I have multiple confidence_class per combination of labels
#     distinct() |>
#     with_groups(c(cell_type_clean, predicted.celltype.l2, blueprint_singler), ~ .x |> arrange(confidence_class) |> slice(1)) |>
#     
#     # Simplify after harmonisation
#     mutate(cell_type_harmonised =	case_when(
#       cell_type_harmonised %in% c("b memory", "b intermediate", "classswitched memory b", "memory b" ) ~ "b memory",
#       cell_type_harmonised %in% c("b naive", "naive b") ~ "b naive",
#       cell_type_harmonised %in% c("nk_cd56bright", "nk", "nk proliferating", "ilc") ~ "ilc",
#       cell_type_harmonised %in% c("mpp", "clp", "hspc", "mep", "cmp", "hsc", "gmp") ~ "stem",
#       cell_type_harmonised %in% c("macrophages",  "macrophages m1", "macrophages m2") ~ "macrophage",
#       cell_type_harmonised %in% c("treg",  "tregs") ~ "treg",
#       cell_type_harmonised %in% c("gdt",  "tgd") ~ "tgd",
#       cell_type_harmonised %in% c("cd8 proliferating",  "cd8 tem") ~ "cd8 tem",
#       cell_type_harmonised %in% c("cd4 proliferating",  "cd4 tem") ~ "cd4 tem",
#       cell_type_harmonised %in% c("eosinophils",  "neutrophils", "granulocyte", "neutrophil") ~ "granulocyte",
#       cell_type_harmonised %in% c("cdc",  "cdc1", "cdc2", "dc") ~ "cdc",
#       
#       TRUE ~ cell_type_harmonised
#     )) |>
#     dplyr::select(cell_type_clean, cell_type_harmonised, predicted.celltype.l2, blueprint_singler, confidence_class) |>
#     distinct()
# 
#   
#   curated_annotation =
#     annotation |>
#     clean_cell_types_deeper() |> 
#     filter(lineage_1=="immune") |>
#     dplyr::select(
#       .cell, .sample, cell_type, cell_type_clean, predicted.celltype.l2, blueprint_singler, monaco_singler) |>
#     left_join(
#       annotation_all ,
#       by = c("cell_type_clean", "predicted.celltype.l2", "blueprint_singler")
#     ) |>
#     dplyr::select(
#       .cell, .sample, cell_type, cell_type_harmonised, confidence_class,
#       cell_annotation_azimuth_l2 = predicted.celltype.l2, cell_annotation_blueprint_singler = blueprint_singler,
#       cell_annotation_monaco_singler = monaco_singler
#     ) |>
#     
#     # Reannotation of generic cell types
#     mutate(cell_type_harmonised = case_when(
#       cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("effector memory") ~ "cd4 tem",
#       cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("mait") ~ "mait",
#       cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("central memory") ~ "cd4 tcm",
#       cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("naive") ~ "cd4 naive",
#       cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("effector memory") ~ "cd8 tem",
#       cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("central memory") ~ "cd8 tcm",
#       cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("naive") ~ "cd8 naive",
#       cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler |> str_detect("non classical") ~ "cd16 mono",
#       cell_type == "nonclassical monocyte" & cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="intermediate monocytes"   ~ "cd16 mono",
#       cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler |> str_detect("^classical") ~ "cd14 mono",
#       cell_type == "classical monocyte" & cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="intermediate monocytes"   ~ "cd14 mono",
#       cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="myeloid dendritic" & str_detect(cell_annotation_azimuth_l2, "cdc")   ~ "cdc",
#       
#       
#       TRUE ~ cell_type_harmonised
#     )) |>
#     
#     # Change CD4 classification for version 0.2.1
#     mutate(confidence_class = if_else(
#       cell_type_harmonised |> str_detect("cd4|mait|treg|tgd") & cell_annotation_monaco_singler %in% c("terminal effector cd4 t", "naive cd4 t", "th2", "th17", "t regulatory", "follicular helper t", "th1/th17", "th1", "nonvd2 gd t", "vd2 gd t"),
#       3,
#       confidence_class
#     )) |>
#     
#     # Change CD4 classification for version 0.2.1
#     mutate(cell_type_harmonised = if_else(
#       cell_type_harmonised |> str_detect("cd4|mait|treg|tgd") & cell_annotation_monaco_singler %in% c("terminal effector cd4 t", "naive cd4 t", "th2", "th17", "t regulatory", "follicular helper t", "th1/th17", "th1", "nonvd2 gd t", "vd2 gd t"),
#       cell_annotation_monaco_singler,
#       cell_type_harmonised
#     )) |>
#     
#     
#     mutate(cell_type_harmonised = cell_type_harmonised |>
#              str_replace("naive cd4 t", "cd4 naive") |>
#              str_replace("th2", "cd4 th2") |>
#              str_replace("^th17$", "cd4 th17") |>
#              str_replace("t regulatory", "treg") |>
#              str_replace("follicular helper t", "cd4 fh") |>
#              str_replace("th1/th17", "cd4 th1/th17") |>
#              str_replace("^th1$", "cd4 th1") |>
#              str_replace("nonvd2 gd t", "tgd") |>
#              str_replace("vd2 gd t", "tgd")
#     ) |>
#     
#     # add immune_unclassified
#     mutate(cell_type_harmonised = if_else(cell_type_harmonised == "monocytes", "immune_unclassified", cell_type_harmonised)) |>
#     mutate(cell_type_harmonised = if_else(is.na(cell_type_harmonised), "immune_unclassified", cell_type_harmonised)) |>
#     mutate(confidence_class = if_else(is.na(confidence_class), 5, confidence_class)) |>
#     
#     # drop uncommon cells
#     mutate(cell_type_harmonised = if_else(cell_type_harmonised %in% c("cd4 t", "cd8 t", "asdc", "cd4 ctl"), "immune_unclassified", cell_type_harmonised))
#   
#   
#   # Further rescue of unannotated cells, manually
#   
#   # curated_annotation |>
#   # 	filter(cell_type_harmonised == "immune_unclassified") |>
#   # 	count(cell_type   ,       cell_type_harmonised ,confidence_class ,cell_annotation_azimuth_l2 ,cell_annotation_blueprint_singler ,cell_annotation_monaco_singler) |>
#   # 	arrange(desc(n)) |>
#   # 	write_csv("curated_annotation_still_unannotated_0.2.csv")
#   
#   
#   curated_annotation =
#     curated_annotation |>
#     left_join(
#       read_csv("~/PostDoc/CuratedAtlasQueryR/dev/curated_annotation_still_unannotated_0.2_manually_labelled.csv") |>
#         select(cell_type, cell_type_harmonised_manually_curated = cell_type_harmonised, confidence_class_manually_curated = confidence_class, everything()),
#       by = join_by(cell_type, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler)
#     ) |>
#     mutate(
#       confidence_class = if_else(cell_type_harmonised == "immune_unclassified", confidence_class_manually_curated, confidence_class),
#       cell_type_harmonised = if_else(cell_type_harmonised == "immune_unclassified", cell_type_harmonised_manually_curated, cell_type_harmonised),
#     ) |>
#     select(-contains("manually_curated"), -n) |>
#     
#     # drop uncommon cells
#     mutate(cell_type_harmonised = if_else(cell_type_harmonised %in% c("cd4 tcm", "cd4 tem"), "immune_unclassified", cell_type_harmonised))
#   
#   
#   
#   # # Recover confidence class == 4
#   
#   # curated_annotation |>
#   # 	filter(confidence_class==4) |>
#   # 	count(cell_type   ,       cell_type_harmonised ,confidence_class ,cell_annotation_azimuth_l2 ,cell_annotation_blueprint_singler ,cell_annotation_monaco_singler) |>
#   # 	arrange(desc(n)) |>
#   # 	write_csv("curated_annotation_still_unannotated_0.2_confidence_class_4.csv")
#   
#   curated_annotation =
#     curated_annotation |>
#     left_join(
#       read_csv("~/PostDoc/CuratedAtlasQueryR/dev/curated_annotation_still_unannotated_0.2_confidence_class_4_manually_labelled.csv") |>
#         select(confidence_class_manually_curated = confidence_class, everything()),
#       by = join_by(cell_type, cell_type_harmonised, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler)
#     ) |>
#     mutate(
#       confidence_class = if_else(confidence_class == 4 & !is.na(confidence_class_manually_curated), confidence_class_manually_curated, confidence_class)
#     ) |>
#     select(-contains("manually_curated"), -n)
#   
#   # Correct fishy stem cell labelling
#   # If stem for the study's annotation and blueprint is non-immune it is probably wrong, 
#   # even because the heart has too many progenitor/stem
#   curated_annotation =
#     curated_annotation |>
#     mutate(confidence_class = case_when(
#       cell_type_harmonised == "stem" & cell_annotation_blueprint_singler %in% c(
#         "skeletal muscle", "adipocytes", "epithelial", "smooth muscle", "chondrocytes", "endothelial"
#       ) ~ 5,
#       TRUE ~ confidence_class
#     ))
#   
#   
#   curated_annotation_merged =
#     
#     # Fix cell ID
#     metadata_df |>
#     dplyr::select(.cell, .sample, cell_type) |>
#     as_tibble() |>
#     
#     # Add cell type
#     left_join(curated_annotation |> dplyr::select(-cell_type), by = c(".cell", ".sample")) |>
#     
#     # Add non immune
#     mutate(cell_type_harmonised = if_else(is.na(cell_type_harmonised), "non_immune", cell_type_harmonised)) |>
#     mutate(confidence_class = if_else(is.na(confidence_class) & cell_type_harmonised == "non_immune", 1, confidence_class)) |>
#     
#     # For some unknown reason
#     distinct()
#   
#   
#   curated_annotation_merged |>
#     
#     # Save
#     saveRDS(file_curated_annotation_merged)
#   
#   metadata_annotated =
#     curated_annotation_merged |>
#     
#     # merge with the rest of metadata
#     left_join(
#       metadata_df |>
#         as_tibble(),
#       by=c(".cell", ".sample", "cell_type")
#     )
#   
#   # Replace `.` with `_` for all column names as it can create difficoulties for MySQL and Python
#   colnames(metadata_annotated) = colnames(metadata_annotated) |> str_replace_all("\\.", "_")
#   metadata_annotated = metadata_annotated |> rename(cell_ = `_cell`, sample_ = `_sample`)
#   
# 
#   dictionary_connie_non_immune = 
#     metadata_annotated |> 
#     filter(cell_type_harmonised == "non_immune") |> 
#     distinct(cell_type) |> 
#     harmonise_names_non_immune() |> 
#     rename(cell_type_harmonised_non_immune = cell_type_harmonised )
#   
#   metadata_annotated = 
#     metadata_annotated |> 
#     left_join(dictionary_connie_non_immune) |> 
#     mutate(cell_type_harmonised = if_else(cell_type_harmonised=="non_immune", cell_type_harmonised_non_immune, cell_type_harmonised)) |> 
#     select(-cell_type_harmonised_non_immune)
#   
#   
# }

remove_files_safely <- function(files) {
  for (file in files) {
    if (file.exists(file)) {
      file.remove(file)
    }
  }
}


#' Get positions of each unique element in a vector
#'
#' This function takes a vector and returns a named list where each unique
#' element of the input vector maps to the positions at which it occurs.
#'
#' @param input_vector A vector of elements.
#' @return A named list where each name is a unique element from the input vector
#' and each value is a vector of positions where that element occurs.
#' @examples
#' input_vector <- c("a", "a", "b", "c", "a")
#' positions_list <- get_positions(input_vector)
#' print(positions_list)
#' @importFrom dplyr tibble
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom purrr set_names
#' @export
get_positions <- function(input_vector) {
  # Create a tibble with the input vector and their positions
  df <- tibble(value = input_vector, position = seq_along(input_vector))
  
  # Group by value and summarise the positions
  result <- df %>%
    group_by(value) %>%
    summarise(positions = list(position), .groups = 'drop')
  
  # Convert the result to a named list
  result_list <- set_names(result$positions, result$value)
  
  return(result_list)
}

#' Convert a vector of integers to its R code form as a character
#'
#' This function takes a vector of integers and returns its R code form as a character string.
#'
#' @param int_vector A vector of integers.
#' @return A character string representing the R code form of the vector.
#' @examples
#' int_vector <- c(1, 2, 3, 4, 5)
#' code_string <- vector_to_code(int_vector)
#' print(code_string)
#' @export
vector_to_code <- function(int_vector) {
  # Convert the vector to its R code form
  code_string <- deparse(int_vector)
  
  # Collapse the output into a single string
  code_string <- paste(code_string, collapse = "")
  
  return(code_string)
}

#' Add Tier Inputs to a Function Call String
#'
#' This function modifies a function call string by appending a tier label to specified arguments.
#'
#' @param command A character string representing a function call, e.g., "a(b, c, d)".
#' @param arguments_to_tier A character vector specifying which arguments should be tiered, e.g., c("b", "c", "d").
#' @param i A character string representing the tier label to be appended, e.g., "_1".
#'
#' @return A character string representing the modified function call with tiered arguments.
#'
#' @importFrom stringr str_subset
#' @importFrom stringr str_replace_all
#' @importFrom glue glue
#'
#' @examples
#' # Example usage:
#' command <- "a(b, c, d)"
#' arguments_to_tier <- c("b", "c", "d")
#' i <- "_1"
#' output <- add_tier_inputs(command, arguments_to_tier, i)
#' print(output)  # Outputs: "a(b_1, c_1, d_1)"
#'
#' @noRd
add_tier_inputs <- function(command, arguments_to_tier, i) {
  
  if(i |> length() > 1) stop("HPCell says: argument i must be of length one")
  
  if(length(arguments_to_tier)==0) return(command)
  
  command = command |> deparse() |> paste(collapse = "")  
  
  # Filter out arguments to be tiered from the input command
  #input = command |> str_extract("[a-zA-Z0-9_]+\\(([a-zA-Z0-9_]+),.*", group=1) 
  #arguments_to_tier <- arguments_to_tier |> str_subset(input, negate = TRUE)
  
  # Create a named vector for replacements
  replacement_regexp <- glue("{arguments_to_tier}_{i}") |> as.character() |> set_names(arguments_to_tier)
  
  # Function to add word boundaries and perform the replacements
  # This because we only replace WHOLE words
  add_word_boundaries_and_replace <- function(command, replacements) {
    for (pattern in names(replacements)) {
      # Create the regex pattern with word boundaries
      pattern_with_boundaries <- paste0("\\b", pattern, "\\b")
      # Perform the replacement for each pattern
      command <- str_replace(command, pattern_with_boundaries, replacements[pattern])
    }
    return(command)
  }
  
  # Replace the specified arguments in the command with their tiered versions
  command |> add_word_boundaries_and_replace(replacement_regexp) |>  rlang::parse_expr()
  
}

#' Divide Features into Chunks
#'
#' This function divides a list of features into chunks of a specified size.
#'
#' @param features A vector of features to be divided into chunks.
#' @param chunk_size The size of each chunk. Defaults to 100.
#' @return A tibble with the features and their corresponding chunk numbers.
#' @importFrom dplyr tibble
#' @importFrom rlang rep_along
#' @importFrom magrittr divide_by
#' 
#' 
#' @export
feature_chunks = function(features, chunk_size = 100){
  
  chunks = 
    features |> 
    length() |> 
    divide_by(chunk_size) |> 
    ceiling() |> 
    seq_len() |>
    rep(each = chunk_size, length.out = length(features))
  
  
  tibble(feature = features, chunk___ = chunks) 
  
}



add_missingh_genes_to_se = function(se, all_genes, missing_genes){
  
  missing_matrix = matrix(rep(0, length(missing_genes) * ncol(se)), ncol = ncol(se))
  
  rownames(missing_matrix) = missing_genes
  colnames(missing_matrix) = colnames(se)
  
  new_se = SummarizedExperiment(assays = list(count = missing_matrix |> DelayedArray::DelayedArray() ),
                                colData = colData(se))
  
  
  empty_rowdata = DataFrame(matrix(NA, ncol = ncol(rowData(se)), nrow = length(missing_genes)),
                            row.names = missing_genes)
  names(empty_rowdata) <- names(rowData(se))
  rowData(new_se) = empty_rowdata
  
  se = SummarizedExperiment(assays = assays(se), colData = colData(se), rowData = rowData(se))
  se = se |> rbind(new_se)
  
  se[all_genes,]
  
}



#' Remove Random Effects from a Mixed Effects Model Formula
#'
#' This function takes a formula used in mixed effects models and removes the random effects part,
#' returning the fixed effects formula.
#'
#' @param formula A formula used in mixed effects models, containing both fixed and random effects.
#' @return A formula containing only the fixed effects part of the input formula.
#' @importFrom lme4 findbars
#' @importFrom stats update
#' @importFrom stats as.formula
#' @examples
#' library(lme4)
#' f1 <- Reaction ~ Days + (Days | Subject)
#' fixed_formula <- remove_random_effects(f1)
#' print(fixed_formula)
#' @noRd
remove_random_effects <- function(formula) {
  # Find the random effects parts of the formula
  random_effects <- findbars(formula)
  
  # Start with the original formula
  fixed_formula <- formula
  
  # Remove each random effects term from the formula
  for (re in random_effects) {
    fixed_formula <- update(fixed_formula, paste(". ~ . - (", deparse(re), ")", sep = ""))
  }
  
  # Convert the updated formula to a character string
  fixed_formula_str <- deparse(fixed_formula)
  
  # Remove the leading ". ~" if present and ensure the string is properly formatted
  fixed_formula_str <- sub("^\\. ~", "~", fixed_formula_str)
  
  # Convert the string back to a formula object
  fixed_formula <- as.formula(fixed_formula_str)
  
  # Return the fixed effects formula
  return(fixed_formula)
}

#' @examples
#' # Example usage:
#' # delete_lines_with_word("your_text_file.txt", "bla")
#'
#' @export
delete_lines_with_word <- function(word, file_path) {
  # Step 1: Read the file into R as a vector of lines
  lines <- readLines(file_path)
  
  # Step 2: Filter out lines that contain the specified word
  word = glue("target_output = \"{word}\"")
  filtered_lines <- lines[!grepl(word, lines)]
  
  # Step 3: Write the modified content back to the file
  writeLines(filtered_lines, file_path)
}

#' Get elements with class 'name'
#'
#' This function takes a list and returns a character vector of elements
#' that have the class 'name', converting them to their character equivalents.
#'
#' @param lst A list of elements to process.
#' @return A character vector of elements with class 'name'.
#' @noRd
get_elements_with_name_class <- function(lst) {
  lapply(lst, function(x) {
    if ("name" %in% class(x)) as.character(x) else NULL
  }) %>%
    Filter(Negate(is.null), .) %>% # Remove NULL elements
    unlist()
}

#' Get Arguments to Tier Based on Iteration Settings
#'
#' This function identifies elements from a list that have the class 'name',
#' converts them to character strings, and returns only those elements that are
#' present in the names of a specified input list (`input_hpc`) and have the
#' `iterate` field set to `"tier"`.
#'
#' @param lst A list containing various elements, some of which may have the class 'name'.
#' @param input_hpc A list whose names are used to filter the elements from `lst`.
#'                  The elements in `input_hpc` should include an `iterate` field with the value `"tier"`.
#' @return A character vector of elements from `lst` that have the class 'name',
#'         are present in the names of `input_hpc`, and have `iterate` set to `"tier"`.
#' @noRd
arguments_to_action <- function(lst, input_hpc, value) {
  matching_elements <- character()
  
  for (arg_name in names(lst)) {
    arg_value <- lst[[arg_name]]
    
    # Skip NULL and complex values, because they cannot be a name of a target
    if (
      arg_value |> length() == 0 |
      is.null(arg_value) | !(
        arg_value |> is("character") | 
        arg_value |> is("name") | 
        arg_value |> is("list")
      )) next
    
    # Convert the argument value to a character string vector
    # arg_value_as_char <- as.character(arg_value)
    
    # This because I cannot loop over a single "name" class, while I can loop over a single "character" class
    # NOT SO ELEGANT 
    if(arg_value |> length() == 1){
      
      if (
        (arg_value |> is("character") | arg_value |> is("name")) &&
        as.character(arg_value) %in% names(input_hpc) && 
        input_hpc[[arg_value]]$iterate %in% value
      ) 
        matching_elements <- c(matching_elements, as.character(arg_value) |> set_names(arg_name))
      
    }
    
    else{
      # Iterate over each element in arg_value_as_char
      for (val in arg_value) {
        
        # Again, skip if the list include complex elements
        if (is.null(arg_value) | !(
          arg_value |> is("character") | 
          arg_value |> is("name") 
        )) next
        
        # Check if the value exists in input_hpc and iterate is equal to the specified value
        if (val %in% names(input_hpc) && input_hpc[[val]]$iterate == value) 
          matching_elements <- c(matching_elements, as.character(arg_value) |> set_names(arg_name))
        
      }
    }
    
  }
  
  return(matching_elements)
}



#' Quote elements with class 'name'
#'
#' This function takes a list and returns a new list where any elements
#' with the class 'name' are converted to their quoted equivalent using `quote()`.
#' This is useful for preserving unevaluated expressions in the list.
#'
#' @param lst A list of elements to process.
#' @return A list where elements with class 'name' are quoted.
#' @noRd
quote_name_classes <- function(lst) {
  lapply(lst, function(x) {
    if ("name" %in% class(x)) {
      # Manually create the quoted expression
      as.call(list(as.name("quote"), x))
    } else {
      x  # Leave as is for other elements
    }
  })
}

#' Safe as.name Wrapper
#'
#' This function wraps `as.name()` to safely handle `NULL` input.
#' If the input is `NULL`, the function returns `NULL`; otherwise,
#' it returns the result of `as.name()`.
#'
#' @param input The input to be converted to a name. If `NULL`, the function returns `NULL`.
#' @return The result of `as.name()` applied to the input, or `NULL` if the input is `NULL`.
#' @noRd
safe_as_name <- function(input) {
  if (is.null(input)) {
    return(NULL)
  } else {
    return(as.name(input))
  }
}

#' Check for Name-Value Conflicts in Arguments
#'
#' This function checks if any argument names in a given function call are identical
#' to any of their corresponding values. If such a conflict is found, an error is thrown.
#' This validation step is crucial for ensuring that arguments do not unintentionally
#' share the same name as their value, which could lead to unexpected behavior or errors
#' in downstream processes.
#'
#' @param ... Arguments to be checked for name-value conflicts.
#' @return The function returns the input arguments as a list if no conflicts are found.
#' @details 
#' The `check_for_name_value_conflicts()` function is designed to catch cases where the name of an argument matches one of its values. For example, if you pass an argument like `sample_id = "sample_id"`, this function will detect that the name and value are identical and throw an error. Such conflicts can cause confusion or unintended behavior in your code, especially in complex workflows or pipelines where argument names are often used to identify specific data or parameters.
#' 
#' The function works by iterating over all arguments passed via `...`, converting each argument's value to a character string, and then checking if the argument's name appears in this string. If a match is found, an error is raised with a clear message indicating the problematic argument.
#' 
#' This function is particularly useful in contexts like data processing pipelines where arguments may be dynamically generated or modified. Ensuring that no argument name matches its value helps maintain clarity and prevent errors in such scenarios.
#'
#' @examples
#' check_for_name_value_conflicts(sample_id = "sample_001", group = "control")
#' 
#' # This will throw an error:
#' # check_for_name_value_conflicts(sample_id = "sample_id")
#' 
#' @importFrom glue glue
#' @noRd
check_for_name_value_conflicts <- function(...) {
  # Capture the arguments passed to the function
  args_list <- list(...)
  
  # Iterate through the list and check for name-value conflicts
  for (arg_name in names(args_list)) {
    arg_value <- args_list[[arg_name]]
    
    # Skip NULL values
    if (is.null(arg_value)) next
    
    # Convert the argument value to a character string
    arg_value_as_char <- as.character(arg_value)
    
    # Check if the argument name matches any of the values in arg_value_as_char
    if (arg_name %in% c(arg_value)) {
      stop(glue::glue("HPCell says: Argument name '{arg_name}' cannot be the same as its value '{arg_value_as_char}'"))
    }
  }
  
  # If no conflicts, return the arguments as is or proceed with the function logic
  return(args_list)
}

#' Expand Tiered Arguments in a List
#'
#' This function takes a list of arguments (`lst`), identifies a specific argument to replace (`argument_to_replace`),
#' and expands it into a list of quoted tiered values. This is particularly useful when you need to dynamically generate
#' tiered versions of an argument within a list structure.
#'
#' @param lst A list of arguments where one argument will be replaced by a list of tiered versions.
#' @param tiers A vector of tiers (e.g., `c("1", "2")`) used to generate the tiered versions of the argument.
#' @param argument_to_replace The name of the argument in `lst` that should be replaced by the tiered versions.
#' @param tiered_args The base name used to create the tiered versions. The tiers will be appended to this base name.
#' @return The modified list where the specified argument is replaced by a list of quoted tiered values.
#' @details 
#' The `expand_tiered_arguments()` function is designed to dynamically generate tiered versions of an argument
#' in a list. For example, if you have an argument `pseudobulk_list` in `lst` that you want to replace with tiered
#' versions like `pseudobulk_se_merge_within_tier_1` and `pseudobulk_se_merge_within_tier_2`, this function will 
#' create those versions, quote them (to prevent evaluation), and replace the original argument in `lst` with a 
#' list of these quoted expressions. This is useful in contexts where arguments need to be programmatically 
#' generated and passed to functions that expect lists of unevaluated symbols.
#'
#' The function works by first checking if the `argument_to_replace` exists in the `lst`. If it does, it constructs 
#' the tiered versions by iterating over the `tiers` vector, creating symbols for each tier, and wrapping each 
#' symbol in a `quote()` to prevent immediate evaluation. The result is a list of quoted expressions that replace 
#' the original argument in `lst`.
#' 
#' @examples
#' args_list <- list(
#'   external_path = "_targets/external",
#'   pseudobulk_list = "pseudobulk_se_iterated",
#'   packages = c("tidySummarizedExperiment", "HPCell")
#' )
#' 
#' name_target_intermediate <- "pseudobulk_se_merge_within_tier"
#' 
#' result <- expand_tiered_arguments(
#'   lst = args_list, 
#'   tiers = c("1", "2"), 
#'   argument_to_replace = "pseudobulk_list",
#'   tiered_args = name_target_intermediate
#' )
#' 
#' # The output will be:
#' # $external_path
#' # [1] "_targets/external"
#' #
#' # $pseudobulk_list
#' # list(quote(pseudobulk_se_merge_within_tier_1), quote(pseudobulk_se_merge_within_tier_2))
#' #
#' # $packages
#' # [1] "tidySummarizedExperiment" "HPCell"
#' 
#' @noRd
expand_tiered_arguments <- function(lst, tiers, argument_to_replace, tiered_args) {
  # Check if the argument to replace exists in the list
  if (argument_to_replace %in% names(lst)) {
    # Fetch the correct value of the tiered argument from the list
    tiered_base <- lst[[argument_to_replace]]
    
    # Create a vector of tiered values by combining tiered_base with tiers
    # If no tier do not add the suffix
    tiered_values <- lapply(tiers, function(tier) paste0(tiered_base, "_", tier) |> as.name() )
    
    # Construct the c(...) call with the tiered values
    c_call <- as.call(c(as.name("c"), tiered_values))
    
    # Wrap the entire c(...) call with quote(quote(...))
    lst[[argument_to_replace]] <- substitute(quote(quote(expr)), list(expr = c_call))
  }
  
  return(lst)
}


build_pattern = function(arguments_to_tier = c(), other_arguments_to_map = c(), index = c()){
  
  pattern = NULL 
  
  if(
    arguments_to_tier |> length() > 0 |
    other_arguments_to_map |> length() > 0
  ){
    
    pattern = as.name("map")
    
    if(arguments_to_tier |> length() > 0)
      pattern = pattern |> c(
        arguments_to_tier |>
          map(
            ~ substitute(
              slice(input, index  = arg ), 
              list(input = as.symbol(.x), arg=index)
            ) 
          )
      )
    
    if(other_arguments_to_map |> length() > 0){
      
      pattern = pattern |> c(other_arguments_to_map |> lapply(as.name))
      
    }
    
    pattern = as.call(pattern)
    
  }
  
  pattern
  
}

write_source = function(user_function_source_path, target_script){
  if(user_function_source_path |> is.null() |> not())
    
    source(s) |> 
    substitute(env = list(s =user_function_source_path )) |> 
    deparse() |> 
    write_lines(target_script, append = TRUE)
}

#' @export
target_append <- function(target_list, ...) {
  # Append the new elements to the list
  target_list <<- c(target_list, list(...))
  
}

write_HDF5_array_safe = function(normalized_rna, name, directory){
  
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
  
  hash = digest(normalized_rna)
  file_name = glue("{directory}/{hash}")
  
  if (
    file.exists(file_name) &&
    name %in% rhdf5::h5ls(file_name)$name
  ) {
    names_to_drop = rhdf5::h5ls(file_name)$name |> str_subset(name)
    names_to_drop |> map(~rhdf5::h5delete(file_name, .x))
  }
  
  normalized_rna |> 
    HDF5Array::writeHDF5Array(
      filepath = file_name,
      name = name,
      as.sparse = TRUE
    ) 
  
}


#' Compute the Mode of a DelayedArray
#'
#' This function computes the mode (most frequent value) of a \code{DelayedArray} without loading the entire array into memory. It processes the array in blocks to maintain memory efficiency, making it suitable for large datasets.
#'
#' @param delayed_array A \code{DelayedArray} object for which the mode is to be computed.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{\code{mode}}{Numeric vector of the most frequent value(s) in the array.}
#'   \item{\code{frequency}}{Integer representing the count of the most frequent value(s).}
#' }
#'
#' @details
#' The function utilizes block processing via \code{blockApply()} from the \code{DelayedArray} package to avoid loading the entire array into memory. It computes partial frequency tables for each block and combines them to find the overall mode.
#'
#' @import DelayedArray
#' @importFrom DelayedArray blockApply
#' @importFrom methods as
#' @importFrom utils capture.output
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(DelayedArray)
#'
#' # Create a DelayedArray from an in-memory matrix
#' set.seed(123)
#' n_rows <- 1000
#' n_cols <- 1000
#' matrix_data <- matrix(sample(0:5, n_rows * n_cols, replace = TRUE,
#'                              prob = c(0.5, 0.1, 0.1, 0.1, 0.1, 0.1)),
#'                       nrow = n_rows)
#' delayed_array <- DelayedArray(matrix_data)
#'
#' # Compute the mode
#' mode_result <- compute_mode_delayedarray(delayed_array)
#'
#' # Output the result
#' cat("Most frequent value(s):", paste(mode_result$mode, collapse = ", "), "\n")
#' cat("Frequency:", mode_result$frequency, "\n")
#' }
#'
#' @export
compute_mode_delayedarray <- function(delayed_array) {
  # Compute the mode (most frequent value) of a DelayedArray without loading the entire array into memory.
  
  # Helper function to compute frequency table for a block
  block_table <- function(block) {
    counts_vector <- as.vector(block)
    counts_table <- table(counts_vector)
    return(counts_table)
  }
  
  # Helper function to combine two frequency tables
  combine_tables <- function(table1, table2) {
    if (length(table1) == 0) return(table2)
    if (length(table2) == 0) return(table1)
    
    # Get all unique values
    all_names <- union(names(table1), names(table2))
    
    # Align counts for all unique values
    counts1 <- as.numeric(table1[all_names])
    counts2 <- as.numeric(table2[all_names])
    
    # Replace NA with 0
    counts1[is.na(counts1)] <- 0
    counts2[is.na(counts2)] <- 0
    
    # Sum counts
    combined_counts <- counts1 + counts2
    names(combined_counts) <- all_names
    
    return(combined_counts)
  }
  
  # Process the DelayedArray in blocks
  block_tables <- blockApply(delayed_array, FUN = block_table)
  
  # Combine all partial frequency tables into a single table
  freq_counts <- Reduce(f = combine_tables, x = block_tables)
  
  # Find the value(s) with the maximum count
  max_count <- max(freq_counts)
  most_frequent_values <- as.numeric(names(freq_counts)[freq_counts == max_count])
  
  # Return the result as a list
  result <- list(
    mode = most_frequent_values,
    frequency = max_count
  )
  
  return(result)
}

#' Check if All Assay Values are Greater Than Zero and Subtract One if True
#'
#' This function, `check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY`, checks if all values 
#' in a specified assay of a `SingleCellExperiment` or `Seurat` object are greater than zero. 
#' If all values are greater than zero, it subtracts one from each value. This operation is useful 
#' in cases where a small adjustment to count data is necessary to standardise the minimum count value.
#' 
#' For `SingleCellExperiment` objects, the function accesses the assay data using the `assay` function. 
#' For `Seurat` objects, it retrieves the data using `GetAssayData` and updates it using `SetAssayData`. 
#' This allows seamless handling of different object types in single-cell analysis workflows.
#'
#' @param input_read_RNA_assay A `SingleCellExperiment` or `Seurat` object containing the assay data.
#' @param assay_name A string specifying the name of the assay to be checked and potentially modified.
#'
#' @return The modified `SingleCellExperiment` or `Seurat` object, where one has been subtracted 
#'         from all values in the specified assay if all values were initially greater than zero. 
#'         If any values are zero or negative, the object is returned unmodified.
#'         
#' @examples
#' # For SingleCellExperiment
#' # sce <- SingleCellExperiment(assays = list(RNA = matrix(1:9, 3, 3)))
#' # modified_sce <- check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(sce, "RNA")
#'
#' # For Seurat
#' # seurat <- CreateSeuratObject(counts = matrix(1:9, 3, 3))
#' # modified_seurat <- check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(seurat, "RNA")
#'
#' @import SingleCellExperiment
#' @import Seurat
#' @importFrom SummarizedExperiment assay
#' 
#' @noRd
check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY <- function(input_read_RNA_assay, assay_name, subset_up_to_number_of_cells = dim(input_read_RNA_assay)[2]) {
  
  # Do now overshoor the number of cells
  subset_up_to_number_of_cells = subset_up_to_number_of_cells |> min(dim(input_read_RNA_assay)[2])
  
  # Check if object is SCE or Seurat
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    # For SingleCellExperiment
    assay_data <- assay(input_read_RNA_assay, assay_name)
    
    my_min = min(assay_data[, 1:subset_up_to_number_of_cells])
    
    # Check if all values are > 0
    if (my_min > 0) {
      # Subtract 1 from each value
      assay(input_read_RNA_assay, assay_name) <- assay_data - my_min
    } else {
      message("Not all values are greater than 0. No subtraction performed.")
    }
    
  } else if (inherits(input_read_RNA_assay, "Seurat")) {
    # For Seurat
    assay_data <- GetAssayData(input_read_RNA_assay, assay = assay_name)
    
    my_min = min(assay_data)
    
    # Check if all values are > 0
    if (my_min > 0) {
      # Subtract 1 from each value
      input_read_RNA_assay <- SetAssayData(input_read_RNA_assay, assay = assay_name, 
                                           new.data = assay_data - my_min)
    } else {
      message("Not all values are greater than 0. No subtraction performed.")
    }
    
  } else {
    stop("The input object is neither a SingleCellExperiment nor a Seurat object.")
  }
  
  # Return the modified object
  return(input_read_RNA_assay)
}

#' Clean Metadata in SingleCellExperiment Object
#'
#' This function takes a SingleCellExperiment (SCE) object and removes columns 
#'    that are completely filled with NA values.
#' The cleaned metadata is then returned as a dataframe.
#'
#' @param sce A SingleCellExperiment object containing metadata to be cleaned.
#' @return A SingleCellExperiment with all completely NA columns removed
clean_sce_metadata <- function(sce) {
  sce <- sce |> select(where(~ any(!is.na(.))))
  sce
}

#' @export
computeCommunProbPathway <- function(object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) {
  if (is.null(net)) {
    net <- object@net
  }
  if (is.null(pairLR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  prob <- net$prob
  prob[net$pval > thresh] <- 0
  
  LR <- dimnames(prob)[[3]]
  LR.sig <- LR[apply(prob, 3, sum) != 0]
  
  pathways <- unique(pairLR.use$pathway_name)
  group <- factor(pairLR.use$pathway_name, levels = pathways)
  
  # STEFANO FIX
  if(length(levels(group))==1){
    xx = apply(prob, c(1, 2), by, group, sum)
    prob.pathways = xx |> array(dim = c(nrow(xx), ncol(xx), 1), dimnames = list(rownames(xx), colnames(xx), levels(group)))
  }
  else
    prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum),
                           c(2, 3, 1))
  
  pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 0]
  prob.pathways.sig <- prob.pathways[,,pathways.sig, drop = FALSE]
  idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing=TRUE, index.return = TRUE)$ix
  pathways.sig <- pathways.sig[idx]
  prob.pathways.sig <- prob.pathways.sig[, , idx]
  
  if (is.null(object)) {
    netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
    return(netP)
  } else {
    object@net$LRs <- LR.sig
    object@netP$pathways <- pathways.sig
    object@netP$prob <- prob.pathways.sig
    return(object)
  }
}
