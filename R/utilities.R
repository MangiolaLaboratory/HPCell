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

#' Identify Empty Droplets in Single-Cell RNA-seq Data
#'
#' @description
#' `empty_droplet_id` distinguishes between empty and non-empty droplets using the DropletUtils package.
#' It excludes mitochondrial and ribosomal genes, calculates barcode ranks, and optionally filters input data
#' based on these criteria. The function returns a tibble containing log probabilities, FDR, and a classification
#' indicating whether cells are empty droplets.
#'
#' @param input_read_RNA_assay SingleCellExperiment object containing RNA assay data.
#' @param filter_empty_droplets Logical value indicating whether to filter the input data.
#'
#' @return A tibble with columns: logprob, FDR, empty_droplet (classification of droplets).
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stringr str_subset
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom tidyr replace_na
#' @importFrom DropletUtils emptyDrops
#' @importFrom DropletUtils barcodeRanks
#' @importFrom S4Vectors metadata
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @noRd
empty_droplet_id <- function(input_read_RNA_assay,
                             filter_empty_droplets,
                             assay = NULL){
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  
  significance_threshold = 0.001
  # Genes to exclude
  location <- mapIds(
    EnsDb.Hsapiens.v86,
    keys=rownames(input_read_RNA_assay),
    column="SEQNAME",
    keytype="SYMBOL"
  )
  mitochondrial_genes = which(location=="MT") |> names()
  ribosome_genes = rownames(input_read_RNA_assay) |> str_subset("^RPS|^RPL")
  
  # if ("originalexp" %in% names(input_file@assays)) {
  #   barcode_ranks <- barcodeRanks(input_file@assays$originalexp@counts[!rownames(input_file@assays$originalexp@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE])
  # } else if ("RNA" %in% names(input_file@assays)) {
  #   barcode_ranks <- barcodeRanks(input_file@assays$RNA@counts[!rownames(input_file@assays$RNA@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE])
  # }
  
  # Calculate bar-codes ranks
  barcode_ranks = barcodeRanks(GetAssayData(input_read_RNA_assay, assay, slot = "counts")[!rownames(GetAssayData(input_read_RNA_assay, assay, slot = "counts")) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE])
  
  # Set the minimum total RNA per cell for ambient RNA
  if(min(barcode_ranks$total) < 100) { lower = 100 } else {
    lower = quantile(barcode_ranks$total, 0.05)
    
    # write_lines(
    #   glue("{input_path} has supposely empty droplets with a lot of RNAm maybe a lot of ambient RNA? Please investigate"),
    #   file = glue("{dirname(output_path_result)}/warnings_emptyDrops.txt"),
    #   append = T
    # )
  }
  
  # Remove genes from input
  if (
    # If filter_empty_droplets
    filter_empty_droplets == "TRUE") {
    barcode_table <- GetAssayData(input_read_RNA_assay, assay, slot = "counts")[!rownames(GetAssayData(input_read_RNA_assay, assay, slot = "counts")) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE] |>
      emptyDrops( test.ambient = TRUE, lower=lower) |>
      as_tibble(rownames = ".cell") |>
      mutate(empty_droplet = FDR >= significance_threshold) |>
      replace_na(list(empty_droplet = TRUE))
  }
  else {
    barcode_table <- select(., .cell) |>
      as_tibble() |>
      mutate( empty_droplet = FALSE)
  } 
  
  # barcode ranks
  barcode_table <- barcode_table |>
    left_join(
      barcode_ranks |>
        as_tibble(rownames = ".cell") |>
        mutate(
          knee =  metadata(barcode_ranks)$knee,
          inflection =  metadata(barcode_ranks)$inflection
        )
    )
  
  
  # barcode_table |>  saveRDS(output_path_result)
  
  # # Plot bar-codes ranks
  # plot_barcode_ranks =
  #   barcode_table %>%
  #   ggplot2::ggplot(aes(rank, total)) +
  #   geom_point(aes(color = empty_droplet, size = empty_droplet )) +
  #   geom_line(aes(rank, fitted), color="purple") +
  #   geom_hline(aes(yintercept = knee), color="dodgerblue") +
  #   geom_hline(aes(yintercept = inflection), color="forestgreen") +
  #   scale_x_log10() +
  #   scale_y_log10() +
  #   scale_color_manual(values = c("black", "#e11f28")) +
  #   scale_size_discrete(range = c(0, 2)) +
  #   theme_bw()
  
  # plot_barcode_ranks |> saveRDS(output_path_plot_rds)
  
  # ggsave(
  #   output_path_plot_pdf,
  #   plot = plot_barcode_ranks,
  #   useDingbats=FALSE,
  #   units = c("mm"),
  #   width = 183/2 ,
  #   height = 183/2,
  #   limitsize = FALSE
  # )
  
  barcode_table
  # return(list(barcode_table, plot_barcode_ranks))
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
#' @noRd
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
#' @noRd
reference_label_coarse_id <- function(tissue) {
  return(
    ifelse(tissue == "pbmc", "monaco_first.labels.coarse",
           ifelse(tissue == "solid", "blueprint_first.labels.coarse",
                  ifelse(tissue == "atypical", "none",
                         ifelse(tissue == "none", "monaco_first.labels.coarse", NA)))))
}

#' Change Default Assay to RNA
#'
#' @description
#' `add_RNA_assay` changes the default assay in a Seurat object to RNA.
#'
#' @param input_read Seurat object.
#' @param RNA_assay_name Name of the RNA assay to set as default.
#'
#' @return Modified Seurat object with the default assay set to RNA.
#'
#' @importFrom Seurat DefaultAssay
#' @export
#' @noRd
add_RNA_assay <- function(input_read, RNA_assay_name){
  
  if (RNA_assay_name != "RNA"){
    input_read[["RNA"]] = input_read[[RNA_assay_name]]
    DefaultAssay(object = input_read) <- "RNA"
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
#' @importFrom glue glue
#' @noRd
seurat_to_variable_features_by_cell_type = function(counts, assay, .cell_group = NULL, features_number_per_cell_group = 300){
  
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
#' @noRd
seurat_to_variable_features = function(
    counts,
    assay,
    .sample,
    .cell_group = NULL,
    features_number_independent_of_cell_groups = 300,
    features_number_per_cell_group = 300
){
  
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
#' @noRd
subset_top_rank_variable_genes_across_batches = function(
    table_across_cell_groups,
    table_within_cell_groups,
    .cell_group,
    .batch,
    features_number_independent_of_cell_groups = 2000,
    features_number_per_cell_group = 300
){
  
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
#' @description
#' Identify variable features, scale the data, run Principal Component Analysis (PCA), find neighbors, identify clusters,
#' and compute UMAP (Uniform Manifold Approximation and Projection) embeddings for 
#' visualization of cell clusters. 
#'
#' @noRd
#' 
calc_UMAP <- function(input_seurat){
  assay_name = input_seurat@assays |> names() |> extract2(1)
  find_var_genes <- FindVariableFeatures(input_seurat)
  var_genes<- find_var_genes@assays[[assay_name]]@var.features
  
  x<- ScaleData(input_seurat) |>
    # Calculate UMAP of clusters
    RunPCA(features = var_genes) |>
    FindNeighbors(dims = 1:30) |>
    FindClusters(resolution = 0.5) |>
    RunUMAP(dims = 1:30, spread    = 0.5,min.dist  = 0.01, n.neighbors = 10L) |> 
    as_tibble()
  return(x)
}
#' Subsetting input dataset into a list of seurat objects by sample/ tissue 
#'  @description
#' Function to subset Seurat object by tissue
get_unique_tissues <- function(seurat_object, sample_column) {
  sample_column<- quo_name(sample_column)
  unique_sample <- seurat_object@meta.data |> pull(sample_column) |> unique()
  
  return(unique_sample)
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

#' Clean and Standardize Cell Types (Deeper)
#'
#' This function takes a vector of cell types and applies a series of transformations
#' to clean and standardize them for better consistency.
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#'
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_remove
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_trim
#'
#' @param x A vector of cell types.
#'
#' @return A cleaned and standardized vector of cell types.
#'
# @examples
# cell_types <- c("CD4 T Cell, AlphaBeta", "NK cell, gammadelta", "Central Memory")
# cleaned_cell_types <- clean_cell_types_deeper(cell_types)

clean_cell_types_deeper = function(x){
  x |> 
    # Annotate
    mutate(cell_type_clean = cell_type_clean |> tolower()) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove_all(",")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("alphabeta")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove_all("positive")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("cd4  t", "cd4")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("regulatory t", "treg")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("thymusderived")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("human")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("igg ")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("igm ")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("iga ")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("group [0-9]")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("common")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("cd45ro")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("type i")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("germinal center")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("iggnegative")) |>
    mutate(cell_type_clean = cell_type_clean |> str_remove("terminally differentiated")) |>
    
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect("macrophage"), "macrophage", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean == "mononuclear phagocyte", "macrophage", cell_type_clean) ) |>
    
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect(" treg"), "treg", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect(" dendritic"), "dendritic", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect(" thelper"), "thelper", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect("thelper "), "thelper", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect("gammadelta"), "tgd", cell_type_clean) ) |>
    mutate(cell_type_clean = if_else(cell_type_clean |> str_detect("natural killer"), "nk", cell_type_clean) ) |>
    
    
    mutate(cell_type_clean = cell_type_clean |> str_replace_all("  ", " ")) |>
    
    
    mutate(cell_type_clean = cell_type_clean |> str_replace("myeloid leukocyte", "myeloid")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("effector memory", "tem")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("effector", "tem")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace_all("cd8 t", "cd8")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("central memory", "tcm")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("gammadelta t", "gdt")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("nonclassical monocyte", "cd16 monocyte")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("classical monocyte", "cd14 monocyte")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("follicular b", "b")) |>
    mutate(cell_type_clean = cell_type_clean |> str_replace("unswitched memory", "memory")) |>
    
    mutate(cell_type_clean = cell_type_clean |> str_trim()) 
}

#' Clean and Standardize Cell Types
#'
#' This function takes a vector of cell types and applies a series of transformations
#' to clean and standardize them for better consistency.
#'
#' @importFrom stringr str_remove_all
#' @importFrom stringr str_trim
#'
#' @param .x A vector of cell types.
#'
#' @return A cleaned and standardized vector of cell types.
#'
#' @examples
#' cell_types <- c("CD4+ T-cells", "NK cells", "Blast-cells")
# cleaned_cell_types <- clean_cell_types(cell_types)

clean_cell_types = function(.x){
  .x |>
    str_remove_all("\\+") |>
    str_remove_all("cells") |>
    str_remove_all("cell") |>
    str_remove_all("blast") |>
    str_remove_all("-") |>
    str_trim()
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


get_manually_curated_immune_cell_types = function(){
  
  library(zellkonverter)
  library(Seurat)
  library(SingleCellExperiment) # load early to avoid masking dplyr::count()
  library(tidySingleCellExperiment)
  library(dplyr)
  library(cellxgenedp)
  library(tidyverse)
  #library(tidySingleCellExperiment)
  library(stringr)
  library(scMerge)
  library(glue)
  library(DelayedArray)
  library(HDF5Array)
  library(tidyseurat)
  library(celldex)
  library(SingleR)
  library(glmGamPoi)
  library(stringr)
  library(purrr)
  
  # # source("utility.R")
  # 
  # metadata_file = "/vast/projects/cellxgene_curated//metadata_0.2.rds"
  # file_curated_annotation_merged = "~/PostDoc/CuratedAtlasQueryR/dev/cell_type_curated_annotation_0.2.3.rds"
  # file_metadata_annotated = "/vast/projects/cellxgene_curated/metadata_annotated_0.2.3.rds"
  # annotation_directory = "/vast/projects/cellxgene_curated//annotated_data_0.2/"
  # 
  # # metadata_file = "/vast/projects/cellxgene_curated//metadata.rds"
  # # file_curated_annotation_merged = "~/PostDoc/CuratedAtlasQueryR/dev/cell_type_curated_annotation.rds"
  # # file_metadata_annotated = "/vast/projects/cellxgene_curated//metadata_annotated.rds"
  # # annotation_directory = "/vast/projects/cellxgene_curated//annotated_data_0.1/"
  # 
  # 
  # annotation_harmonised =
  #   dir(annotation_directory, full.names = TRUE) |>
  #   enframe(value="file") |>
  #   tidyr::extract(  file,".sample", "/([a-z0-9]+)\\.rds", remove = F) |>
  #   mutate(data = map(file, ~ .x |> readRDS() |> select(-contains("score")) )) |>
  #   unnest(data) |>
  #   
  #   # Format
  #   mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler, monaco_singler, ),	tolower	)) |>
  #   mutate(across(c(predicted.celltype.l1, predicted.celltype.l2, blueprint_singler, monaco_singler, ),	clean_cell_types	)) |>
  #   
  #   # Format
  # is_strong_evidence(predicted.celltype.l2, blueprint_singler) |> 
  #   
  # 
  # 
  # 
  # job::job({
  #   annotation_harmonised |>  saveRDS("~/PostDoc/CuratedAtlasQueryR/dev/annotated_data_0.2_temp_table.rds")
  # })
  # 

  annotation_harmonised = readRDS("~/PostDoc/CuratedAtlasQueryR/dev/annotated_data_0.2_temp_table.rds")
  
  # library(CuratedAtlasQueryR)
  metadata_df = readRDS(metadata_file)
  
  # Integrate with metadata
  
  annotation =
    metadata_df |>
    select(.cell, cell_type, file_id, .sample) |>
    as_tibble() |>
    left_join(read_csv("~/PostDoc/CuratedAtlasQueryR/dev/metadata_cell_type.csv"),  by = "cell_type") |>
    left_join(annotation_harmonised, by = c(".cell", ".sample")) |>
    
    # Clen cell types
    mutate(cell_type_clean = cell_type |> clean_cell_types())
  
  # annotation |>
  # 	filter(lineage_1=="immune") |>
  # 	count(cell_type, predicted.celltype.l2, blueprint_singler, strong_evidence) |>
  # 	arrange(!strong_evidence, desc(n)) |>
  # 	write_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm.csv")
  
  
  annotation_crated_confirmed =
    read_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm_manually_curated.csv") |>
    
    # TEMPORARY
    rename(cell_type_clean = cell_type) |>
    
    filter(!is.na(azhimut_confirmed) | !is.na(blueprint_confirmed)) |>
    filter(azhimut_confirmed + blueprint_confirmed > 0) |>
    
    # Format
    mutate(cell_type_harmonised = case_when(
      azhimut_confirmed ~ predicted.celltype.l2,
      blueprint_confirmed ~ blueprint_singler
    )) |>
    
    mutate(confidence_class = 1)
  
  
  
  # To avoid immune cell annotation if very contrasting evidence
  blueprint_definitely_non_immune = c(   "astrocytes" , "chondrocytes"  , "endothelial"  ,  "epithelial" ,  "fibros"  ,  "keratinocytes" ,    "melanocytes"  , "mesangial"  ,  "mv endothelial",   "myocytes" ,  "neurons"  ,  "pericytes" ,  "preadipocytes" , "skeletal muscle"  ,  "smooth muscle"      )
  
  
  
  annotation_crated_UNconfirmed =
    
    # Read
    read_csv("~/PostDoc/CuratedAtlasQueryR/dev/annotation_confirm_manually_curated.csv") |>
    
    # TEMPORARY
    rename(cell_type_clean = cell_type) |>
    
    filter(is.na(azhimut_confirmed) | (azhimut_confirmed + blueprint_confirmed) == 0) |>
    
    clean_cell_types_deeper() |> 
    
    mutate(cell_type_harmonised = "") |>
    
    # Classify strong evidence
    mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("cd8 cytokine secreting tem t") & blueprint_singler == "nk", T, blueprint_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("cd8 cytotoxic t") & blueprint_singler == "nk",  T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("cd8alphaalpha intraepithelial t") & predicted.celltype.l2 == "cd8 tem" & blueprint_singler == "cd8 tem", T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("mature t") & strong_evidence & predicted.celltype.l2  |> str_detect("tem|tcm"), T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean |> str_detect("myeloid") & strong_evidence & predicted.celltype.l2  == "cd16 mono", T, azhimut_confirmed) ) |>
    
    # Classify weak evidence
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("b", "B") & predicted.celltype.l2   == "b memory" & blueprint_singler == "classswitched memory b", T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("b", "B") & predicted.celltype.l2   %in% c("b memory", "b intermediate", "b naive", "plasma") & !blueprint_singler %in% c("classswitched memory b", "memory b", "naive b"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("b", "B") & !predicted.celltype.l2   %in% c("b memory", "b intermediate", "b naive") & blueprint_singler %in% c("classswitched memory b", "memory b", "naive b", "plasma"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "activated cd4" & predicted.celltype.l2  %in% c("cd4 tcm", "cd4 tem", "tregs"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "activated cd4" & blueprint_singler  %in% c("cd4 tcm", "cd4 tem", "tregs"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "activated cd8" & predicted.celltype.l2  %in% c("cd8 tcm", "cd8 tem"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "activated cd8" & blueprint_singler  %in% c("cd8 tcm", "cd8 tem"), T, blueprint_confirmed) ) |>
    
    # Monocyte macrophage
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 cd16 monocyte" & predicted.celltype.l2  %in% c("cd14 mono", "cd16 mono"), T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 cd16negative classical monocyte" & predicted.celltype.l2  %in% c("cd14 mono"), T, azhimut_confirmed) ) |>
    mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14 cd16negative classical monocyte" & blueprint_singler  %in% c("monocytes"), "cd14 mono", cell_type_harmonised) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14 monocyte" & predicted.celltype.l2  %in% c("cd14 mono"), T, azhimut_confirmed) ) |>
    mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14 monocyte" & blueprint_singler  %in% c("monocytes"), "cd14 mono", cell_type_harmonised) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd14low cd16 monocyte" & predicted.celltype.l2  %in% c("cd16 mono"), T, azhimut_confirmed) ) |>
    mutate(cell_type_harmonised = if_else(cell_type_clean == "cd14low cd16 monocyte" & blueprint_singler  %in% c("monocytes"), "cd16 mono", cell_type_harmonised) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd16 monocyte" & predicted.celltype.l2  %in% c("cd16 mono"), T, azhimut_confirmed) ) |>
    mutate(cell_type_harmonised = if_else(cell_type_clean == "cd16 monocyte" & blueprint_singler  %in% c("monocytes"), "cd16 mono", cell_type_harmonised) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "monocyte" & blueprint_singler  |> str_detect("monocyte|macrophage") & !predicted.celltype.l2 |> str_detect(" mono"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "monocyte" & predicted.celltype.l2 |> str_detect(" mono"), T, azhimut_confirmed) ) |>
    
    
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd4" & predicted.celltype.l2 |> str_detect("cd4|treg") & !blueprint_singler  |> str_detect("cd4"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "cd4" & !predicted.celltype.l2 |> str_detect("cd4") & blueprint_singler  |> str_detect("cd4|treg"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>
    
    
    mutate(azhimut_confirmed = if_else(cell_type_clean == "memory t" & predicted.celltype.l2 |> str_detect("tem|tcm") & !blueprint_singler  |> str_detect("tem|tcm"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "memory t" & !predicted.celltype.l2 |> str_detect("tem|tcm") & blueprint_singler  |> str_detect("tem|tcm"), T, blueprint_confirmed) ) |>
    
    
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8alphaalpha intraepithelial t" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8alphaalpha intraepithelial t" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>
    
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd8hymocyte" & predicted.celltype.l2 |> str_detect("cd8") & !blueprint_singler  |> str_detect("cd8"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "cd8hymocyte" & !predicted.celltype.l2 |> str_detect("cd8") & blueprint_singler  |> str_detect("cd8"), T, blueprint_confirmed) ) |>
    
    # B cells
    mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("memory b") & predicted.celltype.l2 =="b memory", T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("memory b") & blueprint_singler |> str_detect("memory b"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "immature b" & predicted.celltype.l2 =="b naive", T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "immature b" & blueprint_singler |> str_detect("naive b"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "mature b" & predicted.celltype.l2 %in% c("b memory", "b intermediate"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "mature b" & blueprint_singler |> str_detect("memory b"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "naive b" & predicted.celltype.l2 %in% c("b naive"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "naive b" & blueprint_singler |> str_detect("naive b"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "transitional stage b" & predicted.celltype.l2 %in% c("b intermediate"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "transitional stage b" & blueprint_singler |> str_detect("naive b") & !predicted.celltype.l2 %in% c("b intermediate"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "memory b" & predicted.celltype.l2 %in% c("b intermediate"), T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & predicted.celltype.l2 %in% c("b naive") & !blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & blueprint_singler |> str_detect("naive b") & predicted.celltype.l2 %in% c("hspc"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & predicted.celltype.l2 %in% c("hspc"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "precursor b", "prob") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
    
    # Plasma cells
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "plasma") & predicted.celltype.l2 == "plasma" , T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "plasma") & predicted.celltype.l2 == "plasma" , T, blueprint_confirmed) ) |>
    
    mutate(azhimut_confirmed = case_when(
      cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & predicted.celltype.l2 == "cd4 ctl" & blueprint_singler != "cd4 tcm" ~ T,
      cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & predicted.celltype.l2 == "cd4 tem" & blueprint_singler != "cd4 tcm" ~ T,
      TRUE ~ azhimut_confirmed
    ) ) |>
    mutate(blueprint_confirmed = case_when(
      cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & blueprint_singler == "cd4 tem" & predicted.celltype.l2 != "cd4 tcm" ~ T,
      cell_type_clean %in% c("cd4 cytotoxic t", "cd4 helper t") & blueprint_singler == "cd4 t" & predicted.celltype.l2 != "cd4 tcm" ~ T,
      TRUE ~ blueprint_confirmed
    ) ) |>
    
    mutate(azhimut_confirmed = if_else(cell_type_clean == "cd4hymocyte" & predicted.celltype.l2 |> str_detect("cd4|treg") & !blueprint_singler  |> str_detect("cd4"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "cd4hymocyte" & !predicted.celltype.l2 |> str_detect("cd4") & blueprint_singler  |> str_detect("cd4|treg"), T, blueprint_confirmed) ) |>
    
    mutate(azhimut_confirmed = case_when(
      cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 == "cd8 tem" & blueprint_singler != "cd8 tcm" ~ T,
      cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 == "cd8 tcm" & blueprint_singler != "cd8 tem" ~ T,
      TRUE ~ azhimut_confirmed
    ) ) |>
    mutate(blueprint_confirmed = case_when(
      cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 != "cd8 tem" & blueprint_singler == "cd8 tcm" ~ T,
      cell_type_clean %in% c("cd8 memory t") & predicted.celltype.l2 != "cd8 tcm" & blueprint_singler == "cd8 tem" ~ T,
      TRUE ~ blueprint_confirmed
    ) ) |>
    
    mutate(azhimut_confirmed = case_when(
      cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 == "cd4 tem" & blueprint_singler != "cd8 tcm" ~ T,
      cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 == "cd4 tcm" & blueprint_singler != "cd8 tem" ~ T,
      TRUE ~ azhimut_confirmed
    ) ) |>
    mutate(blueprint_confirmed = case_when(
      cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 != "cd4 tem" & blueprint_singler == "cd4 tcm" ~ T,
      cell_type_clean %in% c("cd4 memory t") & predicted.celltype.l2 != "cd4 tcm" & blueprint_singler == "cd4 tem" ~ T,
      TRUE ~ blueprint_confirmed
    ) ) |>
    
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "t") & blueprint_singler =="cd8 t" & predicted.celltype.l2 |> str_detect("cd8"), T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "t") & blueprint_singler =="cd4 t" & predicted.celltype.l2 |> str_detect("cd4|treg"), T, azhimut_confirmed) ) |>
    
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "treg") & blueprint_singler %in% c("tregs"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "treg") & predicted.celltype.l2 == "treg", T, azhimut_confirmed) ) |>
    
    
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tcm cd4") & blueprint_singler %in% c("cd4 tcm"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tcm cd4") & predicted.celltype.l2 == "cd4 tcm", T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tcm cd8") & blueprint_singler %in% c("cd8 tcm"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tcm cd8") & predicted.celltype.l2 == "cd8 tcm", T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tem cd4") & blueprint_singler %in% c("cd4 tem"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tem cd4") & predicted.celltype.l2 == "cd4 tem", T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "tem cd8") & blueprint_singler %in% c("cd8 tem"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tem cd8") & predicted.celltype.l2 == "cd8 tem", T, azhimut_confirmed) ) |>
    
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "tgd") & predicted.celltype.l2 == "gdt", T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "activated cd4") & predicted.celltype.l2 == "cd4 proliferating", T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "activated cd8") & predicted.celltype.l2 == "cd8 proliferating", T, azhimut_confirmed) ) |>
    
    
    
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("naive cd4", "naive t") & predicted.celltype.l2 %in% c("cd4 naive"), T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("naive cd8", "naive t") & predicted.celltype.l2 %in% c("cd8 naive"), T, azhimut_confirmed) ) |>
    
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("cd4 naive") & !blueprint_singler |> str_detect("clp|hcs|mpp|cd8"), T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("cd8 naive") & !blueprint_singler |> str_detect("clp|hcs|mpp|cd4"), T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "prot") & predicted.celltype.l2 %in% c("hspc"), T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "prot") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
    
    mutate(azhimut_confirmed = if_else(cell_type_clean == "dendritic" & predicted.celltype.l2 %in% c("asdc", "cdc2", "cdc1", "pdc"), T, azhimut_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean == "double negative t regulatory" & predicted.celltype.l2 == "dnt", T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "early t lineage precursor", "immature innate lymphoid") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "early t lineage precursor", "immature innate lymphoid") & predicted.celltype.l2 == "hspc" & blueprint_singler != "clp", T, azhimut_confirmed) ) |>
    
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("ilc1", "ilc2", "innate lymphoid") & blueprint_singler == "nk", T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c("ilc1", "ilc2", "innate lymphoid") & predicted.celltype.l2 %in% c( "nk", "ilc", "nk proliferating"), T, azhimut_confirmed) ) |>
    
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "immature t") & blueprint_singler %in% c("naive t"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "immature t") & predicted.celltype.l2 == "t naive", T, azhimut_confirmed) ) |>
    
    mutate(cell_type_harmonised = if_else(cell_type_clean == "fraction a prepro b", "naive b", cell_type_harmonised))  |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == "granulocyte" & blueprint_singler %in% c("eosinophils", "neutrophils"), T, blueprint_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c("immature neutrophil", "neutrophil") & blueprint_singler %in% c( "neutrophils"), T, blueprint_confirmed) ) |>
    
    mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("megakaryocyte") & blueprint_singler |> str_detect("megakaryocyte"), T, blueprint_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean |> str_detect("macrophage") & blueprint_singler |> str_detect("macrophage"), T, blueprint_confirmed) ) |>
    
    mutate(blueprint_confirmed = if_else(cell_type_clean %in% c( "nk") & blueprint_singler %in% c("nk"), T, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean %in% c( "nk") & predicted.celltype.l2 %in% c("nk", "nk proliferating", "nk_cd56bright", "ilc"), T, azhimut_confirmed) ) |>
    
    
    # If identical force
    mutate(azhimut_confirmed = if_else(cell_type_clean == predicted.celltype.l2 , T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean == blueprint_singler , T, blueprint_confirmed) ) |>
    
    # Perogenitor
    mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("progenitor|hematopoietic|precursor") & predicted.celltype.l2  == "hspc", T, azhimut_confirmed) ) |>
    mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("progenitor|hematopoietic|precursor") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, blueprint_confirmed) ) |>
    
    # Generic original annotation and stem for new annotations
    mutate(azhimut_confirmed = if_else(
      cell_type_clean  %in% c("T cell", "myeloid cell", "leukocyte", "myeloid leukocyte", "B cell") &
        predicted.celltype.l2  == "hspc" &
        blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), T, azhimut_confirmed) ) |>
    
    # Omit mature for stem
    mutate(blueprint_confirmed = if_else(cell_type_clean  |> str_detect("mature") & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), F, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean  |> str_detect("mature") & predicted.celltype.l2  == "hspc", F, azhimut_confirmed) ) |>
    
    # Omit megacariocyte for stem
    mutate(blueprint_confirmed = if_else(cell_type_clean  == "megakaryocyte" & blueprint_singler %in% c("clp","hcs", "mpp", "gmp"), F, blueprint_confirmed) ) |>
    mutate(azhimut_confirmed = if_else(cell_type_clean  == "megakaryocyte" & predicted.celltype.l2  == "hspc", F, azhimut_confirmed) ) |>
    
    # Mast cells
    mutate(cell_type_harmonised = if_else(cell_type_clean == "mast", "mast", cell_type_harmonised))  |>
    
    
    # Visualise
    #distinct(cell_type_clean, predicted.celltype.l2, blueprint_singler, strong_evidence, azhimut_confirmed, blueprint_confirmed) |>
    arrange(!strong_evidence, cell_type_clean) |>
    
    # set cell names
    mutate(cell_type_harmonised = case_when(
      cell_type_harmonised == "" & azhimut_confirmed ~ predicted.celltype.l2,
      cell_type_harmonised == "" & blueprint_confirmed ~ blueprint_singler,
      TRUE ~ cell_type_harmonised
    )) |>
    
    # Add NA
    mutate(cell_type_harmonised = case_when(cell_type_harmonised != "" ~ cell_type_harmonised)) |>
    
    # Add unannotated cells because datasets were too small
    mutate(cell_type_harmonised = case_when(
      is.na(cell_type_harmonised) & cell_type_clean  |> str_detect("progenitor|hematopoietic|stem|precursor") ~ "stem",
      
      is.na(cell_type_harmonised) & cell_type_clean == "cd14 monocyte" ~ "cd14 mono",
      is.na(cell_type_harmonised) & cell_type_clean == "cd16 monocyte" ~ "cd16 mono",
      is.na(cell_type_harmonised) & cell_type_clean %in% c("cd4 cytotoxic t", "tem cd4") ~ "cd4 tem",
      is.na(cell_type_harmonised) & cell_type_clean %in% c("cd8 cytotoxic t", "tem cd8") ~ "cd8 tem",
      is.na(cell_type_harmonised) & cell_type_clean |> str_detect("macrophage") ~ "macrophage",
      is.na(cell_type_harmonised) & cell_type_clean %in% c("mature b", "memory b", "transitional stage b") ~ "b memory",
      is.na(cell_type_harmonised) & cell_type_clean == "mucosal invariant t" ~ "mait",
      is.na(cell_type_harmonised) & cell_type_clean == "naive b" ~ "b naive",
      is.na(cell_type_harmonised) & cell_type_clean == "nk" ~ "nk",
      is.na(cell_type_harmonised) & cell_type_clean == "naive cd4" ~"cd4 naive",
      is.na(cell_type_harmonised) & cell_type_clean == "naive cd8" ~"cd8 naive",
      is.na(cell_type_harmonised) & cell_type_clean == "treg" ~ "treg",
      is.na(cell_type_harmonised) & cell_type_clean == "tgd" ~ "tgd",
      TRUE ~ cell_type_harmonised
    )) |>
    
    mutate(confidence_class = case_when(
      !is.na(cell_type_harmonised) & strong_evidence ~ 2,
      !is.na(cell_type_harmonised) & !strong_evidence ~ 3
    )) |>
    
    # Lowest grade annotation UNreliable
    mutate(cell_type_harmonised = case_when(
      
      # Get origincal annotation
      is.na(cell_type_harmonised) & cell_type_clean %in% c("neutrophil", "granulocyte") ~ cell_type_clean,
      is.na(cell_type_harmonised) & cell_type_clean %in% c("conventional dendritic", "dendritic") ~ "cdc",
      is.na(cell_type_harmonised) & cell_type_clean %in% c("classical monocyte") ~ "cd14 mono",
      
      # Get Seurat annotation
      is.na(cell_type_harmonised) & predicted.celltype.l2 != "eryth" & !is.na(predicted.celltype.l2) ~ predicted.celltype.l2,
      is.na(cell_type_harmonised) & !blueprint_singler %in% c(
        "astrocytes", "smooth muscle", "preadipocytes", "mesangial", "myocytes",
        "doublet", "melanocytes", "chondrocytes", "mv endothelial", "fibros",
        "neurons", "keratinocytes", "endothelial", "epithelial", "skeletal muscle", "pericytes", "erythrocytes", "adipocytes"
      ) & !is.na(blueprint_singler) ~ blueprint_singler,
      TRUE ~ cell_type_harmonised
      
    )) |>
    
    # Lowest grade annotation UNreliable
    mutate(cell_type_harmonised = case_when(
      
      # Get origincal annotation
      !cell_type_harmonised %in% c("doublet", "platelet") ~ cell_type_harmonised
      
    )) |>
    
    mutate(confidence_class = case_when(
      is.na(confidence_class) & !is.na(cell_type_harmonised) ~ 4,
      TRUE ~ confidence_class
    ))
  
  # Another passage
  
  # annotated_samples = annotation_crated_UNconfirmed |> filter(!is.na(cell_type_harmonised)) |>  distinct( cell_type, .sample, file_id)
  #
  # annotation_crated_UNconfirmed |>
  # 	filter(is.na(cell_type_harmonised))  |>
  # 	count(cell_type ,    cell_type_harmonised ,predicted.celltype.l2 ,blueprint_singler) |>
  # 	arrange(desc(n)) |>
  # 	print(n=99)
  
  
  annotation_all =
    annotation_crated_confirmed |>
    clean_cell_types_deeper() |> 
    bind_rows(
      annotation_crated_UNconfirmed
    ) |>
    
    # I have multiple confidence_class per combination of labels
    distinct() |>
    with_groups(c(cell_type_clean, predicted.celltype.l2, blueprint_singler), ~ .x |> arrange(confidence_class) |> slice(1)) |>
    
    # Simplify after harmonisation
    mutate(cell_type_harmonised =	case_when(
      cell_type_harmonised %in% c("b memory", "b intermediate", "classswitched memory b", "memory b" ) ~ "b memory",
      cell_type_harmonised %in% c("b naive", "naive b") ~ "b naive",
      cell_type_harmonised %in% c("nk_cd56bright", "nk", "nk proliferating", "ilc") ~ "ilc",
      cell_type_harmonised %in% c("mpp", "clp", "hspc", "mep", "cmp", "hsc", "gmp") ~ "stem",
      cell_type_harmonised %in% c("macrophages",  "macrophages m1", "macrophages m2") ~ "macrophage",
      cell_type_harmonised %in% c("treg",  "tregs") ~ "treg",
      cell_type_harmonised %in% c("gdt",  "tgd") ~ "tgd",
      cell_type_harmonised %in% c("cd8 proliferating",  "cd8 tem") ~ "cd8 tem",
      cell_type_harmonised %in% c("cd4 proliferating",  "cd4 tem") ~ "cd4 tem",
      cell_type_harmonised %in% c("eosinophils",  "neutrophils", "granulocyte", "neutrophil") ~ "granulocyte",
      cell_type_harmonised %in% c("cdc",  "cdc1", "cdc2", "dc") ~ "cdc",
      
      TRUE ~ cell_type_harmonised
    )) |>
    dplyr::select(cell_type_clean, cell_type_harmonised, predicted.celltype.l2, blueprint_singler, confidence_class) |>
    distinct()

  
  curated_annotation =
    annotation |>
    clean_cell_types_deeper() |> 
    filter(lineage_1=="immune") |>
    dplyr::select(
      .cell, .sample, cell_type, cell_type_clean, predicted.celltype.l2, blueprint_singler, monaco_singler) |>
    left_join(
      annotation_all ,
      by = c("cell_type_clean", "predicted.celltype.l2", "blueprint_singler")
    ) |>
    dplyr::select(
      .cell, .sample, cell_type, cell_type_harmonised, confidence_class,
      cell_annotation_azimuth_l2 = predicted.celltype.l2, cell_annotation_blueprint_singler = blueprint_singler,
      cell_annotation_monaco_singler = monaco_singler
    ) |>
    
    # Reannotation of generic cell types
    mutate(cell_type_harmonised = case_when(
      cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("effector memory") ~ "cd4 tem",
      cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("mait") ~ "mait",
      cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("central memory") ~ "cd4 tcm",
      cell_type_harmonised=="cd4 t" & cell_annotation_monaco_singler |> str_detect("naive") ~ "cd4 naive",
      cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("effector memory") ~ "cd8 tem",
      cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("central memory") ~ "cd8 tcm",
      cell_type_harmonised=="cd8 t" & cell_annotation_monaco_singler |> str_detect("naive") ~ "cd8 naive",
      cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler |> str_detect("non classical") ~ "cd16 mono",
      cell_type == "nonclassical monocyte" & cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="intermediate monocytes"   ~ "cd16 mono",
      cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler |> str_detect("^classical") ~ "cd14 mono",
      cell_type == "classical monocyte" & cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="intermediate monocytes"   ~ "cd14 mono",
      cell_type_harmonised=="monocytes" & cell_annotation_monaco_singler =="myeloid dendritic" & str_detect(cell_annotation_azimuth_l2, "cdc")   ~ "cdc",
      
      
      TRUE ~ cell_type_harmonised
    )) |>
    
    # Change CD4 classification for version 0.2.1
    mutate(confidence_class = if_else(
      cell_type_harmonised |> str_detect("cd4|mait|treg|tgd") & cell_annotation_monaco_singler %in% c("terminal effector cd4 t", "naive cd4 t", "th2", "th17", "t regulatory", "follicular helper t", "th1/th17", "th1", "nonvd2 gd t", "vd2 gd t"),
      3,
      confidence_class
    )) |>
    
    # Change CD4 classification for version 0.2.1
    mutate(cell_type_harmonised = if_else(
      cell_type_harmonised |> str_detect("cd4|mait|treg|tgd") & cell_annotation_monaco_singler %in% c("terminal effector cd4 t", "naive cd4 t", "th2", "th17", "t regulatory", "follicular helper t", "th1/th17", "th1", "nonvd2 gd t", "vd2 gd t"),
      cell_annotation_monaco_singler,
      cell_type_harmonised
    )) |>
    
    
    mutate(cell_type_harmonised = cell_type_harmonised |>
             str_replace("naive cd4 t", "cd4 naive") |>
             str_replace("th2", "cd4 th2") |>
             str_replace("^th17$", "cd4 th17") |>
             str_replace("t regulatory", "treg") |>
             str_replace("follicular helper t", "cd4 fh") |>
             str_replace("th1/th17", "cd4 th1/th17") |>
             str_replace("^th1$", "cd4 th1") |>
             str_replace("nonvd2 gd t", "tgd") |>
             str_replace("vd2 gd t", "tgd")
    ) |>
    
    # add immune_unclassified
    mutate(cell_type_harmonised = if_else(cell_type_harmonised == "monocytes", "immune_unclassified", cell_type_harmonised)) |>
    mutate(cell_type_harmonised = if_else(is.na(cell_type_harmonised), "immune_unclassified", cell_type_harmonised)) |>
    mutate(confidence_class = if_else(is.na(confidence_class), 5, confidence_class)) |>
    
    # drop uncommon cells
    mutate(cell_type_harmonised = if_else(cell_type_harmonised %in% c("cd4 t", "cd8 t", "asdc", "cd4 ctl"), "immune_unclassified", cell_type_harmonised))
  
  
  # Further rescue of unannotated cells, manually
  
  # curated_annotation |>
  # 	filter(cell_type_harmonised == "immune_unclassified") |>
  # 	count(cell_type   ,       cell_type_harmonised ,confidence_class ,cell_annotation_azimuth_l2 ,cell_annotation_blueprint_singler ,cell_annotation_monaco_singler) |>
  # 	arrange(desc(n)) |>
  # 	write_csv("curated_annotation_still_unannotated_0.2.csv")
  
  
  curated_annotation =
    curated_annotation |>
    left_join(
      read_csv("~/PostDoc/CuratedAtlasQueryR/dev/curated_annotation_still_unannotated_0.2_manually_labelled.csv") |>
        select(cell_type, cell_type_harmonised_manually_curated = cell_type_harmonised, confidence_class_manually_curated = confidence_class, everything()),
      by = join_by(cell_type, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler)
    ) |>
    mutate(
      confidence_class = if_else(cell_type_harmonised == "immune_unclassified", confidence_class_manually_curated, confidence_class),
      cell_type_harmonised = if_else(cell_type_harmonised == "immune_unclassified", cell_type_harmonised_manually_curated, cell_type_harmonised),
    ) |>
    select(-contains("manually_curated"), -n) |>
    
    # drop uncommon cells
    mutate(cell_type_harmonised = if_else(cell_type_harmonised %in% c("cd4 tcm", "cd4 tem"), "immune_unclassified", cell_type_harmonised))
  
  
  
  # # Recover confidence class == 4
  
  # curated_annotation |>
  # 	filter(confidence_class==4) |>
  # 	count(cell_type   ,       cell_type_harmonised ,confidence_class ,cell_annotation_azimuth_l2 ,cell_annotation_blueprint_singler ,cell_annotation_monaco_singler) |>
  # 	arrange(desc(n)) |>
  # 	write_csv("curated_annotation_still_unannotated_0.2_confidence_class_4.csv")
  
  curated_annotation =
    curated_annotation |>
    left_join(
      read_csv("~/PostDoc/CuratedAtlasQueryR/dev/curated_annotation_still_unannotated_0.2_confidence_class_4_manually_labelled.csv") |>
        select(confidence_class_manually_curated = confidence_class, everything()),
      by = join_by(cell_type, cell_type_harmonised, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler)
    ) |>
    mutate(
      confidence_class = if_else(confidence_class == 4 & !is.na(confidence_class_manually_curated), confidence_class_manually_curated, confidence_class)
    ) |>
    select(-contains("manually_curated"), -n)
  
  # Correct fishy stem cell labelling
  # If stem for the study's annotation and blueprint is non-immune it is probably wrong, 
  # even because the heart has too many progenitor/stem
  curated_annotation =
    curated_annotation |>
    mutate(confidence_class = case_when(
      cell_type_harmonised == "stem" & cell_annotation_blueprint_singler %in% c(
        "skeletal muscle", "adipocytes", "epithelial", "smooth muscle", "chondrocytes", "endothelial"
      ) ~ 5,
      TRUE ~ confidence_class
    ))
  
  
  curated_annotation_merged =
    
    # Fix cell ID
    metadata_df |>
    dplyr::select(.cell, .sample, cell_type) |>
    as_tibble() |>
    
    # Add cell type
    left_join(curated_annotation |> dplyr::select(-cell_type), by = c(".cell", ".sample")) |>
    
    # Add non immune
    mutate(cell_type_harmonised = if_else(is.na(cell_type_harmonised), "non_immune", cell_type_harmonised)) |>
    mutate(confidence_class = if_else(is.na(confidence_class) & cell_type_harmonised == "non_immune", 1, confidence_class)) |>
    
    # For some unknown reason
    distinct()
  
  
  curated_annotation_merged |>
    
    # Save
    saveRDS(file_curated_annotation_merged)
  
  metadata_annotated =
    curated_annotation_merged |>
    
    # merge with the rest of metadata
    left_join(
      metadata_df |>
        as_tibble(),
      by=c(".cell", ".sample", "cell_type")
    )
  
  # Replace `.` with `_` for all column names as it can create difficoulties for MySQL and Python
  colnames(metadata_annotated) = colnames(metadata_annotated) |> str_replace_all("\\.", "_")
  metadata_annotated = metadata_annotated |> rename(cell_ = `_cell`, sample_ = `_sample`)
  

  dictionary_connie_non_immune = 
    metadata_annotated |> 
    filter(cell_type_harmonised == "non_immune") |> 
    distinct(cell_type) |> 
    harmonise_names_non_immune() |> 
    rename(cell_type_harmonised_non_immune = cell_type_harmonised )
  
  metadata_annotated = 
    metadata_annotated |> 
    left_join(dictionary_connie_non_immune) |> 
    mutate(cell_type_harmonised = if_else(cell_type_harmonised=="non_immune", cell_type_harmonised_non_immune, cell_type_harmonised)) |> 
    select(-cell_type_harmonised_non_immune)
  
  
}
