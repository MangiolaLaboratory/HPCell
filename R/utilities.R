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
#' @export
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
