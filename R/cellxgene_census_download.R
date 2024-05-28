#' Download `[SingleCellExperiment::SingleCellExperiment-class]` from cellxgene-census
#' 
#' @description
#' This function downloads Census data, renames Ensembl features to gene names, and
#' save the data to a defined path.
#' 
#' @param df A data frame of metadadta from Census.
#' @param save_path A character vector of path to save Census SingleCellExperiment objects.
#' @param experiment A character vector of experiment name. Either sce for `SingleCellExperiment-class`,
#'   or seurat for `Seurat-class`.
#' @param data_output_type A character vector of data type in the output. The data type 
#'   need to be one of the following: rds, hdf5, or anndata.
#' @importFrom purrr pmap
#' @importFrom glue glue
#' @importFrom dplyr mutate 
#' @importFrom zellkonverter writeH5AD
#' @importFrom cellxgene.census get_single_cell_experiment get_seurat
#' @noRd
get_census_data <- function(df, 
                            save_path,
                            experiment,
                            data_output_type) {
  # add a file name column for file saving
  df <- df %>%
    mutate(filename = paste(gsub("[- ']", "_", assay), gsub("[- ]", "_", disease), gsub("[- ]", "_", donor_id),
                            gsub("[- ]", "_", sex), gsub("[- ]", "_", self_reported_ethnicity), 
                            gsub("[- ]", "_", tissue), gsub("[- ]", "_", development_stage), sep = "___")) 
  
  # map each row in the data frame
  results <- pmap(df, function(...) {
    args <- list(...)
    # construct the filter string following Census syntax
    obs_value_filter <- 
      glue(
        "donor_id == '{args$donor_id}' & 
        assay == '{gsub(\"'\", \"\\\\'\", args$assay)}' & 
        disease == '{args$disease}' & sex == '{args$sex}' & 
        self_reported_ethnicity == '{args$self_reported_ethnicity}' & 
        tissue == '{args$tissue}' & 
        development_stage == '{args$development_stage}'"
      )
    
    print(glue("Processing: {args$filename}"))
    if (experiment == "sce") {
      data <-
        cellxgene.census::get_single_cell_experiment(census = census,
                                                     organism = "Homo sapiens",
                                                     obs_value_filter = obs_value_filter) |>
        rename_features()
    } else if (experiment == "seurat") {
      data <- cellxgene.census::get_seurat(census = census,
                                           organism = "Homo sapiens",
                                           obs_value_filter = obs_value_filter) |>
        rename_features()
    }
    
    if (data_output_type == "rds") {
      file_path <- file.path(save_path, glue("{args$filename}.rds"))
      saveRDS(data, file = file_path)
      
    } else if (data_output_type == "hdf5") {
      file_path <- file.path(save_path, glue("{args$filename}"))
      saveHDF5SummarizedExperiment(data, dir = file_path, verbose = TRUE)
      
    } else if (data_output_type == "anndata") {
      file_path <- file.path(save_path, glue("{args$filename}.h5ad"))
      zellkonverter::writeH5AD(sce_object, file = file_path)
      
    }
  })
}

#' Rename features from ensembl to gene names
#' @param data A `[Seurat::Seurat-class]` or `[SingleCellExperiment::SingleCellExperiment-class]` object.
#' @param assay An optional character vector of lenth one. If provided, it should indicate 
#'   the assay to be used for the analysis.
#' @return A `[Seurat::Seurat-class]` or `[SingleCellExperiment::SingleCellExperiment-class]` object
#'   with renamed genes
#' @importFrom Seurat GetAssayData GetAssay
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom tibble tibble
#' @noRd
rename_features <- function(data,
                            assay = NULL) {
  if (inherits(data, "Seurat")) {
    if (is.null(assay))
      assay <- data@assays |> names()
    
    meta_features <-
      GetAssay(data, assay = assay) |> slot("meta.features")
    counts <- GetAssayData(data, assay = assay, layer = "counts")
    assay_data <- GetAssayData(data, assay = assay, layer = "data")
    
    gene_names_map <- tibble(
      ensembl = meta_features |> rownames(),
      # replace _ to - in features
      gene_names = gsub("_", "-", meta_features$feature_name)
    )
    # ensure features order are the same
    if (identical(counts |> rownames(),
                  gene_names_map$ensembl)) {
      rownames(counts) <- gene_names_map$gene_names
    }
    if (identical(assay_data |> rownames(),
                  gene_names_map$ensembl)) {
      rownames(assay_data) <- gene_names_map$gene_names
    }
    
    # modify meta_features slot
    meta_features$feature_name <- NULL
    rownames(meta_features) <- gene_names_map$gene_names
    meta_features <- meta_features
    
    data
    
  } else if (inherits(data, "SingleCellExperiment")) {
    current_features <- rowData(data)$feature_name
    
    gene_names_map <- tibble(ensembl = rownames(data),
                             gene_names = gsub("_", "-", current_features))
    
    if (identical(rownames(data), gene_names_map$ensembl)) {
      rownames(data) <- gene_names_map$gene_names
      rownames(data@assays@data$counts) <-
        gene_names_map$gene_names
    }
    
    rowData(data)$feature_name <- NULL
    rowData(data)$gene_name <- gene_names_map$gene_names
    
    data
  }
}
