
# Read arguments
args = commandArgs(trailingOnly=TRUE)
code_directory = args[[1]]
reference_label_fine = args[[2]]
input_path_preprocessing_output = args[3:(length(args)-2)]
output_path_sample_cell_type = args[[length(args)-1]]
output_path_sample = args[[length(args)]]

renv::load(project = code_directory)

library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(glue)
library(tidyseurat)
library(tidySingleCellExperiment)
library(tidysc)
library(tidybulk)
library(rlang)
library(stringr)
library(purrr)

# # Parallelise internally
# library(furrr)
# options(future.globals.maxSize = 49 * 1024 ^ 3) # 50Gb
# plan(strategy = "multisession", workers = 10)

# Do we have RNA and also ADT
assays = readRDS(input_path_preprocessing_output[1])@assays |> names() |> intersect(c("RNA", "ADT"))

# Create dir
output_path_sample |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

pseudobulk =
  input_path_preprocessing_output |>

  # Aggregate
  map(~ {
    library(rlang)
    readRDS(.x) |>
      aggregate_cells(c(sample, !!as.symbol(reference_label_fine)), slot = "counts", assays=assays) |>

      # Reshape to make RNA and ADT both features
      pivot_longer(
        cols = assays,
        names_to = "data_source",
        values_to = "count"
      ) |>
      filter(!count |> is.na()) |>

      # Some manipulation to get unique feature because RNA and ADT
      # both can have sma name genes
      rename(symbol = feature) |>
      mutate(data_source = str_remove(data_source, "abundance_")) |>
      unite( ".feature", c(symbol, data_source), remove = FALSE) |>

      # Covert
      as_SummarizedExperiment(
        .sample = c( sample,!!as.symbol(reference_label_fine)),
        .transcript = .feature,
        .abundance = count
      )

  })

# pseudobulk |> saveRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/PPCG_tumour_microenvironment/PPCG_deconvolution_signatures_single_cell/PPCG_deconvolution_signatures_single_cell_PROCESSED_v3/preprocessing_results/pseudobulk_preprocessing/temp_pseudobulk_preprocessing_sample_output.rds")

# This should not be needed if I create count files weel form the beginning
# Select only common column
common_columns =
  pseudobulk |>
  map(~ .x |> tidySummarizedExperiment::as_tibble() |> colnames()) |>
  unlist() |>
  table() %>%
  .[.==max(.)] |>
  names()

all_genes =
  pseudobulk |>
  map(~ .x |> rownames()) |>
  unlist() |>
  unique() |>
  as.character()

# Select and save
pseudobulk |>

  # Add missing genes
  map(~{

    missing_genes = all_genes |> setdiff(rownames(.x))

    missing_matrix = matrix(rep(0, length(missing_genes) * ncol(.x)), ncol = ncol(.x))

    rownames(missing_matrix) = missing_genes
    colnames(missing_matrix) = colnames(.x)

    new_se = SummarizedExperiment(assay = list(count = missing_matrix))
    colData(new_se) = colData(.x)
    #rowData(new_se) =  DataFrame(symbol = missing_genes, row.names = missing_genes)
    rowData(.x) = NULL
    .x = .x |> rbind(new_se)

    .x[all_genes,]

  }) |>

  map(~ .x |> tidySummarizedExperiment::select(any_of(common_columns)))   %>%

  do.call(cbind, .) |>

  # Save
  saveRDS(output_path_sample_cell_type)

gc()

# ONLY SAMPLE
pseudobulk =
  input_path_preprocessing_output |>

  # Aggregate
  map(~
                   readRDS(.x) |>
                   aggregate_cells(c(sample), slot = "counts", assays=assays) |>

            # Reshape to make RNA and ADT both features
            pivot_longer(
              cols = assays,
              names_to = "data_source",
              values_to = "count"
            ) |>
            filter(!count |> is.na()) |>

            # Some manipulation to get unique feature because RNa and ADT both can have sma name genes
            rename(symbol = feature) |>
            mutate(data_source = str_remove(data_source, "abundance_")) |>
            unite( ".feature", c(symbol, data_source), remove = FALSE) |>

            # Covert
            as_SummarizedExperiment(
              .sample = c( sample),
              .transcript = .feature,
              .abundance = count
            )
  )

# This should not be needed if I create count files weel form the beginning
# Select only common column
common_columns =
  pseudobulk |>
  map(~ .x |> tidySummarizedExperiment::as_tibble() |> colnames()) |>
  unlist() |>
  table() %>%
  .[.==max(.)] |>
  names()

# Select and save
pseudobulk |>
  # Add missing genes
  map(~{

    missing_genes = all_genes |> setdiff(rownames(.x))

    missing_matrix = matrix(rep(0, length(missing_genes) * ncol(.x)), ncol = ncol(.x))

    rownames(missing_matrix) = missing_genes
    colnames(missing_matrix) = colnames(.x)

    new_se = SummarizedExperiment(assay = list(count = missing_matrix))
    colData(new_se) = colData(.x)
    #rowData(new_se) =  DataFrame(symbol = missing_genes, row.names = missing_genes)
    rowData(.x) = NULL
    .x = .x |> rbind(new_se)

    .x[all_genes,]

  }) |>

  map(~ .x |> tidySummarizedExperiment::select(any_of(common_columns)))   %>%

  do.call(cbind, .) |>

  # Save
  saveRDS(output_path_sample)

