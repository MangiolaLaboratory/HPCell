
set.seed(42)

number_features_overall = 300
number_features_per_cell_type = 300


library(dplyr); library(tidyr); library(ggplot2)
library(Seurat)
library(tidyseurat)
library(glue)
library(purrr)
library(magrittr)
library(tibble)
library(jascap)

variable_gene_identification<-function(code_directory,tissue, input_paths,cell_type_column_for_subsetting,output_path){
 
    # Divide into 5 chunks
    input_paths_chunks = input_paths |> split( rep_len(1:5, length(input_paths)) |> sort()) 
    
    # Create dir
    output_path |> dirname() |> dir.create( showWarnings = FALSE)

    # Filter data
    counts =
      tibble(
        input_data_paths = input_paths_chunks[[1]],
        input_empty_droplets_paths = input_paths_chunks[[2]],
        input_dead_cells_paths = input_paths_chunks[[3]],
        input_doublets_paths = input_paths_chunks[[4]],
        input_cell_type_annotation_paths = input_paths_chunks[[5]]
      ) |>
      mutate(seurat = pmap(
        list(input_data_paths, input_empty_droplets_paths, input_dead_cells_paths, input_doublets_paths, input_cell_type_annotation_paths),
        ~ readRDS(..1) |>

          # Filter empty droplets
          left_join(readRDS(..2) |> select(.cell, empty_droplet), by = ".cell") |>
          filter(!empty_droplet) |>

          # Filter dead cells
          left_join(readRDS(..3) |> select(.cell, alive), by = ".cell") |>
          filter(alive) |>

          # Filter doublets
          left_join(readRDS(..4) |> select(.cell, scDblFinder.class), by = ".cell") |>
          filter(scDblFinder.class=="singlet") |>

          # Filter Red blood cells and platelets using pbmc seurat for any tissue
          left_join(readRDS(..5), by = ".cell") |>

          # If PBMC filter for ertthrocytes
          when(
            tissue |>
              tolower() |>
              equals("pbmc") ~
              filter(., !predicted.celltype.l2 %in% c("Eryth", "Platelet")),
            ~ (.)
          )


      )) |>
       pull(seurat) |>
       purrr::reduce(merge)

    all_features_df =
      counts |>
      GetAssay("RNA") |>   #origianl version Assays("RNA")
      rownames() |>
      enframe(value = "feature") |>
      select(-name) |>
      mutate(group= "all_features")

# I BROKE cell_type_column_for_subsetting = "none"
  seurat_to_variable_features(
    counts,
    "RNA",
    sample,
    sym(cell_type_column_for_subsetting),
    features_number_independent_of_cell_groups = 300,
    features_number_per_cell_group = 300
  )|>
  bind_rows(all_features_df) |>
  saveRDS(output_path)

}


#实验！！！！
#input_files<-input_files_variable_gene_identification(input_demultiplexed,emptyDroplet,alive,doublet,annotation)
#variable_gene_identification(code_directory,'solid',input_files,reference_label_coarse,output_variable)
#Rscript function__variable_gene_identification.R aa bb input1 input2 input3 input4 inutp5 dd ee
#Rscript function__variable_gene_identification.R aa bb /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/third_party_software/target_pipeline/test_data/JD1800154SL.rds /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory/preprocessing_results/empty_droplet_identification/JD1800154SL__empty_droplet_identification_output.rds /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory/preprocessing_results/alive_identification/JD1800154SL__alive_identification_output.rds /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory/preprocessing_results/doublet_identification/JD1800154SL__doublet_identification_output.rds /stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory/preprocessing_results/annotation_label_transfer/JD1800154SL__annotation_label_transfer_output.rds dd ee





