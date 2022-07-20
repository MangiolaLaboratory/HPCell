# Create makeflow file for doublet identification
# ~/third_party_sofware/cctools-7.2.0-x86_64-centos7/bin/makeflow -T slurm -J 200 analysis/communication_fast/communication_fast.makeflow 

# This script slipts the dataset and creates the list of files in a specific directory
library(tidyverse)
library(Seurat)
library(tidyseurat)
library(glue)
library(here)

# Read arguments
master_input_directory = commandArgs(trailingOnly=TRUE)[[1]]

tab = "\t"
code_directory = "analysis/communication"
cores = 4

# Preprocessing
input_directory_preprocessing_output = glue("{master_input_directory}/preprocessing_results/preprocessing_output")
input_files_preprocessing_output = dir(input_directory_preprocessing_output, pattern = "output.rds") 
input_path_preprocessing_output = glue("{input_directory_preprocessing_output}/{input_files_preprocessing_output}")

# Output
output_directory = glue("{master_input_directory}/fast_pipeline_results/communication")
output_path_result =   glue("{output_directory}/communication_output.rds")

# Output
samples = input_files_preprocessing_output |> str_remove("__preprocessing_output_output.rds")
output_directory = glue("{master_input_directory}/fast_pipeline_results/communication")
output_path_result =   glue("{output_directory}/{samples}_ligand_receptor_count_output.rds")


# Create input
commands = 
  c(
  "CATEGORY=create_output\nMEMORY=10024\nCORES=4\nWALL_TIME=9000",
  glue("{output_path_result}:{input_path_preprocessing_output}\n{tab}Rscript {code_directory}/run__ligand_receptor_count.R {input_path_preprocessing_output} {output_path_result}")
  
) 

# CellChat results
output_path_plot_overall =  glue("{output_directory}/plot_communication_overall.pdf")
output_path_plot_heatmap =  glue("{output_directory}/plot_communication_heatmap.pdf")
output_path_plot_circle =  glue("{output_directory}/plot_communication_circle.pdf")
output_path_values_communication =  glue("{output_directory}/values_communication.rds")


commands |> 
  c(
    "CATEGORY=create_output\nMEMORY=50024\nCORES=2\nWALL_TIME=18000",
    glue("{output_path_plot_overall} {output_path_plot_heatmap} {output_path_plot_circle} {output_path_values_communication}:{paste(output_path_result, collapse=\" \")}\n{tab}Rscript {code_directory}/run__differential_communication.R {paste(output_path_result, collapse=\" \")} {master_input_directory}/metadata.rds {output_path_plot_overall} {output_path_plot_heatmap} {output_path_plot_circle} {output_path_values_communication}")
    
  ) %>% 
  
  write_lines(glue("{master_input_directory}/pipeline.makeflow"), append = TRUE)