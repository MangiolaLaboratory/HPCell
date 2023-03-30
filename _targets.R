# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint



# Set target options:
tar_option_set(
  packages = c("tidyr", "dplyr", "ggplot2",'purrr','Seurat','tidyseurat','glue','scater',
               'DropletUtils','EnsDb.Hsapiens.v86','here','stringr','readr','celldex'), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
# options(clustermq.scheduler = "slurm")
# options(clustermq.template = "clustermq.tmpl")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
tar_source('R/function__empty_droplet_identification.R')
tar_source('R/function_path.R')
# source("other_functions.R") # Source other scripts as needed. # nolint

# # Read arguments
# code_directory = '/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/third_party_software/target_pipeline/R'
# input_directory='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/input_directory'
# input_path = '/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/input_directory/JD1800154SL.rds'  # change here to file
# filtered = 'filtered'
# output_path_result = '/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory/preprocessing_results/empty_droplet_identification/JD1800154SL__empty_droplet_identification_output.rds'
# output_path_plot_pdf = '/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory/preprocessing_results/empty_droplet_identification/JD1800154SL__empty_droplet_identification_plot.pdf'
# output_path_plot_rds = '/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory/preprocessing_results/empty_droplet_identification/JD1800154SL__empty_droplet_identification_plot.rds'
# input_directory_demultiplexed = input_directory
# input_files_demultiplexed = dir(input_directory_demultiplexed, pattern = ".rds")
# input_path_demultiplexed = glue("{input_directory_demultiplexed}/{input_files_demultiplexed}")
# suffix = "__annotation_label_transfer"
# output_directory_label_transfer = glue("{result_directory}/preprocessing_results/annotation_label_transfer")
# output_paths_annotation_label_transfer =   glue("{output_directory_label_transfer}/{samples}{suffix}_output.rds")
# # Messy!!!!!!!! wrtie a script to change this situation!!!!!!!!

master_directory='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test'
result_directory='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory'
input_directory='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/third_party_software/target_pipeline/test_data'
code_directory='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/third_party_software/renv-optional'
metadata_path='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/metadata.rds'
reference_azimuth_path='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_azimuth.rds'
reports_directory='/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/cellsig/dev/xinpu_datascript/jascap_root_test/result_directory/preprocessing_results/reports'
filtered='filtered'


# Replace the target list below with your own:
list(
  tar_target(input_path,input_path_demultiplexed(input_directory)),
  tar_target(input_demultiplexed, input_path[[1]]),
  tar_target(samples,input_path[[2]]),
  tar_target(output_emptyDroplet_result,output_emptyDroplet(result_directory,samples,reports_directory)),
  # #tar_target(code_directory,'/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/xinpu/master_project/third_party_software/target_pipeline/R',format='file'),
  tar_target(
  emptyDroplet,
  run_empty_droplet_identification(code_directory,input_demultiplexed,filtered,output_emptyDroplet_result[[1]],output_emptyDroplet_result[[2]],output_emptyDroplet_result[[3]]),
  #tar_target()
    )
  # tar_target(
  #   annotation_label,annotation_label_transfer(code_directory,input_path_demultiplexed,empty_droplets,reference_azimuth_path,output_paths_annotation_label_transfer)
  # )
    
  )

# a<-input_path_demultiplexed(input_directory)
# print(output_emptyDroplet(result_directory,a[[2]],reports_directory))
#output_emptyDroplet(result_directory)

# delete later
#.  list(
#   tar_target(file, "data.csv", format = "file"),
#   tar_target(data, get_data(file)),
#   tar_target(model, fit_model(data)),
#   tar_target(plot, plot_model(model, data))
# )
