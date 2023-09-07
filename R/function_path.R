library(glue)
library(stringr)

#check modality
reference_label_fine<-function(tissue){
  reference_fine<- tissue |> when(
    (.) == "pbmc" ~ "monaco_first.labels.fine",
    (.) =="solid" ~ "blueprint_first.labels.fine",
    (.) == "atypical" ~ "none"
  )
}

reference_label_coarse<-function(tissue){
  reference_label_coarse<- tissue |> when(
    (.) == "pbmc" ~ "monaco_first.labels.coarse",
    (.) == "solid" ~ "blueprint_first.labels.coarse",
    (.) == "atypical" ~ "none"
  )
}

#input_path_demultiplexed
#return the input path
input_path_demultiplexed<-function(input_directory){ #input directory is the data file directory
  input_directory_demultiplexed = input_directory
  input_files_demultiplexed = dir(input_directory_demultiplexed, pattern = ".rds")
  input_path_demultiplexed = glue("{input_directory_demultiplexed}/{input_files_demultiplexed}")
  samples = input_files_demultiplexed |> str_remove(".rds")
  return(list(input_path_demultiplexed,samples))
}

#input_path_demultiplexed
#return the sample name from the demultiplexed file path
sample_name<-function(input_path_demultiplexed){ # input_path_demultiplexed is the path of each rds file 
  input_files_demultiplexed = str_extract(input_path_demultiplexed,'[^"/"]*\\.rds')
  samples = input_files_demultiplexed |> str_remove(".rds")
  return(samples)
}


# Empty droplet
output_emptyDroplet<-function(result_directory,samples,reports_directory){
  suffix = "__empty_droplet_identification"
  output_directory_empty_droplets = glue("{result_directory}/preprocessing_results/empty_droplet_identification")
  output_path_empty_droplets =   glue("{output_directory_empty_droplets}/{samples}{suffix}_output.rds")
  output_path_plot_pdf =  glue("{output_directory_empty_droplets}/{samples}{suffix}_plot.pdf")
  output_path_plot_rds =  glue("{output_directory_empty_droplets}/{samples}{suffix}_plot.rds")
  report_empty_droplets = glue("{reports_directory}/report{suffix}.md")
  return(list(output_path_empty_droplets,output_path_plot_pdf,output_path_plot_rds,report_empty_droplets))
}

# annotation label transfer
output_annotation_label_transfer<-function(result_directory,samples){
  suffix = "__annotation_label_transfer"
  output_directory_label_transfer = glue("{result_directory}/preprocessing_results/annotation_label_transfer")
  output_paths_annotation_label_transfer =glue("{output_directory_label_transfer}/{samples}{suffix}_output.rds")
}

# alive identification
output_alive_identification<-function(result_directory,samples){
  suffix = "__alive_identification"
  output_directory = glue("{result_directory}/preprocessing_results/alive_identification")
  output_path_alive =   glue("{output_directory}/{samples}{suffix}_output.rds")
}

# doublet identification
output_doublet_identification<-function(result_directory,samples){
  suffix = "__doublet_identification"
  output_directory_doublet = glue("{result_directory}/preprocessing_results/doublet_identification")
  output_path_doublet_identification =  glue("{output_directory_doublet}/{samples}{suffix}_output.rds")
}

# variable gene identification
input_files_variable_gene_identification<-function(input_path_demultiplexed,output_path_empty_droplets,
                                                   output_path_alive, output_path_doublet_identification, 
                                                   output_paths_annotation_label_transfer){
  input_files=c(input_path_demultiplexed,output_path_empty_droplets,
                      output_path_alive, output_path_doublet_identification, 
                      output_paths_annotation_label_transfer)
}

output_variable_gene_identification<-function(result_directory,metadata_path,samples){
  suffix = "__variable_gene_identification"
  output_directory = glue("{result_directory}/preprocessing_results/variable_gene_identification")
  #add batch later
  metadata= readRDS(metadata_path)
  batch=samples
  output_file = glue("{output_directory}/{batch}{suffix}_output.rds")
}


#cell cycle scoring
output_cell_cycle_scoring<-function(result_directory,samples){
  suffix = "__cell_cycle_scoring"
  output_directory_cell_cycle_scoring = glue("{result_directory}/preprocessing_results/cell_cycle_scoring")
  output_cell_cycle_scoring =   glue("{output_directory_cell_cycle_scoring}/{samples}{suffix}_output.rds")
}

#normalisation, non batch variation removal
output_non_batch_variation_removal<-function(result_directory,samples){
  suffix = "__non_batch_variation_removal"
  output_directory_non_batch_variation_removal = glue("{result_directory}/preprocessing_results/non_batch_variation_removal")
  output_path_non_batch_variation_removal =   glue("{output_directory_non_batch_variation_removal}/{samples}{suffix}_output.rds")
}


