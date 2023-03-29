library(glue)

input_path_demultiplexed<-function(input_directory){
  input_directory_demultiplexed = input_directory
  input_files_demultiplexed = dir(input_directory_demultiplexed, pattern = ".rds")
  input_path_demultiplexed = glue("{input_directory_demultiplexed}/{input_files_demultiplexed}")
}


# Empty droplet
output_emptyDroplet<-function(result_directory){
  suffix = "__empty_droplet_identification"
  output_directory_empty_droplets = glue("{result_directory}/preprocessing_results/empty_droplet_identification")
  output_path_empty_droplets =   glue("{output_directory_empty_droplets}/{samples}{suffix}_output.rds")
  output_path_plot_pdf =  glue("{output_directory_empty_droplets}/{samples}{suffix}_plot.pdf")
  output_path_plot_rds =  glue("{output_directory_empty_droplets}/{samples}{suffix}_plot.rds")
  report_empty_droplets = glue("{reports_directory}/report{suffix}.md")
  return(list(output_path_empty_droplets,output_path_plot_pdf,output_path_plot_rds,report_empty_droplets))
}

# output_path_plot_pdf
# 
# output_path_plot_rds