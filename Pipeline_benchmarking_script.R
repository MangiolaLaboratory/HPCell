install.packages("lobstr")
library(lobstr)
setwd("/vast/scratch/users/si.j")
initial_file_count <- 2
files <- list.files(pattern = "\\.txt$")
for (i in initial_file_count:length(files)) {
  # Memory usage before pipeline execution
  mem_before <- obj_size(get("preprocessed_seurat", envir = globalenv()))
  
  # Select the subset of files to process in this iteration
  file_subset <- files[1:i] 
  max_workers <- 100  
  workers_per_sample <- 4
  number_of_samples <- length(file_subset)
  total_workers <- min(number_of_samples * workers_per_sample, max_workers)
  
  # Initialize computing resources for all files, assuming this is intended once per loop
  computing_resources = crew_controller_slurm(
    name = "my_controller",
    slurm_memory_gigabytes_per_cpu = 20,
    slurm_cpus_per_task = 1, 
    workers = total_workers, 
    verbose = FALSE
  )
  
  # Time and run your pipeline function
  time_taken <- system.time({
    preprocessed_seurat <- run_targets_pipeline(
      input_data = file_subset,
      tissue = "pbmc",
      computing_resources = computing_resources,
      sample_column = "sampleName", 
      store = store,  # Ensure 'store' is defined outside the loop or appropriately initialized
      input_reference = NULL, 
      cell_type_annotation_column = "cellAnno"
    )
  })
  
  # Memory usage after pipeline execution
  mem_after <- obj_size(get("preprocessed_seurat", envir = globalenv()))
  mem_used_this_run <- mem_after - mem_before
  
  # Output the time and memory used for this run
  cat("Sample size:", length(file_subset), "\n",
      "Time taken: User time =", time_taken["user.self"], 
      "System time =", time_taken["sys.self"], 
      "Elapsed time =", time_taken["elapsed"], "seconds\n",
      "Memory used:", format(mem_used_this_run, units = "Mb"), "\n\n")
}
