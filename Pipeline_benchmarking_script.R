install.packages("lobstr")
library(lobstr)

# Defining resources
Cores <- c(3)
Sample_size <- c(1,2, 5, 10, 20, 50, 100, 137)
#Sample_size <- c(2, 5, 10, 20, 50)
#Sample_size <- c(2, 5, 10, 20, 50, 100, 138)



setwd("/vast/scratch/users/si.j/susan_fibrosis")
#setwd("/stornext/General/scratch/GP_Transfer/susan_fibrosis")
#initial_file_count <- 2
files <- list.files()
length(files)
store <- "/stornext/General/scratch/GP_Transfer/si.j/benchmark_store"

#store_contents <- list.files(store)
#need_invalidate <- length(store_contents) > 0
# for (i in initial_file_count:length(files)) {
#   
#   tar_invalidate(names = everything(), store = store)
#   
#   # Memory usage before pipeline execution
#   if(exists("preprocessed_seurat", envir = globalenv())) {
#     mem_before <- obj_size(get("preprocessed_seurat", envir = globalenv()))
#   } else {
#     mem_before <- 0
#   }
for(core in Cores) {
  for(sample_size in Sample_size) {
    if(length(files) < sample_size) {
      break # Break if the sample_size exceeds the available files
    }
    #setwd("~/HPCell")
  tar_invalidate(names = everything(), store = store)
  
  # Select the subset of files to process in this iteration
  #file_subset <- files[1:i] 
  file_subset <- files[1:sample_size]
  max_workers <- 100  
  workers_per_sample <- 4
  number_of_samples <- length(file_subset)
  #total_workers <- min(number_of_samples * workers_per_sample, max_workers)
  total_workers <- min(core * sample_size, length(file_subset))
  # Initialize computing resources for all files
  computing_resources = crew_controller_slurm(
    name = "my_controller",
    slurm_memory_gigabytes_per_cpu = 20,
    slurm_cpus_per_task = 1, 
    workers = total_workers, 
    verbose = FALSE, 
    seconds_idle = 30
  )
    # Time and run your pipeline function
    #setwd("~/HPCell")
    time_taken <- system.time({
      preprocessed_seurat <- run_targets_pipeline(
        input_data = file_subset,
        tissue = "pbmc",
        computing_resources = computing_resources,
        sample_column = "sampleName", 
        store = store,  
        input_reference = NULL, 
        cell_type_annotation_column = "cellAnno"
      )
    })
    
    # Memory usage after pipeline execution
    #mem_after <- obj_size(get("preprocessed_seurat", envir = globalenv()))
    #mem_used_this_run <- mem_after - mem_before
    
    # Output the time and memory used for this run
    cat("Running with", core, "cores for", sample_size, "samples, using", total_workers, "workers\n")
    cat("Sample size:", length(file_subset), "\n",
        "Time taken: User time =", time_taken["user.self"], 
        "System time =", time_taken["sys.self"], 
        "Elapsed time =", time_taken["elapsed"], "seconds\n")
    #"Memory used:", format(mem_used_this_run, units = "Mb"), "\n\n")
  }
}


###PLOTTING

library(ggplot2)

# Data
data <- data.frame(
  SampleSize = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
  UserTime = c(23.944, 26.71, 36.797, 43.466, 46.846, 54.226, 59.966, 69.946, 76.142, 86.875),
  SystemTime = c(5.856, 6.216, 8.143, 9.805, 11.088, 13.302, 15.419, 18.439, 21.315, 23.56),
  ElapsedTime = c(278.592, 280.922, 310.607, 306.815, 312.007, 333.49, 333.517, 384.848, 380.266, 396.131)
)

# Melting data for ggplot
data_melted <- reshape2::melt(data, id.vars = "SampleSize")

# Plotting
ggplot(data_melted, aes(x = SampleSize, y = value, colour = variable)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Sample Size", y = "Time (seconds)", title = "Performance Metrics by Sample Size", color = "Metric") +
  scale_colour_manual(values = c("UserTime" = "cornflowerblue", "SystemTime" = "slategrey", "ElapsedTime" = "coral"))



### Rewriting alternative script 
# Defining resources
Cores <- c(16)
Sample_size <- c(1, 2, 5, 10, 20, 50)

setwd("/stornext/General/scratch/GP_Transfer/susan_fibrosis/")
files <- list.files()
store <- "/stornext/General/scratch/GP_Transfer/si.j/store_pipeline_benchmark_fibrosis_all_data_3"

# Initialize results dataframe
results <- data.frame(SampleNumber = integer(), DataSize = numeric(), RunningTimeMin = numeric(), Cores = integer())

for(core in Cores) {
  for(sample_size in Sample_size) {
    if(length(files) < sample_size) {
      break # Break if the sample_size exceeds the available files
    }
    
    #tar_invalidate(names = everything(), store = store)
    
    # Select the subset of files to process in this iteration
    file_subset <- files[1:sample_size]
    total_workers <- min(core, length(file_subset))
    
    # Initialize computing resources for all files
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
        store = store,  
        input_reference = NULL, 
        cell_type_annotation_column = "cellAnno"
      )
    })
    
    
    # Output the results
    results <- rbind(results, data.frame(
      SampleNumber = sample_size,
      DataSize = data_size,
      RunningTimeMin = time_taken["elapsed"] / 60, # Convert seconds to minutes
      Cores = core
    ))
    
    cat("Running with", core, "cores for", sample_size, "samples, using", total_workers, "workers\n")
    cat("Sample size:", length(file_subset), "\n",
        "Time taken: User time =", time_taken["user.self"], 
        "System time =", time_taken["sys.self"], 
        "Elapsed time =", time_taken["elapsed"], "seconds\n")
  }
}

# Save the results to a CSV file
write.csv(results, "benchmark_results.csv", row.names = FALSE)

### PLOTTING

library(ggplot2)
library(reshape2)

# Melting data for ggplot
data_melted <- melt(results, id.vars = "SampleNumber")

# Plotting
ggplot(data_melted, aes(x = SampleNumber, y = value, colour = variable)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Sample Number", y = "Value", title = "Performance Metrics by Sample Number", color = "Metric") +
  scale_colour_manual(values = c("DataSize" = "cornflowerblue", "RunningTimeMin" = "slategrey", "Cores" = "coral", "TotalMemoryGB" = "purple"))





















