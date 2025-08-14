library(tidyverse)
library(HPCell)
library(targets)
library(glue)
library(tidySummarizedExperiment)
library(tidybulk)
library(tictoc)
library(microbenchmark)

# Small dataset
se =
  tidySummarizedExperiment::se |>
  tidybulk::keep_abundant()
  
se |>
  test_differential_abundance_hpc(
    ~ dex + (1 | cell), 
    computing_resources = crew.cluster::crew_controller_slurm(
      name = "slurm",
      slurm_memory_gigabytes_per_cpu = 5,
      slurm_cpus_per_task = 2,
      workers = 200,
      verbose = T, 
      seconds_idle = 30
    )
  )


se |>
  test_differential_abundance_hpc(   ~ dex + (1 | cell)    )
  
# Big dataset

# Ethnicity
samples_NOT_complete_confounders_for_sex_assay = function(se){
  
  
  
  se = 
    se |> 
    # distinct(assay_simplified, sex, .sample) |>
    # 
    nest(se_data = -c(assay_simplified, sex)) |>
    
    # How many ethnicity per assay
    nest(data = -assay_simplified) |> 
    mutate(n1 = map_int(data, ~ .x |> distinct(sex) |> nrow())) |> 
    unnest(data) |> 
    
    # How many assay per ethnicity
    nest(data = - sex) |> 
    mutate(n2 = map_int(data, ~ .x |> distinct(assay_simplified) |> nrow())) |> 
    unnest(data) 
  
  # Replace ethnicity
  dummy_assay = se |> arrange(desc(n1 + n2)) |> dplyr::slice(1) |> pull(assay_simplified)
  
  se |>
    mutate(assay_simplified = if_else(n1 + n2 < 3, dummy_assay, assay_simplified)) 	|>
    
    # # Filter
    # filter(!(n1==1 & n2==1)) |>
    select(-n1, -n2) |>
    
    unnest_summarized_experiment(se_data) 
  # |>
  # 	pull(.sample) |>
  # 	unique()
}

samples_NOT_complete_confounders_for_sex_disease = function(se){
  
  
  
  se = 
    se |> 
    #distinct(disease, sex, .sample) |>
    
    nest(se_data = -c(disease, sex)) |>
    
    # How many ethnicity per assay
    nest(data = -disease) |> 
    mutate(n1 = map_int(data, ~ .x |> distinct(sex) |> nrow())) |> 
    unnest(data) |> 
    
    # How many assay per ethnicity
    nest(data = - sex) |> 
    mutate(n2 = map_int(data, ~ .x |> distinct(disease) |> nrow())) |> 
    unnest(data) 
  
  
  # Replace ethnicity
  dummy_ethnicity = se |> arrange(desc(n1 + n2)) |> dplyr::slice(1) |> pull(sex)
  
  se |>
    mutate(sex = if_else(n1 + n2 < 3, dummy_ethnicity, sex)) 	|>
    
    # # Filter
    # filter(!(n1==1 & n2==1)) |>
    select(-n1, -n2) |>
    
    unnest_summarized_experiment(se_data)
  # |>
  # 	pull(.sample) |>
  # 	unique()
}

samples_NOT_complete_confounders_for_assay_ethnicity = function(se){
  se |> 
    #distinct(ethnicity_simplified, assay_simplified, .sample) |>
    
    nest(se_data = -c(ethnicity_simplified, assay_simplified)) |>
    
    # How many ethnicity per assay
    nest(data = -ethnicity_simplified) |> 
    mutate(n1 = map_int(data, ~ .x |> distinct(assay_simplified) |> nrow())) |> 
    unnest(data) |> 
    
    # How many assay per ethnicity
    nest(data = - assay_simplified) |> 
    mutate(n2 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |> 
    unnest(data) |>
    
    filter(n1+n2>2) |>
    select(-n1, -n2) |>
    unnest_summarized_experiment(se_data)
}

samples_NOT_complete_confounders_for_disease_ethnicity = function(se){
  se |> 
    #distinct(ethnicity_simplified, assay_simplified, .sample) |>
    
    nest(se_data = -c(ethnicity_simplified, disease)) |>
    
    # How many ethnicity per assay
    nest(data = -ethnicity_simplified) |> 
    mutate(n1 = map_int(data, ~ .x |> distinct(disease) |> nrow())) |> 
    unnest(data) |> 
    
    # How many assay per ethnicity
    nest(data = - disease) |> 
    mutate(n2 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |> 
    unnest(data) |>
    
    filter(n1+n2>2) |>
    select(-n1, -n2) |>
    unnest_summarized_experiment(se_data)
}

samples_NOT_complete_confounders_for_age_ethnicity = function(se){
  
  clean = 
    se |> 
    #distinct(ethnicity_simplified, assay_simplified, .sample) |>
    
    nest(se_data = -c(ethnicity_simplified, age_days)) |>
    
    # How many ethnicity per assay
    nest(data = -ethnicity_simplified) |> 
    mutate(n1 = map_int(data, ~ .x |> distinct(age_days) |> nrow())) |> 
    unnest(data) |> 
    
    # How many assay per ethnicity
    nest(data = - age_days) |> 
    mutate(n2 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |> 
    unnest(data) |>
    
    filter(n1+n2>2) |>
    select(-n1, -n2) 
  if(
    nrow(clean) == 0 || clean |>
    distinct(ethnicity_simplified) |>
    nrow() == 1
  )
    clean = 
      se |> 
      #distinct(ethnicity_simplified, assay_simplified, .sample) |>
      mutate(age_days = age_days > 1) |>
      nest(se_data = -c(ethnicity_simplified, age_days)) |>
      
      # How many ethnicity per assay
      nest(data = -ethnicity_simplified) |> 
      mutate(n1 = map_int(data, ~ .x |> distinct(age_days) |> nrow())) |> 
      unnest(data) |> 
      
      # How many assay per ethnicity
      nest(data = - age_days) |> 
      mutate(n2 = map_int(data, ~ .x |> distinct(ethnicity_simplified) |> nrow())) |> 
      unnest(data) |>
      
      filter(n1+n2>2) |>
      select(-n1, -n2) 
  
  clean |> unnest_summarized_experiment(se_data)
}

nest_detect_complete_confounder = function(.data, .col1, .col2){
  
  .col1 = enquo(.col1)
  .col2 = enquo(.col2)
  
  .data |>
    
    nest(se_data = -c(!!.col1, !!.col2)) |>
    
    # How many ethnicity per assay
    nest(data = -!!.col1) |>
    mutate(n1 = map_int(data, ~ .x |> distinct(!!.col2) |> nrow())) |>
    unnest(data) |>
    
    # How many assay per ethnicity
    nest(data = - !!.col2) |>
    mutate(n2 = map_int(data, ~ .x |> distinct(!!.col1) |> nrow())) |>
    unnest(data)
}


# se_big = 
#   readRDS("R_scripts/pseudobulk_big.rds") |> 
#   nest(data = -cell_type_harmonised) |> 
#   mutate(data = map(data, ~ .x |> tidybulk::as_SummarizedExperiment(.sample, .feature, RNA))) |> 
#   mutate(data = map(
#     data,
#     ~ {
#       
#       if(ncol(.x) == 0) return(.x)
#      
#       # Filter
#       se = 
#         .x |> 
#         
#         # Eliminate complete confounders
#         samples_NOT_complete_confounders_for_sex_assay() |>
#         samples_NOT_complete_confounders_for_sex_disease()
#       
#       rm(.x)
#       gc()
#       
#       # Filter disease
#       se =
#         se |>
#         filter(disease %in% (
#           se |>
#             distinct(disease, ethnicity_simplified) |>
#             count(disease) |>
#             filter(n>1) |>
#             pull(disease)
#         ))
#       
#       # Return prematurely
#       if(ncol(se) == 0) return(se)
#       if(se |> distinct(sex, ethnicity_simplified) |> count(sex) |> pull(n) |> max() == 1) return(se)
#       
#       # Vell types with enough samples
#       tissues_to_keep =
#         se |>
#         distinct(sample_, tissue_harmonised) |>
#         count(  tissue_harmonised) |>
#         filter(n > 3) |>
#         pull(tissue_harmonised)
#       
#       se =
#         se |>
#         
#         # Scale contninuous variables
#         mutate(age_days = scale(age_days) |> as.numeric()) |>
#         
#         # Filter cell types to keep
#         filter(tissue_harmonised %in% tissues_to_keep) |>
#         
#         # otherwise I get error for some reason
#         mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised")), as.character)) |>
#         mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European")) 
#       
#       # Drop random effect grouping with no enough data
#       combinations_to_keep = 
#         se |>
#         distinct(sex, ethnicity_simplified, tissue_harmonised) |>
#         add_count(tissue_harmonised) |>
#         filter(n>1)
#       
#         se |>
#         right_join(combinations_to_keep)
#       
#     }, 
#     .progress = TRUE
#   )) |> 
#   filter(map_int(data, ncol) > 0) |> 
#   mutate(formula = map(
#     data,
#     ~ {
#       
#       se = .x
#       # Build the formula
#       factors = 
#         c("age_days", "sex", "ethnicity_simplified", "assay_simplified", ".aggregated_cells", "disease") |>
#         enframe(value = "factor") |>
#         mutate(n = map_int(
#           factor, ~ se |> select(.x) |> distinct() |> nrow()
#         )) |>
#         filter(n>1) |>
#         pull(factor) |>
#         str_c(collapse = " + ")
#       
#       random_effects =
#         c("age_days", "sex", "ethnicity_simplified") |>
#         enframe(value = "factor") |>
#         mutate(n = map_int(
#           factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
#         ))   |>
#         filter(n>1) |>
#         pull(factor) |>
#         str_c(collapse = " + ")
#       
#       # The default
#       my_formula = glue("~ {factors}")
#       
#       if( 
#         se |> distinct(tissue_harmonised) |> nrow() > 1 &
#         length(random_effects) > 0
#       ) 
#         my_formula = glue("{my_formula} + (1 + {random_effects} | tissue_harmonised)")
#       
#       if( 	se |> distinct(file_id) |> nrow() > 1	)
#         my_formula = glue("{my_formula} + (1 | file_id)")
#       
#       # Add the interaction
#       if(se |> 
#          nest_detect_complete_confounder(age_days, sex) |> 
#          filter(n1 + n2 <= 2) |> 
#          nrow() == 0
#       ) my_formula = my_formula |> str_replace_all("age_days \\+ sex", "age_days * sex")
#       
#       as.formula(my_formula)
#     }
#   )) |> 
# mutate(data = map(data, tidybulk::identify_abundant, factor_of_interest = ethnicity_simplified )) |> 
# slice(1:22)
# 
# se_big |> saveRDS("R_scripts/se_big.rds", compress = "xz")

se_big = readRDS("~/PostDoc/HPCell/R_scripts/se_big.rds")

tic()
se_big |> 
  pull(data) %>%
  .[[24]] |> 
  tidybulk::identify_abundant(factor_of_interest = ethnicity_simplified) |> 
  tidybulk::test_differential_abundance(
    ~ age_days * sex + ethnicity_simplified + assay_simplified + .aggregated_cells + (1 | file_id), 
    method = "glmmSeq_lme4", cores = 2
  )
time_parallel_local = toc()


slurm = crew.cluster::crew_controller_slurm(
  name = "slurm",
  slurm_memory_gigabytes_per_cpu = 5,
  slurm_cpus_per_task = 1,
  workers = 200,
  verbose = T, 
  seconds_idle = 30
)


microbenchmark(
    se_big |>
    mutate(data = map2(
      data,
      formula ,
      ~ tidybulk::test_differential_abundance(
        .x, .y, method = "glmmSeq_lme4", cores = 1
      )
    )), 
    se_big |>
      mutate(data = map2_test_differential_abundance_hpc(
        data,
        formula ,
        computing_resources = slurm
      )),
  times = 1
)




se =
  tidySummarizedExperiment::se |>
  tidybulk::keep_abundant()

  se |>
    test_differential_abundance_hpc(
      ~ dex + (1 | cell),
      computing_resources = slurm
    )


  
  # Test of multilevel models
  
  # Load necessary libraries
  library(CuratedAtlasQueryR)  # For accessing curated atlas data
  library(tidyverse)          # For data manipulation and visualization
  
  counts = 
  # Retrieve and process the metadata
  get_metadata() |> 
    # Filter for specific conditions
    filter(
      tissue_harmonised == "blood",          # Select only 'blood' tissue
      cell_type_harmonised == "b memory",    # Select only 'b memory' cell type
      disease %in% c("normal", "COVID-19")   # Filter for 'normal' and 'COVID-19' diseases
    ) |> 
    # Convert the data frame to a tibble for better handling
    as_tibble() |> 
    # Nest data excluding sample and disease columns
    nest(data_cells = -c(sample_, disease)) |> 
    # Add a new column 'n' that contains the row count of each nested dataframe
    mutate(n = map_int(data_cells, nrow)) |> 
    # Filter out groups with less than 10 rows
    filter(n > 9) |> 
    # Nest the data again, this time excluding the disease column
    nest(data_samples = -disease) |> 
    # Add columns for lower and upper count thresholds
    mutate(count_low = c(30, 30), count_high = c(500, 30)) |> 
    # Apply function to each row using parallel mapping
    mutate(data_samples = pmap(
      list(data_samples, count_low, count_high),
      ~ bind_rows(
        # Select nine samples closest to the lower count threshold
        ..1 |> 
          arrange(abs(n-..2)) |> 
          head(9),
        
        # Select one sample closest to the higher count threshold
        ..1 |> 
          arrange(abs(n-..3)) |> 
          head(1)
      ) |> 
        distinct()
    )) |> 
    # Unnest the nested 'data_samples' dataframe
    unnest(data_samples) |> 
    # Finally, unnest the 'data_cells' to expand the nested data
    unnest(data_cells) |> 
    get_single_cell_experiment()
  
  assay(counts) = assay(counts) |> as.matrix()
  
  
 de_results = 
   counts |> 
   mutate(disease = if_else(disease == "normal", "a_normal", "b_COVID-19")) |> 
   tidybulk::keep_abundant(factor_of_interest = disease) |> 
    tidybulk::test_differential_abundance(
      ~disease, 
      scaling_method = "TMMwsp"
    ) 
 
 de_results |> 
   as("SummarizedExperiment") |> 
   tidybulk::scale_abundance(method = "TMMwsp") |> 
   _[rowData(de_results)$FDR<0.05,] |> 
   _[7, , drop=FALSE] |> 
   # _[
   #   de_results |> 
   #     rowData() |> 
   #     as_tibble(rownames = "feature") |> 
   #     arrange(desc(abs(logFC))) |> 
   #     head(10) |> 
   #     pull(feature)
   #  ,] |> 
   
   ggplot(aes(sample_, counts_scaled)) +
   geom_boxplot(aes(fill = disease), varwidth = TRUE, outlier.shape = NA) + 
   geom_jitter(shape = ".") +
   facet_wrap(~.feature) +
   scale_y_log10()
 
 de_results_lme4 = 
   counts |> 
   mutate(disease = if_else(disease == "normal", "a_normal", "b_COVID-19")) |> 
   tidybulk::keep_abundant(factor_of_interest = disease) |> 
   HPCell::test_differential_abundance_hpc(
     ~disease + (1|sample_), 
     scaling_method = "TMMwsp", 
     computing_resources = 
       crew.cluster::crew_controller_slurm(
         name = "slurm",
         slurm_memory_gigabytes_per_cpu = 5,
         slurm_cpus_per_task = 1,
         workers = 200,
         verbose = T, 
         seconds_idle = 30
       )
   ) 
 
   as("SummarizedExperiment") |> 
   as_tibble()
    filter(FDR<0.05) |> 
    nest(data = .feature)
    tidybulk::pivot_transcript(feature)
  
  counts = 
    tibble(
    sample = letters[1:10],
    condition  = rep(c("a_untreated", "b_treated"), each = 5), 
    number_of_cells = 10,
    mrna_abundance = rep(c(100, 10), each =5)
  ) |> 
    mutate(
      number_of_cells = if_else(sample==max(sample), 500, number_of_cells),
      mrna_abundance = if_else(sample==max(sample), 500, mrna_abundance),
    ) |> 
    mutate(counts = map2(
      mrna_abundance,
      number_of_cells,
      ~ rnbinom(mu = .x,n = .y, size = 5) |> 
        enframe(value = "counts", name = "cell")
    )) |> 
    unnest(counts) |> 
    unite("cell_name", cell, sample, remove = FALSE)  
  
  counts |> 
    ggplot(aes(fct_reorder(sample, condition), counts)) + 
    geom_boxplot(aes(fill = condition), varwidth = TRUE, outlier.shape = NA) + 
    geom_jitter(shape = ".") +
    scale_y_log10()
  
  
  # single cell with fixed effect models
  counts |> 
    mutate(log_counts = log1p(counts)) |> 
    lm(log_counts~condition, data = _) |> 
    summary()
  
  counts |> 
    mutate(log_counts = log1p(counts)) |> 
    with_groups(c(sample, condition, number_of_cells), ~ .x |> summarise(log_counts = mean(log_counts))) |> 
    lm(log_counts~condition + number_of_cells, data = _, ) |> 
    summary()
  
  counts |> 
    mutate(log_counts = log1p(counts)) |> 
    lmerTest::lmer(log_counts~condition + (1|sample), data = _) |> 
    summary()
    
  counts |> 
    mutate(feature = "gene_x") |> 
    filter(counts>0) |> 
    tidybulk::test_differential_abundance(
      ~condition, 
      cell_name, feature, counts, 
      scaling_method = "none" 
    ) |> 
    pivot_transcript(feature)
  
  
    lmerTest::lmer(log_counts~condition + (1|sample), data = _) |> 
    summary()