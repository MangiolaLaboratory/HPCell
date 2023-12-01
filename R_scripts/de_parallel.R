library(tidyverse)
library(HPCell)
library(targets)
library(glue)
library(tidySummarizedExperiment)
library(tidybulk)

# Get input

# my_data = 
#   tidySummarizedExperiment::se |> 
#   tidybulk::keep_abundant()
# my_data_df = 
#   tibble(name = "my_data", se = list(my_data)) 
# my_store = tempfile(tmpdir = ".")
# 
# 
# my_data_df |> 
#   hpcell_test_differential_abundance(store = my_store)

# Small dataset
se =
  tidySummarizedExperiment::se |>
  tidybulk::keep_abundant()
  
se |>
  hpcell_test_differential_abundance(
    ~ dex + (1 | cell), 
    computing_resources = crew.cluster::crew_controller_slurm(
      name = "slurm",
      slurm_memory_gigabytes_per_cpu = 5,
      slurm_cpus_per_task = 2,
      workers = 200,
      verbose = T
    )
  )


se |>
  hpcell_test_differential_abundance(   ~ dex + (1 | cell)    )
  
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


se_big = 
  readRDS("R_scripts/pseudobulk_big.rds") |> 
  nest(data = -cell_type_harmonised) |> 
  mutate(data = map(data, ~ .x |> tidybulk::as_SummarizedExperiment(.sample, .feature, RNA))) |> 
  mutate(data = map(
    data,
    ~ {
      
      if(ncol(.x) == 0) return(.x)
     
      # Filter
      se = 
        .x |> 
        
        # Eliminate complete confounders
        samples_NOT_complete_confounders_for_sex_assay() |>
        samples_NOT_complete_confounders_for_sex_disease()
      
      rm(.x)
      gc()
      
      # Filter disease
      se =
        se |>
        filter(disease %in% (
          se |>
            distinct(disease, ethnicity_simplified) |>
            count(disease) |>
            filter(n>1) |>
            pull(disease)
        ))
      
      # Return prematurely
      if(ncol(se) == 0) return(se)
      if(se |> distinct(sex, ethnicity_simplified) |> count(sex) |> pull(n) |> max() == 1) return(se)
      
      # Vell types with enough samples
      tissues_to_keep =
        se |>
        distinct(sample_, tissue_harmonised) |>
        count(  tissue_harmonised) |>
        filter(n > 3) |>
        pull(tissue_harmonised)
      
      se =
        se |>
        
        # Scale contninuous variables
        mutate(age_days = scale(age_days) |> as.numeric()) |>
        
        # Filter cell types to keep
        filter(tissue_harmonised %in% tissues_to_keep) |>
        
        # otherwise I get error for some reason
        mutate(across(any_of(c("sex", "ethnicity_simplified", "assay_simplified", "file_id", "tissue_harmonised")), as.character)) |>
        mutate(ethnicity_simplified = ethnicity_simplified |> str_replace("European", "aaa_European")) 
      
      # Drop random effect grouping with no enough data
      combinations_to_keep = 
        se |>
        distinct(sex, ethnicity_simplified, tissue_harmonised) |>
        add_count(tissue_harmonised) |>
        filter(n>1)
      
        se |>
        right_join(combinations_to_keep)
      
    }, 
    .progress = TRUE
  )) |> 
  filter(map_int(data, ncol) > 0) |> 
  mutate(formula = map(
    data,
    ~ {
      
      se = .x
      # Build the formula
      factors = 
        c("age_days", "sex", "ethnicity_simplified", "assay_simplified", ".aggregated_cells", "disease") |>
        enframe(value = "factor") |>
        mutate(n = map_int(
          factor, ~ se |> select(.x) |> distinct() |> nrow()
        )) |>
        filter(n>1) |>
        pull(factor) |>
        str_c(collapse = " + ")
      
      random_effects =
        c("age_days", "sex", "ethnicity_simplified") |>
        enframe(value = "factor") |>
        mutate(n = map_int(
          factor, ~ se |> select(all_of(.x)) |> distinct() |> nrow()
        ))   |>
        filter(n>1) |>
        pull(factor) |>
        str_c(collapse = " + ")
      
      # The default
      my_formula = glue("~ {factors}")
      
      if( 
        se |> distinct(tissue_harmonised) |> nrow() > 1 &
        length(random_effects) > 0
      ) 
        my_formula = glue("{my_formula} + (1 + {random_effects} | tissue_harmonised)")
      
      if( 	se |> distinct(file_id) |> nrow() > 1	)
        my_formula = glue("{my_formula} + (1 | file_id)")
      
      # Add the interaction
      if(se |> 
         nest_detect_complete_confounder(age_days, sex) |> 
         filter(n1 + n2 <= 2) |> 
         nrow() == 0
      ) my_formula = my_formula |> str_replace_all("age_days \\+ sex", "age_days * sex")
      
      as.formula(my_formula)
    }
  ))

se_big |> 
  pull(data) %>%
  .[[1]] |> 
  tidybulk::identify_abundant(factor_of_interest = ethnicity_simplified) |> 
  tidybulk::test_differential_abundance(
    ~ age_days * sex + ethnicity_simplified + assay_simplified + .aggregated_cells + (1 | file_id), 
    method = "glmmSeq_lme4", cores = 4
  )

se_big |>
  slice(1) |> 
  hpcell_map_test_differential_abundance(
    ~ age_days * sex + ethnicity_simplified + assay_simplified + .aggregated_cells + (1 | file_id) ,
    data,
    cell_type_harmonised,
    computing_resources = crew.cluster::crew_controller_slurm(
      name = "slurm",
      slurm_memory_gigabytes_per_cpu = 5,
      slurm_cpus_per_task = 2,
      workers = 200,
      verbose = T
    )
  )

