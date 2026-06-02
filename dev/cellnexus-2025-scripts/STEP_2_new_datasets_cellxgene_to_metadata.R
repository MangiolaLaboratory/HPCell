# Step2
.rs.restartR()
library(tidyverse)
library(targets)
library(glue)
library(arrow)
library(duckdb)

result_directory = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026"

tar_script({
  
  #-----------------------#
  # Input
  #-----------------------#
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(qs)
  library(crew)
  library(crew.cluster)
  
  #-----------------------#
  # Packages
  #-----------------------#
  tar_option_set(
    packages = c(
      "zellkonverter", "cellxgenedp", "CuratedAtlasQueryR", "stringr", "tibble", "tidySingleCellExperiment", "dplyr", "Matrix",
      "glue", "qs",  "purrr", "tidybulk", "tidySummarizedExperiment",  "crew", "magrittr", "digest", "readr", "forcats"
    ),
    
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    #  debug = "dataset_id_sce_b5312463451d7ee3", 
    cue = tar_cue(mode = "never"),
    format = "qs",
    
    
    #-----------------------#
    # SLURM
    #-----------------------#
    
    controller = crew_controller_group(
      
      
      crew_controller_slurm(
        name = "slurm_1_5",
        slurm_memory_gigabytes_per_cpu = 5,
        slurm_cpus_per_task = 1,
        workers = 200,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_10",
        slurm_memory_gigabytes_per_cpu = 10,
        slurm_cpus_per_task = 1,
        workers = 100,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_20",
        slurm_memory_gigabytes_per_cpu = 20,
        slurm_cpus_per_task = 1,
        workers = 100,
        tasks_max = 5,
        verbose = T, , 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_40",
        slurm_memory_gigabytes_per_cpu = 40,
        slurm_cpus_per_task = 1,
        workers = 50,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_80",
        slurm_memory_gigabytes_per_cpu = 80,
        slurm_cpus_per_task = 1,
        workers = 30,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      ),
      crew_controller_slurm(
        name = "slurm_1_200",
        slurm_memory_gigabytes_per_cpu = 200,
        slurm_cpus_per_task = 1,
        workers = 5,
        tasks_max = 5,
        verbose = T, 
        seconds_idle = 30
      )
    ),
    resources = tar_resources(crew = tar_resources_crew("slurm_1_10"))  
    
  )
  
  sample_heuristics = function(col_data){
    
    col_data |>
      
      # Sort sample ID
      # Fix some sample id missing
      when(unique(.$dataset_id) %in% c(
        "11b86bc3-6d4d-4e28-903a-0361ea8f6bdf",
        "492b0613-ff5b-4fca-a585-503fc4102e4f",
        "11b86bc3-6d4d-4e28-903a-0361ea8f6bdf",
        "0e8f9ce4-46e5-434e-9ca0-e769d1dd27ea"
      ) ~ mutate(., PatientID = glue("{sample} {replicate} {time_point} {target}") |> as.character()) , ~ (.)) %>%
      when(unique(.$dataset_id) %in% c(
        "0273924c-0387-4f44-98c5-2292dbaab11e",
        "a16bec18-5c9f-40ad-8169-12c5199c7506",
        "556bb449-bbef-43d3-9487-87031fc0decb"
      ) ~ mutate(., PatientID = glue("{Collection.ID} {Genotype} {Location}")|> as.character()) , ~ (.)) %>%
      when(unique(.$dataset_id) %in% c(
        "b83afdc1-baa1-42c0-bd5b-cb607084757d"
      ) ~ mutate(., PatientID = glue("{sex} {development_stage} {disease}")|> as.character()) , ~ (.)) %>%
      when(unique(.$dataset_id) %in% c(
        "3fe53a40-38ff-4f25-b33b-e4d60f2289ef",
        "5c1cc788-2645-45fb-b1d9-2f43d368bba8"
      ) ~ mutate(., PatientID = glue("{Batch} {Fetus_id} {Development_day} {sex} {tissue} {disease}")|> as.character()) , ~ (.)) |>
      
      
      mutate_if(is.factor, as.character) |>
      type_convert(guess_integer = TRUE) |>
      mutate_if(is.integer, as.character) |>
      
      # Convert types
      when("donor_id" %in% colnames(.) ~ mutate(., donor_id = donor_id |> as.character() ), ~(.)) |>
      when("Cluster" %in% colnames(.) ~ mutate(., Cluster = Cluster |> as.character() ), ~(.)) |>
      when("cluster_id" %in% colnames(.) ~ mutate(., cluster_id = cluster_id |> as.character() ), ~(.)) |>
      when("Batch" %in% colnames(.) ~ mutate(., Batch = Batch |> as.character() ), ~(.)) |>
      when("batch" %in% colnames(.) ~ mutate(., batch = batch |> as.character() ), ~(.)) |>
      when("age" %in% colnames(.) ~ mutate(., age = age |> as.character() ), ~(.)) |>
      when("BMI" %in% colnames(.) ~ mutate(., BMI = BMI |> as.character() ), ~(.)) |>
      when("donor_BMI" %in% colnames(.) ~ mutate(., donor_BMI = donor_BMI |> as.character() ), ~(.)) |>
      when("author_cell_type" %in% colnames(.) ~ mutate(., author_cell_type = author_cell_type |> as.character() ), ~(.)) |>
      when("time_point" %in% colnames(.) ~ mutate(., time_point = time_point |> as.character() ), ~(.)) |>
      when("cluster" %in% colnames(.) ~ mutate(., cluster = cluster |> as.character() ), ~(.)) |>
      when("ClusterID" %in% colnames(.) ~ mutate(., ClusterID = ClusterID |> as.character() ), ~(.)) |>
      when("Stage" %in% colnames(.) ~ mutate(., Stage = Stage |> as.character() ), ~(.)) |>
      when("individual" %in% colnames(.) ~ mutate(., individual = individual |> as.character() ), ~(.)) |>
      when("recurrent_cluster" %in% colnames(.) ~ mutate(., recurrent_cluster = recurrent_cluster |> as.character() ), ~(.)) |>
      when("PatientID" %in% colnames(.) ~ mutate(., PatientID = PatientID |> as.character() ), ~(.)) |>
      when("PMI" %in% colnames(.) ~ mutate(., PMI = PMI |> as.character() ), ~(.)) |>
      when("n_genes" %in% colnames(.) ~ mutate(., n_genes = n_genes |> as.numeric() ), ~(.)) |>
      when("n_counts" %in% colnames(.) ~ mutate(., n_counts = n_counts |> as.numeric() ), ~(.)) |>
      when("n_genes_by_counts" %in% colnames(.) ~ mutate(., n_genes_by_counts = n_genes_by_counts |> as.numeric() ), ~(.)) |>
      when("nUMI" %in% colnames(.) ~ mutate(., nUMI = nUMI |> as.numeric() ), ~(.)) |>
      when("percent.cortex" %in% colnames(.) ~ mutate(., percent.cortex = percent.cortex |> as.character() ), ~(.)) |>
      when("percent.medulla" %in% colnames(.) ~ mutate(., percent.medulla = percent.medulla |> as.character() ), ~(.)) |>
      when("Age" %in% colnames(.) ~ mutate(., Age = Age |> as.numeric() ), ~(.)) |>
      when("nCount_RNA" %in% colnames(.) ~ mutate(., nCount_RNA = nCount_RNA |> as.numeric() ), ~(.)) |>
      when("is_primary_data" %in% colnames(.) ~ mutate(., is_primary_data = is_primary_data |> as.character() ), ~(.)) |>
      
      #mutate(across(contains("cluster", ignore.case = TRUE), ~ as.character)) |>
      select(-one_of('PCW')) %>%
      
      # Sort sample ID. It works but not elegant.
      # Based on observation of strangely behaving datasets, where sample ID is not clear
      when("sampleID" %in% colnames(.) & !"PatientID" %in% colnames(.) ~
             mutate(., PatientID = as.character(sampleID )) |>  select(-sampleID), ~(.)) %>%
      when("Patient" %in% colnames(.) ~ mutate(., Sample = NA |> as.character()), ~(.)) |>
      mutate(sample_placeholder = NA |> as.character()) %>%
      when(unique(.$dataset_id)=="e40591e7-0e5a-4bef-9b60-7015abe5b17f" ~ mutate(., sample_placeholder = glue("{batch} {development_stage}") |> as.character()), ~ (.)) %>%
      when(unique(.$dataset_id)=="39b6cc45-8c5c-4f7b-944c-58f66da5efb1" ~ mutate(., sample_placeholder =sample_id), ~ (.))  %>%
      when(unique(.$dataset_id)=="443d6a0e-dbcb-4002-8af0-628e7d4a18fa" ~ mutate(., sample_placeholder =sample_id), ~ (.))  %>%
      when(unique(.$dataset_id)=="a91f075b-52d5-4aa3-8ecc-86c4763a49b3" ~ mutate(., sample_placeholder =sample), ~ (.))  %>%
      when(unique(.$dataset_id)=="0af763e1-0e2f-4de6-9563-5abb0ad2b01e" ~ mutate(., sample_placeholder ="only_one_culture"), ~ (.))  %>%
      when(unique(.$dataset_id)=="5c64f247-5b7c-4842-b290-65c722a65952" ~ mutate(., sample_placeholder ="only_one_culture"), ~ (.))  %>%
      when(unique(.$dataset_id)=="d6f92754-e178-4202-b86f-0f430e965d72" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
      when(unique(.$dataset_id)=="c790ef7a-1523-4627-8603-d6a02f8f4877" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
      when(unique(.$dataset_id)=="1e81a742-e457-4fc6-9c39-c55189ec9dc2" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
      when(unique(.$dataset_id)=="351ef284-b59e-43a5-83ba-0eb907dc282c" ~ mutate(., sample_placeholder =orig.ident), ~ (.))  %>%
      when(unique(.$dataset_id)=="f498030e-246c-4376-87e3-90b28c7efb00" ~ mutate(., sample_placeholder =Name), ~ (.))  %>%
      
      # These are the datasets with too few cells per inferred samples, therefore simplifying
      when(unique(.$dataset_id)=="e3a56e00-8417-4d82-9d35-3fab3aac12f2" ~ mutate(., SpecimenID =NA), ~ (.))  %>%
      when(unique(.$dataset_id)=="17b34e42-bbd2-494b-bf32-b9229344a3f6" ~ mutate(., Sample =NA), ~ (.))  %>%
      
      
      # Fix huge samples for plate experiments
      tidyr::extract(cell_, "experiment___", "(^expr?[0-9]+)", remove = F) |>
      tidyr::extract(cell_, c("run_from_cell_id"), "(run[[:alnum:]_]+?)(?=_[ACGT]{5,}|-.+|$)", remove = FALSE) |> 
      
      mutate(experiment___ = if_else(dataset_id=="3fe53a40-38ff-4f25-b33b-e4d60f2289ef", experiment___, "")) |>
      mutate(run_from_cell_id = if_else(dataset_id=="d7f5d8d0-6150-48d7-b094-c34286ad11a1", "" , run_from_cell_id)) |>
      # If run-based embrio study get sample ID from cell ID
      
      
      # Empirically infer samples from many characteristics
      unite("sample_heuristic", one_of(
        "sample_placeholder",
        "Sample",
        "SampleID",
        "sample_uuid",
        "Sample_ID",
        "scRNASeq_sample_ID",
        "Sample_Tag",
        "Sample.ID",
        "sample_names",
        "Short_Sample",
        "Sample.ID.short",
        "Sample.name",
        "patient",
        "Donor.ID",
        "donor_id",
        "donor",
        "tech_sample",
        "PatientID",
        "donor_uuid",
        "library_uuid",
        "suspension_uuid",
        "Patient",
        "tissue_section_uuid",
        "DonorID",
        "specimen",
        "SpecimenID",
        "Fetus_id",
        "individual",
        "tissue",
        "development_stage",
        "assay",
        "experiment___",
        "disease",
        "run_from_cell_id",
        "is_primary_data"
      ), na.rm = TRUE, sep = "___", remove = F) |>
      
      
      # Add sample hash
      mutate(sample_ = getVDigest(algo="md5")(glue("{sample_heuristic}{dataset_id}"))) |>
      
      mutate(sample_id = sample_) |> 
      
      # make lighter
      mutate_if(is.character, as.factor)
    
  }
  
  get_metadata = function(.x){
    
    cache.path = "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/h5ad/2025-11-08/"
    
    dataset_id = as.character(unique(.x$dataset_id))
    
    h5_path = file.path(cache.path, paste0(dataset_id, ".h5ad"))
    
    sce = 
      h5_path |> 
      readH5AD(use_hdf5 = TRUE,  raw = FALSE, skip_assays = TRUE, layers=FALSE, reader = "R"	) 
    
    
    if(is.null(sce) || !"donor_id" %in% colnames(colData(sce)))
      sce = 
      h5_path |> 
      readH5AD(use_hdf5 = TRUE, raw = FALSE, skip_assays = TRUE, layers=FALSE	)
    
    metadata = 
      sce |> 
      mutate(dataset_id = dataset_id) |> 
      as_tibble() |>
      select(-contains("X_pca"), -contains("X_umap"), -contains("X_scVI")) |> 
      dplyr::rename(cell_ = .cell)
    
    metadata
    
    
    # # join the file metadata
    # column_to_omit_becuse_duplicated = 
    #   colnames(.x) |> 
    #   intersect(colData(sce) |> colnames()) |> 
    #   str_subset("donor_id", negate = TRUE) |> 
    #   c("embedding")
    # 
    # rm(sce)
    # gc(verbose = FALSE)
    # 
    # metadata = 
    #   metadata |> 
    #   left_join(
    #     .x |> 
    #       select(!any_of(column_to_omit_becuse_duplicated)) |> 
    #       unnest(donor_id) |> 
    #       unnest(donor_id) |> 
    #       select_if(negate(is.list)) ,
    #     by = join_by(donor_id)
    #   ) 
    
    # metadata = 
    #   metadata |> 
    #   sample_heuristics() 
    
    # # delete raw data
    # sample_column_to_preserve = 
    #   metadata |> 
    #   slice_sample(n = 500, by = donor_id) |> 
    #   tidybulk::pivot_sample(.sample = sample_) |> 
    #   colnames()
    
    # # Select only sample_ columns
    # metadata = 
    #   metadata |> 
    #   select(sample_, any_of(sample_column_to_preserve)) |> 
    #   distinct()
    
    
    
    # metadata
  }
  
  #-----------------------#
  # Pipeline
  #-----------------------#
  list(
    
    # Track file changes
    tar_file(
      census_meta_file,
      "/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/census_new_datasets_2025-11-08.parquet",
      deployment = "main"
    ),
    
    # New dataset ids
    tar_target(
      new_dataset_ids,
      # Metadata is created in ~/git_control/cellNexus/dev/STEP_1_census_dataset_ids_download_path.R
      arrow::read_parquet(census_meta_file) |> 
        pull(dataset_id) |> unique() |> as.character(),
      resources = tar_resources(crew = tar_resources_crew("slurm_1_10"))
    ),
    
    # Get rownames
    tar_target(
      my_db,
      db(overwrite=TRUE),
      resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
    ),
    
    # Get cellxgene annotation
    tar_target(
      dataset_df,
      datasets(my_db) |>
        # Only update new datasets 
        filter(dataset_id %in% new_dataset_ids) |>
        left_join(
          files(my_db) |> filter(filetype=="H5AD"),
          by = "dataset_id"
        ) |> 
        unnest(donor_id),
      resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
    ),
    
    tar_target(
      files_dataset_id,
      arrow::read_parquet(census_meta_file) |> 
        
        # # TESTING PURPOSE ONLY
        # #filter(dataset_id %in% c("d7f5d8d0-6150-48d7-b094-c34286ad11a1", "72955cdb-bd92-4135-aa52-21f33f9640db")) |>
        # filter(dataset_id %in% c("149b2c3f-ee11-47a7-984b-923570280bd7", "601c7454-17ff-432c-8e0d-902bfddb8013")) |>
        
        group_split(dataset_id),
      iteration = "list",
      resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
    ),
    
    # Get dataset_id cell index dictionary
    tar_target(
      dataset_cell_dict,
      get_metadata(files_dataset_id) |>
        mutate(cell_index = row_number(),
               new_cell_id = paste(dataset_id, cell_index, sep = "___"),
               cell_id_in_sample = paste(cell_, dataset_id, sep = "___")) |> 
        select(cell_, dataset_id, new_cell_id, cell_id_in_sample),
      pattern = map(files_dataset_id),
      iteration = "list",
      resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
    ),
    
    # Get SCE SMALL
    tar_target(
      metadata_dataset_id,
      get_metadata(files_dataset_id) |> 
        sample_heuristics() ,
      pattern = map(files_dataset_id),
      iteration = "list",
      resources = tar_resources(crew = tar_resources_crew("slurm_1_20"))
    ),
    
    # # select column that are present in half of the datasets at least, so the common column
    # tar_target(
    #   common_columns,
    #   metadata_dataset_id |>
    #     map_dfr(~ .x |> colnames() |> as_tibble()) |> 
    #     dplyr::count(value) |>
    #     mutate(n_datasets = length(metadata_dataset_id)) |>
    #     filter(n > (n_datasets / 2)) |>
    #     pull(value) ,
    #   resources = tar_resources(crew = tar_resources_crew("slurm_1_200"))
    # ),
    # 
    # tar_target(
    #   metadata_dataset_id_common_sample_columns,
    #   metadata_dataset_id |> 
    #     
    #     # Only get primary data
    #     # filter(is_primary_data=="TRUE") |> 
    #     
    #     mutate(cell_ =  as.character(cell_)) |> 
    #     select(any_of(common_columns)) |> 
    #     
    #     # Drop some clearly cell-wise columns
    #     select(-any_of(c("observation_joinid", "cell_")), -contains("cell_type"), -contains("X_pca")) |> 
    #     
    #     select_sample_columns(),
    #   pattern = map(metadata_dataset_id),
    #   resources = tar_resources(crew = tar_resources_crew("slurm_1_80") )
    # ),
    
    tar_target(
      metadata_dataset_id_cell_to_sample_mapping,
      metadata_dataset_id |> 
        
        # Only get primary data
        # filter(is_primary_data=="TRUE") |> 
        
        mutate(
          cell_ =  as.character(cell_), 
          observation_joinid = as.character(observation_joinid)
        ) |>
        select(
          cell_,
          sample_,
          sample_id,
          sample_heuristic,
          dataset_id,
          observation_joinid,
          assay,
          assay_ontology_term_id,
          disease,
          disease_ontology_term_id,
          donor_id,
          sex,
          sex_ontology_term_id,
          self_reported_ethnicity,
          self_reported_ethnicity_ontology_term_id,
          tissue,
          tissue_ontology_term_id,
          tissue_type,
          development_stage,
          development_stage_ontology_term_id,
          is_primary_data,
          suspension_type,
          cell_type,
          cell_type_ontology_term_id
        ),
      pattern = map(metadata_dataset_id),
      resources = tar_resources(crew = tar_resources_crew("slurm_1_80"))
    )
    
  )
  
  
}, 
ask = FALSE, 
script = glue("{result_directory}/_targets.R")
)

job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "summary",
    script = glue("{result_directory}/_targets.R"),
    store = glue("{result_directory}/_targets")
  )
  
})

# tar_workspace(metadata_dataset_id_a2a343b326a16a34, script = glue("{result_directory}/_targets.R"),
#               store = glue("{result_directory}/_targets"))
# x =  get_metadata(files_dataset_id)
# debugonce(sample_heuristics)
# x |> sample_heuristics()
# x |> mutate(cell_index = row_number(),
#             new_cell_id = paste(dataset_id, cell_index, sep = "___"),
#             cell_id_in_sample = paste(cell_, dataset_id, sep = "___")) 
# metadata_dataset_id = tar_read(metadata_dataset_id, store = glue("{result_directory}/_targets"))
# metadata_dataset_id$metadata_dataset_id_f9463ee2b53034d1 |> colnames()
#   


cellxgene_annotations <- readRDS("~/cellxgene_curated/metadata_cellxgenedp_Dec_2025/cellxgene_annotations_to_keep.rds")
file.copy("~/cellxgene_curated/metadata_cellxgenedp_Dec_2025/cellxgene_annotations_to_keep.rds", 
          "~/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cellxgene_annotations_to_keep.rds")

# dataset metadata
tar_read(dataset_df, store = glue("{result_directory}/_targets"))|>
  select(any_of(cellxgene_annotations)) |>
  unnest(donor_id) |> 
  # One donor in a dataset have more than one assay, get all assays.
  unnest(assay) |>
  mutate(assay = purrr::map_chr(assay, "label"))  |>
  write_parquet("~/cellxgene_curated/metadata_cellxgenedp_Jan_2026/dataset.parquet",
                compression = "zstd")

# # Sample to cell link
tar_read(metadata_dataset_id_cell_to_sample_mapping, store = glue("{result_directory}/_targets"))|>
  write_parquet("~/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_ids_for_metadata.parquet",
                compression = "zstd")

tar_read(dataset_cell_dict, store = glue("{result_directory}/_targets"))|>bind_rows() |>
  write_parquet("~/cellxgene_curated/metadata_cellxgenedp_Jan_2026/dataset_cell_dict.parquet",
                compression = "zstd")

get_tissue_grouped = function(tissue){
  
  list(
    
    # Respiratory System
    "respiratory system" = c(
      "lung", "lung parenchyma", "alveolus of lung",  "bronchus",
      "respiratory airway", "pleura", "pleural effusion", "middle lobe of right lung",
      "upper lobe of left lung", "lower lobe of left lung", "upper lobe of right lung",
      "lower lobe of right lung", "lingula of left lung", "right lung", "left lung"
    ),
    
    trachea = c( "epithelium of trachea", "trachea"),
    
    # Cardiovascular System
    "cardiovascular system" = c(
      "heart", "heart left ventricle", "heart right ventricle", "cardiac ventricle",
      "cardiac atrium", "right cardiac atrium", "left cardiac atrium", "apex of heart",
      "aorta", "coronary artery", 
      "venous blood", "anterior wall of left ventricle", "myocardium", "interventricular septum", "ventricular tissue", "basal zone of heart"
    ),
    
    vasculature = c("kidney blood vessel", "artery", "vein", "vasculature", "mesenteric artery"),
    # Umbilical Cord Blood
    "umbilical cord blood" = "umbilical cord blood",
    
    # Oesophagus
    "oesophagus" = c(
      "esophagus", "lower esophagus", "esophagus muscularis mucosa",
      "submucosal esophageal gland"
    ),
    
    # Stomach
    "stomach" = c(
      "stomach", "body of stomach", "cardia of stomach"
    ),
    
    # Small Intestine
    "small intestine" = c(
      "small intestine", "duodenum", "jejunum", "ileum"
    ),
    
    # Large Intestine
    "large intestine" = c(
      "large intestine", "colon", "left colon", "right colon",
      "sigmoid colon", "descending colon", "transverse colon",
      "ascending colon", "hepatic flexure of colon", "caecum",
      "rectum", "appendix", "vermiform appendix"
    ),
    
    # Digestive System (General)
    "digestive system (general)" = c(
      "intestine", "hindgut"
    ),
    
    # Nasal, Oral, and Pharyngeal Regions
    "nasal, oral, and pharyngeal regions" = c(
      "nasal cavity", "nasopharynx", "oral mucosa", "tongue", "anterior part of tongue",
      "posterior part of tongue", "gingiva", "nose", "saliva"
    ),
    
    # Cerebral Lobes and Cortical Areas
    "cerebral lobes and cortical areas" = c(
      "frontal lobe", "left frontal lobe", "right frontal lobe", "primary motor cortex",
      "dorsolateral prefrontal cortex", "superior frontal gyrus", "orbitofrontal cortex",
      "medial orbital frontal cortex", "Broca's area", "prefrontal cortex",
      "temporal lobe", "left temporal lobe", "right temporal lobe", 
      "angular gyrus", "entorhinal cortex",
      "parietal lobe", "left parietal lobe", "right parietal lobe", "primary somatosensory cortex",
      "occipital lobe", "right occipital lobe", "primary visual cortex",
      "occipital cortex", "insular cortex", "parietal cortex", "temporal cortex",
      "frontal cortex", "Brodmann (1909) area 4", "temporoparietal junction",
      "middle temporal gyrus", "cingulate cortex", "brain", "brain white matter", "cerebral cortex", "cerebral nuclei"
    ),
    
    # Limbic and Basal Systems
    "limbic and basal systems" = c(
      "anterior cingulate cortex", "anterior cingulate gyrus", "hippocampal formation",
      "hypothalamus", "thalamic complex", "dentate nucleus", "basal ganglion",
      "caudate nucleus", "putamen", "substantia nigra pars compacta",
      "lateral ganglionic eminence", "medial ganglionic eminence",
      "caudal ganglionic eminence", "ganglionic eminence"
    ),
    
    # Brainstem and Cerebellar Structures
    "brainstem and cerebellar structures" = c(
      "pons", "midbrain", "myelencephalon", "telencephalon", "forebrain",
      "cerebellum", "cerebellum vermis lobule", "cerebellar cortex",
      "hemisphere part of cerebellar posterior lobe", "white matter of cerebellum"
    ),
    
    # General Brain and Major Structures
    "general brain and major structures" = c(
      "spinal cord", "neural tube", "cervical spinal cord white matter"
    ),
    
    # Muscular System (Skeletal Muscles)
    "muscular system (skeletal muscles)" = c(
      "rectus abdominis muscle", "gastrocnemius", "muscle of abdomen", "muscle organ",
      "muscle tissue", "pelvic diaphragm muscle", "skeletal muscle tissue", "muscle of pelvic diaphragm"
    ),
    
    # Connective Tissue
    "connective tissue" = c(
      "connective tissue", "tendon of semitendinosus", "vault of skull", "bone spine",
      "rib"
    ),
    
    # Adipose Tissue
    "adipose tissue" = c(
      "adipose tissue", "subcutaneous adipose tissue", "visceral abdominal adipose tissue",
      "perirenal fat", "omental fat pad", "subcutaneous abdominal adipose tissue",
      "abdominal adipose tissue"
    ),
    
    # Endocrine System
    "endocrine system" = c(
      "thyroid gland", "adrenal tissue", "adrenal gland", "islet of Langerhans",
      "endocrine pancreas", "pineal gland"
    ),
    
    # Lymphatic System
    "lymphatic system" = c(
      "lymph node", "mesenteric lymph node", "thoracic lymph node",
      "cervical lymph node", "bronchopulmonary lymph node", "tonsil", "inguinal lymph node"
    ),
    
    # Integumentary System (Skin)
    "integumentary system (skin)" = c(
      "skin of abdomen", "skin of forearm", "skin of scalp", "skin of face", "skin of leg",
      "skin of chest", "skin of back", "skin of hip", "skin of body", "skin of cheek",
      "skin of temple", "skin of shoulder", "skin of external ear", "skin of trunk",
      "skin of prepuce of penis", "skin epidermis", "arm skin", "lower leg skin",
      "hindlimb skin", "zone of skin", "dermis", "skin of nose", "skin of forehead",
      "skin of pes", "axilla"
    ),
    
    # Gastrointestinal Accessory Organs
    "gallbladder" =  "gallbladder",
    
    # Gastrointestinal Accessory Organs
    "pancreas" = c( "pancreas", "exocrine pancreas" ),
    
    # Gastrointestinal Accessory Organs
    "liver" = c( "liver", "caudate lobe of liver", "hepatic cecum" ),
    
    # Spleen
    "spleen" = "spleen",
    
    # Thymus
    "thymus" = "thymus",
    
    # Blood
    "blood" = "blood",
    
    # Bone Marrow
    "bone marrow" = "bone marrow",
    
    # Female Reproductive System
    "female reproductive system" = c(
      "uterus", "myometrium", "fallopian tube", "ampulla of uterine tube",
      "fimbria of uterine tube", "uterine cervix", "endometrium",
      "decidua", "decidua basalis", "placenta", "yolk sac", "isthmus of fallopian tube"
    ),
    "ovary" = "ovary", 
    
    # Male Reproductive System
    "male reproductive system (other)" = c(
      "testis", "gonad"
    ),
    
    # Prostate
    "prostate" = c(
      "prostate gland", "transition zone of prostate", "peripheral zone of prostate"
    ),
    
    # Renal System
    "renal system" = c(
      "kidney", "cortex of kidney", "renal medulla", "renal papilla",
      "renal pelvis", "ureter", "bladder organ"
    ),
    
    # Miscellaneous Glands
    "miscellaneous glands" = c(
      "parotid gland", "lacrimal gland", "sublingual gland", "mammary gland",
      "chorionic villus"
    ),
    
    # Epithelium and Mucosal Tissues
    "epithelium and mucosal tissues" = c(
      "epithelium of small intestine", "epithelium of esophagus", "caecum epithelium",
      "jejunal epithelium", "ileal epithelium", "colonic epithelium",
      "submucosa of ascending colon", "submucosa of ileum", "lamina propria",
      "lamina propria of large intestine", "lamina propria of small intestine",
      "mucosa", "mucosa of colon", "lamina propria of mucosa of colon"
    ),
    
    # Eye and Visual-Related Structures
    "sensory-related structures" = c(
      "retina",
      "retinal neural layer",
      "macula lutea",
      "macula lutea proper",
      "sclera",
      "trabecular meshwork",
      "conjunctiva",
      "pigment epithelium of eye",
      "cornea",
      "iris",
      "ciliary body",
      "peripheral region of retina",
      "eye trabecular meshwork",
      "perifoveal part of retina",
      "choroid plexus",
      "lens of camera-type eye",
      "corneo-scleral junction",
      "fovea centralis",
      "eye",
      "inner ear",
      "vestibular system",
      "primary auditory cortex"
    ),
    
    # Digestive Tract Junctions and Connections
    "digestive tract junctions and connections" = c(
      "esophagogastric junction", "duodeno-jejunal junction", "hepatopancreatic ampulla",
      "hepatopancreatic duct", "pyloric antrum"
    ),
    
    # Peritoneal and Abdominal Cavity Structures
    "peritoneal and abdominal cavity structures" = c(
      "peritoneum", "omentum", "retroperitoneum", "mesentery"
    ),
    
    # Breast
    "breast" = c(
      "breast", "upper outer quadrant of breast"
    )
  ) |> 
    enframe(name ="tissue_groups") |> 
    distinct() |> 
    unnest(value) |> 
    dplyr::rename(tissue = value) |> 
    mutate()
  
  # #check
  # distinct_tissue = 
  #   tissue |> 
  #   enframe(name = "tissue") |> 
  #   distinct(tissue)
  # 
  # if(nrow(distinct_tissue) != distinct_tissue |> left_join(tissue_grouped_df, copy = TRUE))
  # 
  # 
  # tissue |> 
  #   enframe(name = "tissue") |> 
  #   left_join(tissue_grouped_df)
}

convert_age_labels_to_days <- function(labels) {
  age_days <- rep(NA, length(labels))
  
  carnegie_stages <- c(
    '1' = 1, '2' = 2, '3' = 4, '4' = 6, '5' = 8,
    '6' = 12, '7' = 16, '8' = 18, '9' = 20,
    '10' = 22, '11' = 24, '12' = 26, '13' = 28,
    '14' = 32, '15' = 35, '16' = 37, '17' = 41,
    '18' = 44, '19' = 46, '20' = 49, '21' = 51,
    '22' = 53, '23' = 56
  )
  
  word_to_num <- c(
    'first' = 0, 'second' = 10, 'third' = 20, 'fourth' = 30,
    'fifth' = 40, 'sixth' = 50, 'seventh' = 60, 'eighth' = 70,
    'ninth' = 80, 'tenth' = 90
  )
  
  word_ordinal_to_num <- c(
    'first' = 1, 'second' = 2, 'third' = 3, 'fourth' = 4,
    'fifth' = 5, 'sixth' = 6, 'seventh' = 7, 'eighth' = 8,
    'ninth' = 9, 'tenth' = 10, 'eleventh' = 11, 'twelfth' = 12,
    'thirteenth' = 13, 'fourteenth' = 14, 'fifteenth' = 15,
    'sixteenth' = 16, 'seventeenth' = 17, 'eighteenth' = 18,
    'nineteenth' = 19, 'twentieth' = 20, 'twenty-first' = 21,
    'twenty-second' = 22, 'twenty-third' = 23, 'twenty-fourth' = 24,
    'twenty-fifth' = 25, 'twenty-sixth' = 26, 'twenty-seventh' = 27,
    'twenty-eighth' = 28, 'twenty-ninth' = 29, 'thirtieth' = 30,
    'thirty-first' = 31, 'thirty-second' = 32, 'thirty-third' = 33,
    'thirty-fourth' = 34
  )
  
  stage_to_age <- list(
    'newborn human' = 0,
    'newborn' = 14,           # midpoint 0-28 days
    'infant' = 0.5 * 365,
    'nursing' = 0.5 * 365,    # 0-11 months
    'child' = 2.5 * 365,      # 1-4 yo midpoint
    'pediatric' = 6 * 365,
    'juvenile' = 9.5 * 365,   # 5-14 yo midpoint
    'adolescent' = 15 * 365,
    'young adult' = 25 * 365,
    'human early adulthood' = 25 * 365,
    'prime adult' = 30 * 365,
    'adult' = 40 * 365,
    'human adult' = 40 * 365,
    'human middle aged' = 50 * 365,
    'middle aged' = 50 * 365,
    'mature' = 40 * 365,
    'immature' = 1 * 365,
    'late adult' = 70 * 365,
    'human late adulthood' = 70 * 365,
    'human aged' = 75 * 365,
    'postnatal' = 1 * 365,
    'prenatal' = 28,
    'embryonic human' = 28,
    'embryonic' = 28,
    'organogenesis' = 28,
    'blastula' = 5,
    'unknown' = NA
  )
  
  for (i in seq_along(labels)) {
    label <- trimws(labels[i])
    age <- NA
    
    # 1. Unknown
    if (grepl("^unknown$", label, ignore.case = TRUE)) {
      age <- NA
      
      # 2. [N]-year-old (human )stage  — with or without "human"
    } else if (grepl("^(\\d+)-year-old( human)? stage$", label)) {
      num <- as.numeric(sub("^(\\d+)-year-old.*$", "\\1", label))
      age <- num * 365
      
      # 3. [N]-month-old (human )stage
    } else if (grepl("^(\\d+)-month-old( human)? stage$", label)) {
      num <- as.numeric(sub("^(\\d+)-month-old.*$", "\\1", label))
      age <- num * 30
      
      # 4. [N]th week post-fertilization (human )stage  — numeric ordinal suffix
    } else if (grepl("^(\\d+)(?:st|nd|rd|th) week post-fertilization( human)? stage$", label, perl = TRUE)) {
      num <- as.numeric(sub("^(\\d+)(?:st|nd|rd|th).*$", "\\1", label, perl = TRUE))
      age <- num * 7
      
      # 5. [word] week post-fertilization (human )stage  — ordinal word
    } else if (grepl("^([a-z-]+) week post-fertilization( human)? stage$", label, ignore.case = TRUE)) {
      ord_word <- tolower(sub("^([a-z-]+) week post-fertilization.*$", "\\1", label))
      if (ord_word %in% names(word_ordinal_to_num)) {
        age <- word_ordinal_to_num[[ord_word]] * 7
      }
      
      # 6. Carnegie stage [N] (handles leading zeros by coercing to integer)
    } else if (grepl("^Carnegie stage (\\d+)$", label)) {
      num <- as.character(as.integer(sub("^Carnegie stage (\\d+)$", "\\1", label)))
      age <- if (num %in% names(carnegie_stages)) carnegie_stages[[num]] else NA
      
      # 7. [N]-[N] year-old (human |child )?stage
    } else if (grepl("^(\\d+)-(\\d+) year-old( human| child)? stage$", label)) {
      num1 <- as.numeric(sub("^(\\d+)-.*$", "\\1", label))
      num2 <- as.numeric(sub("^\\d+-(\\d+).*$", "\\1", label))
      age <- ((num1 + num2) / 2) * 365
      
      # 8. [N]-[N] year-old (no "stage" suffix, e.g. "15-19 year-old")
    } else if (grepl("^(\\d+)-(\\d+) year-old$", label)) {
      num1 <- as.numeric(sub("^(\\d+)-.*$", "\\1", label))
      num2 <- as.numeric(sub("^\\d+-(\\d+).*$", "\\1", label))
      age <- ((num1 + num2) / 2) * 365
      
      # 9. under-1-year-old (human )stage
    } else if (grepl("^under-1-year-old( human)? stage$", label)) {
      age <- 0.5 * 365
      
      # 10. [ordinal word] LMP month (human )stage
    } else if (grepl("^([a-z-]+) LMP month( human)? stage$", label, ignore.case = TRUE)) {
      ord_word <- tolower(sub("^([a-z-]+) LMP month.*$", "\\1", label))
      if (ord_word %in% names(word_ordinal_to_num)) {
        age <- word_ordinal_to_num[[ord_word]] * 30
      }
      
      # 11. [ordinal word] decade (human )stage
    } else if (grepl("^([a-z]+) decade( human)? stage$", label, ignore.case = TRUE)) {
      decade_word <- tolower(sub("^([a-z]+) decade.*$", "\\1", label))
      if (decade_word %in% names(word_to_num)) {
        num1 <- word_to_num[[decade_word]]
        num2 <- num1 + 9
        age <- ((num1 + num2) / 2) * 365
      }
      
      # 12. [N] year-old and over (human )stage
    } else if (grepl("^(\\d+) year-old and over( human)? stage$", label)) {
      num <- as.numeric(sub("^(\\d+).*$", "\\1", label))
      age <- num * 365
      
      # 13. Named stages with parenthetical ranges, e.g. "newborn stage (0-28 days)"
      #     Strip the qualifier and look up the base name
    } else if (grepl("^(.*?) stage \\(.*\\)$", label)) {
      stage <- tolower(sub("^(.*?) stage \\(.*\\)$", "\\1", label))
      age <- if (stage %in% names(stage_to_age)) stage_to_age[[stage]] else NA
      
      # 14. Plain "[name] stage"
    } else if (grepl("^(.*) stage$", label)) {
      stage <- tolower(sub("^(.*) stage$", "\\1", label))
      age <- if (stage %in% names(stage_to_age)) stage_to_age[[stage]] else NA
      
      # 15. Default
    } else {
      age <- NA
    }
    
    age_days[i] <- age
  }
  
  return(as.integer(age_days))
}

age_days_tbl = 
  tbl(
    dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
    sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_ids_for_metadata.parquet')")
  ) |> 
  select(-contains("joinid"), -contains("cell_")) |>
  distinct(development_stage) |> 
  as_tibble() |> 
  mutate(age_days = convert_age_labels_to_days(development_stage)) 

age_days_tbl |>
  write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/age_days.parquet")


tissues_grouped = get_tissue_grouped() 

tissues_grouped |>
  write_parquet("/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/tissue_grouped.parquet")

# #
# cell_ids_for_metadata <- tbl(
#   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_ids_for_metadata.parquet')")
# )
# 
# # this metadata is generated in ~/git_control/cellNexus/dev/STEP_3_split_large_datasets_create_samples.R
# cell_to_refined_sample_from_Mengyuan <- tbl(
#   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   sql(glue::glue("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/{Date}_census_samples_to_download.parquet')"))
# ) |>
#   select(cell_, observation_joinid, dataset_id, sample_id) 
# 
# 
# # Establish a connection to DuckDB in memory
# job::job({
#   
#   con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
#   
#   # Create views for each of the datasets in DuckDB
# #   dbExecute(con, "
# #   CREATE VIEW cell_to_refined_sample_from_Mengyuan AS
# #   SELECT cell_, observation_joinid, dataset_id, sample_id, cell_type, cell_type_ontology_term_id
# #   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/2025-11-08_census_samples_to_download.parquet')
# # ")
#   
#   dbExecute(con, "
#   CREATE VIEW cell_ids_for_metadata AS
#   SELECT *
#   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_ids_for_metadata.parquet')
# ")
#   
#   dbExecute(con, "
#   CREATE VIEW sample_metadata AS
#   SELECT *
#   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/dataset.parquet')
# ")
#   
#   dbExecute(con, "
#   CREATE VIEW age_days_tbl AS
#   SELECT development_stage, age_days
#   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/age_days.parquet')
# ")
#   
#   dbExecute(con, "
#   CREATE VIEW tissue_grouped AS
#   SELECT tissue, tissue_groups
#   FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/tissue_grouped.parquet')
# ")
#   
#   # Perform optimised joins within DuckDB
#   copy_query <- "
# COPY (
#   SELECT 
#     cell_ids_for_metadata.*,
#     sample_metadata.*,
#     age_days_tbl.age_days,
#     tissue_grouped.tissue_groups,
#     'cellxgene' AS atlas_id
#   
#   FROM cell_ids_for_metadata
#   
#   -- LEFT JOIN cell_ids_for_metadata
#   --   ON cell_ids_for_metadata.cell_ = cell_to_refined_sample_from_Mengyuan.cell_
#   --   AND cell_ids_for_metadata.observation_joinid = cell_to_refined_sample_from_Mengyuan.observation_joinid
#   --   AND cell_ids_for_metadata.dataset_id = cell_to_refined_sample_from_Mengyuan.dataset_id
#   
#   LEFT JOIN sample_metadata
#     ON sample_metadata.dataset_id = cell_ids_for_metadata.dataset_id
#     AND sample_metadata.donor_id = cell_ids_for_metadata.donor_id
#     
#   LEFT JOIN age_days_tbl
#     ON age_days_tbl.development_stage = cell_ids_for_metadata.development_stage
# 
#   LEFT JOIN tissue_grouped
#     ON tissue_grouped.tissue = cell_ids_for_metadata.tissue
#     
# ) TO '/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet'
# (FORMAT PARQUET, COMPRESSION 'gzip');
# "
#   
#   # Execute the final query to write the result to a Parquet file
#   dbExecute(con, copy_query)
#   
#   # Disconnect from the database
#   dbDisconnect(con, shutdown = TRUE)
#   
# })
# 
# # system("~/bin/rclone copy /vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/")
# 
# 
# cell_metadata = tbl(
#   dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
#   sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/cell_metadata.parquet')")
# ) 



##############
#  PLOTS     #
##############

cell_metadata |> 
  left_join(tissues_grouped, copy = TRUE) |> 
  mutate(age_days = convert_age_labels_to_days(development_stage))
distinct(donor_id, tissue_groups) |> 
  ggplot(aes(fct_infreq(tissue_groups))) +
  geom_bar() +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

cell_metadata |> 
  left_join(tissues_grouped, copy = TRUE) |> 
  left_join(days_df, copy = TRUE) |> 
  filter(age_days > 365) |> 
  distinct(sample_id, donor_id, tissue_groups) |> 
  as_tibble() |> 
  nest(data = -tissue_groups) |> 
  mutate(n_sample = map_int(data, ~ .x |> distinct(sample_id) |> nrow())) |> 
  mutate(n_donor = map_int(data, ~ .x |> distinct(donor_id) |> nrow())) |> 
  ggplot(aes(n_sample, n_donor)) +
  geom_point() +
  geom_text(aes(label = tissue_groups)) +
  scale_y_log10() +
  scale_x_log10() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

cell_metadata |> 
  left_join(tissues_grouped, copy = TRUE) |> 
  left_join(days_df, copy = TRUE) |> 
  write_parquet("/vast/scratch/users/shen.m/cellNexus_run/cell_metadata_temp.parquet")



# Disconnect from the database
dbDisconnect(con, shutdown = TRUE)


# Get Dharmesh metadata consensus
# system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/consensus_output_new.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/ ")
# system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/data_driven_consensus_new.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/ ")
# system("~/bin/rclone copy box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/data_driven_consensus.parquet /vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/ ")


# (Mengyuan: UNSURE ABOUT THE PURPOSE OF CODE BELOW)
# Non immune harmonisation to Dharmesh immune harmonisation
non_immune_harmonisation = 
  read_csv("/vast/projects/mangiola_immune_map/PostDoc/CuratedAtlasQueryR/dev/cell_type_harmonisation_non_immune.csv") 

# system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/CuratedAtlasQueryR/dev/cell_type_harmonisation_non_immune.csv box_adelaide:/Mangiola_ImmuneAtlas/reannotation_consensus/")


tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/consensus_output_new.parquet')")
) |> 
  
  # Add non immune harmonisation to Dharmesh immune harmonisation
  mutate(is_immune = reannotation_consensus == "non immune") |> 
  left_join(non_immune_harmonisation, copy = TRUE
  ) |> 
  mutate(reannotation_consensus = case_when(reannotation_consensus=="non immune" ~ non_immune_harmonised_cell_type, TRUE ~ reannotation_consensus)) |> 
  select(-non_immune_harmonised_cell_type) |> 
  cellNexus:::duckdb_write_parquet(path = "/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/consensus_output_new_plus_non_immune_harmonisation.parquet")


annotation_with_harmonised <- tbl(
  dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
  sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgenedp_Jan_2026/consensus_output_new_plus_non_immune_harmonisation.parquet')")
)


# 
# # a test to see whether donor ID is present in the new metadata
# test = 
#   files_metadata |> 
#   slice(1:50) |> 
#   nest(data = c(dataset_id, dataset_version_id, filetype, url)) |> 
#   mutate(has_donor_id = map_lgl(
#     data,
#     ~ {
#       browser()
#       h5_path = .x |> files_download(dry.run = FALSE)
#       has_donor_id = 
#         h5_path |> 
#         readH5AD(use_hdf5 = TRUE	) |> 
#         colData() |> 
#         as_tibble() |> 
#         select(any_of("donor_id")) |> 
#         ncol() >
#         0
#       file.remove(h5_path)
#       has_donor_id
#     }
#   )) |> 
#   unnest(data) |> 
#   select(dataset_version_id, has_donor_id)
# 
# 


# files_metadata |>
#   
#   # Get organism list and filter human
#   mutate(organism_name = map_chr(organism, ~ .x |> map(~.x$label) |> paste(collapse=", ") )) |>
#   filter(organism_name |> str_detect("Homo sapiens")) |>
#   
#   # Download
#   files_download(dry.run = FALSE, cache_path = "{root_directory}/raw_data/") |>
#   
#   # Save file list
#   saveRDS("{root_directory}/file_location.rds")
# 

