summary_store = "/vast/scratch/users/shen.m/calculate_census_raw_counts_target_store"

tar_script({
  library(dplyr)
  library(SummarizedExperiment)
  library(zellkonverter)
  library(crew)
  library(crew.cluster)
  library(duckdb)
  
  # Helper to avoid repetition
  new_elastic <- function(name, mem_gb, time_min, workers, crashes_max, cpus_per_task = 2, backup = NULL) {
    crew_controller_slurm(
      name = name,
      workers = workers,
      crashes_max = crashes_max,
      seconds_idle = 30,
      options_cluster = crew_options_slurm(
        memory_gigabytes_required = mem_gb,
        cpus_per_task = cpus_per_task,
        time_minutes = time_min
      ),
      backup = backup
    )
  }
  
  elastic_300      <- new_elastic("elastic_300",      300, 60 * 24, workers = 8,   crashes_max = 2)
  elastic_160      <- new_elastic("elastic_160",      160, 60 * 24, workers = 8,   crashes_max = 2, backup = elastic_300)
  elastic_120      <- new_elastic("elastic_120",      120, 60 * 4,  workers = 16,  crashes_max = 1, cpus_per_task = 1, backup = elastic_160)
  elastic_80       <- new_elastic("elastic_80",        80, 60 * 4,  workers = 24,  crashes_max = 1, cpus_per_task = 1, backup = elastic_120)
  elastic_40       <- new_elastic("elastic_40",        40, 60 * 4,  workers = 32,  crashes_max = 1, cpus_per_task = 1, backup = elastic_80)
  elastic_20       <- new_elastic("elastic_20",        20, 60 * 4,  workers = 48,  crashes_max = 1, cpus_per_task = 1, backup = elastic_40)
  elastic_10       <- new_elastic("elastic_10",        10, 60 * 4,  workers = 150, crashes_max = 2, cpus_per_task = 1, backup = elastic_20)
  elastic_5_minimal <- new_elastic("elastic_5_minimal", 5, 60 * 4,  workers = 300, crashes_max = 2, cpus_per_task = 1, backup = elastic_10)
  
  controllers <- crew_controller_group(
    elastic_10, elastic_20, elastic_40, elastic_80, elastic_120, elastic_160, elastic_300, elastic_5_minimal
  )
  
  tar_option_set(
    memory             = "transient",
    garbage_collection = 100,
    storage            = "worker",
    retrieval          = "worker",
    error              = "continue",
    cue                = tar_cue(mode = "thorough"),
    format             = "qs",
    workspace_on_error = TRUE,
    controller         = controllers,
    trust_object_timestamps = TRUE,
    resources = tar_resources(
      crew = tar_resources_crew(controller = "elastic_5_minimal")
    )
  )
  
  # ── helpers ────────────────────────────────────────────────────────────────
  
  pos_min_med_ratio <- function(x) {
    pos <- x[x > 0]
    min(pos) / median(pos)
  }
  
  pos_min_mean_ratio <- function(x) {
    pos <- x[x > 0]
    min(pos) / mean(pos)
  }
  
  get_positive_mode <- function(x) {
    sort(table(x[x > 0]), decreasing = TRUE)[1] |> names() |> as.numeric()
  }
  
  # Stage 1 – read SCE and return a minimal list with only what downstream needs.
  # Keeping the raw matrix avoids re-reading the file in the metrics stage, while
  # returning a plain list (not a full SCE) keeps the serialised object small.
  read_sce_counts <- function(file) {
    sce        <- readH5AD(file, reader = "R", use_hdf5 = TRUE)
    if (ncol(sce) == 0) return(NULL)
    
    assay_name <- names(sce@assays)[1]
    counts_mat <- as.matrix(assay(sce, assay_name))
    
    list(
      sample_id  = basename(file),
      counts_mat = counts_mat,
      n_cells    = ncol(counts_mat),
      n_genes    = nrow(counts_mat)
    )
  }
  
  # Stage 2 – pure computation on the pre-loaded matrix; no I/O.
  calc_counts_metrics <- function(sce_data) {
    if (is.null(sce_data)) return(NULL)
    
    counts_vec <- as.numeric(sce_data$counts_mat)
    tol        <- 1e-4
    all_int    <- all(counts_vec == floor(counts_vec), na.rm = TRUE)
    
    tibble::tibble(
      sample_id          = sce_data$sample_id,
      min_val            = min(counts_vec,    na.rm = TRUE),
      median_val         = median(counts_vec, na.rm = TRUE),
      max_val            = max(counts_vec,    na.rm = TRUE),
      counts_gap_min_med = pos_min_med_ratio(counts_vec),
      counts_gap_min_mean= pos_min_mean_ratio(counts_vec),
      positive_mode      = get_positive_mode(counts_vec),
      has_negative       = min(counts_vec, na.rm = TRUE) < 0,
      max_gt_10          = max(counts_vec, na.rm = TRUE) > 10,
      all_integer        = all_int,
      has_floating       = !all_int && all(abs(counts_vec - round(counts_vec)) < tol, na.rm = TRUE),
      n_cells            = sce_data$n_cells,
      n_genes            = sce_data$n_genes
    ) |>
      left_join(
        tbl(
          dbConnect(duckdb::duckdb(), dbdir = ":memory:"),
          sql("SELECT * FROM read_parquet('/vast/projects/cellxgene_curated/metadata_cellxgene_mengyuan/census_samples_to_download_groups_MODIFIED.parquet')")
        ) |>
          select(-observation_joinid, list_length) |>
          distinct() |>
          mutate(sample_id = paste0(sample_2, ".h5ad")) |>
          collect(),
        by   = "sample_id",
        copy = TRUE
      )
  }
  
  # ── pipeline ───────────────────────────────────────────────────────────────
  
  list(
    tar_target(
      files,
      # {
      #   test_samples <- c(
      #     "dc30c21e4ca4563aef9493cbd9e4b586.h5ad",
      #     "087dd99243884e159d4ac04da9b4f6f3.h5ad",
      #     "6b4de4d5b05dd96303d0f979815cc142___-1_2_1_1_1_1_1_1_1_1_1_1_1.h5ad",
      #     "38937b107daa51abc815831ce3e8dd24.h5ad",
      #     "534095645762eaeb261a55aed943480b___0109_cv_vax_1_1.h5ad",
      #     "aee0db17365492b18701a47660f3e440.h5ad",
      #     "a52d2d1130700a9b8d8d124b6476376e.h5ad",
      #     "7960a558de32e9878edaa899c97c80dd.h5ad",
      #     "de3aad3f8fe99e636e82295650ab1016.h5ad",
      #     "b7ddd9ea2220a46492819907b1d0ea19.h5ad",
      #     "0e1bc953c54f0140b83a2b19c6aece78.h5ad",
      #     "674f201eb6a1de59241c2dde6eaee50a.h5ad",
      #     "41f21a039b475ed214e639e6e508c458.h5ad",
      #     "1a53aa2f40a4b0516eb018bf7892d503.h5ad",
      #     "35b1ca9b2912083fc9d8379dcaff0906.h5ad",
      #     "70461d10b92e386b22d759d81707bb76.h5ad",
      #     "a39b7c9ed193a2ca7b124a0c677d70cc.h5ad",
      #     "ade147d91712a9217066cb68c8b2606e.h5ad",
      #     "6765c944aeb997e6570d154f23852380___-23L8TX_191105_01_D09.h5ad",
      #     "ee71b036aa9944b4809615152e49d542___A57_.h5ad",
      #     "1fd673dad737d111f68441ae1ed1aa9e.h5ad",
      #     "5d61ebd7daac1a8e96c454adecb19e8c.h5ad",
      #     "990182fecfe9aab5f97d6e394a595919.h5ad",
      #     "9e8b5ab0f935b61fead3eee13687d395.h5ad",
      #     "3f2358dd1d90fe5d8afee467e0f0ab49.h5ad",
      #     "26798f69434f7fce86212348d4392d93.h5ad",
      #     "42711051a1fdad5a059a700ff34b4253.h5ad",
      #     "44fc2f063662d954dec29d9a48f8a0e0.h5ad",
      #     "c1d4e349bd6614ee903691ebc3123ac4.h5ad",
      #     "1cadcd59e1a5f3d456ee7c23f88b8036____12.h5ad",
      #     "8f1e8028f50bef82122f32950354d2ff.h5ad",
      #     "2eca1037df233359bb1767cd4f334ceb.h5ad",
      #     "850e404ead05edb4d1fbc31977db67ae.h5ad",
      #     "11e339b4543c982c9c41f5daa15f343a___-Per_R_L3-HUMAN.h5ad",
      #     "c4105c78b7383a7ff41c5dfc4942c283.h5ad",
      #     "2f13234a080a8926e896f7161e3a6961.h5ad",
      #     "81ceba2757ba91991add932a95c04cd6.h5ad",
      #     "df10139747ddb1ce024d34f922d4c376___-LKTX_200924_01_G01-1078639226.h5ad",
      #     "1dca04a3ae829a506968f1a950b89871.h5ad",
      #     "cb3b213ef867e166880261e6cb44ea3f.h5ad",
      #     "5d88e53d4d2a58ae1176bd9121aef59c.h5ad",
      #     "ea212026be0a3aaa2cbaec7793af4d15___E2L4_.h5ad",
      #     "94ad3ca9fbcab46a5907393253a77745___-L8TX_200709_01_F03-1047273139.h5ad",
      #     "218cf4ff55bc6102c80afc33e67357d2.h5ad",
      #     "8f91fabe44a45ad86d61a01645c97674___FCAImmP7352191-.h5ad"
      #   )
      #   all_files <- list.files(
      #     "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01",
      #     full.names = TRUE, pattern = "\\.h5ad$"
      #   )
      #   all_files[basename(all_files) %in% test_samples]
      # },
      list.files(
        "/vast/scratch/users/shen.m/Census/split_h5ad_based_on_sample_id/2024-07-01",
        full.names = TRUE, pattern = "\\.h5ad$"
      ),
      deployment = "main"
    ),
    
    # Stage 1 – I/O-bound; needs memory for the full matrix
    tar_target(
      sce_counts,
      read_sce_counts(files),
      pattern   = map(files),
      iteration = "list",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_5_minimal")
      )
    ),
    
    # Stage 2 – CPU-bound, no disk I/O; can run on smaller workers
    tar_target(
      sample_summary_df,
      calc_counts_metrics(sce_counts),
      pattern   = map(sce_counts),
      iteration = "list",
      resources = tar_resources(
        crew = tar_resources_crew(controller = "elastic_10")
      )
    )
  )
  
}, ask = FALSE, script = glue("{summary_store}/_targets.R"))

job::job({
  
  tar_make(
    # callr_function = NULL,
    reporter = "summary",
    script = glue("{summary_store}/_targets.R"),
    store = glue("{summary_store}/_targets")
  )
  
})

sample_summary_df = tar_read(sample_summary_df,  store = glue("{summary_store}/_targets")) |> bind_rows() |> 
  mutate(max_gt_20 = ifelse(max_val > 20, TRUE, FALSE))
sample_summary_df |> saveRDS("~/test_sample_summary_df.rds")
sample_summary_df <- readRDS("~/test_sample_summary_df.rds")

plot_raw_counts_hist_pdf <- function(
    samples_to_plot,
    pdf_file,
    width = 8,
    height = 8,
    nrow = 3,
    ncol = 3,
    max_cells = 5e3,
    hist_ylim = c(0, 1e6)
) {
  
  pdf(pdf_file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  
  par(mfrow = c(nrow, ncol), mar = c(3, 3, 2, 1))
  
  samples_to_plot |>
    dplyr::transmute(
      file_name,
      counts_gap = as.character(counts_gap_min_mean)
    ) |>
    purrr::pmap(
      \(file_name, counts_gap) {
        
        sce <- zellkonverter::readH5AD(
          file = file_name,
          reader = "R",
          use_hdf5 = TRUE
        )
        
        if (ncol(sce) == 0) return(NULL)
        
        if (ncol(sce) > max_cells) {
          sce <- sce[, sample(ncol(sce), max_cells)]
        }
        
        assay_name <- names(SummarizedExperiment::assays(sce))[1]
        
        x <- SummarizedExperiment::assay(sce, assay_name) |>
          as.numeric()
        
        hist(
          x,
          main = counts_gap,
          ylim = hist_ylim,
          breaks = 100
        )
        
        invisible(NULL)
      },
      .progress = TRUE
    )
  
  par(mfrow = c(1, 1))
}

# use counts_gap_min_mean 0.25
samples_to_plot1 = sample_target_tbl |> left_join(sample_summary_df, by = c("sample_id" = "sample_2")) |> 
  filter(counts_gap_min_mean >= 0.25) |>
  group_by(dataset_id) |> 
  slice_head(n=1) |>
  ungroup()

samples_to_plot1 = samples_to_plot1 |>left_join(sample_tbl, c("sample_id" = "sample_2"))

plot_raw_counts_hist_pdf(samples_to_plot1, "~/census_sample_counts_gap_min_mean_0.25.pdf")


impute_x_approximate_distribution <- function(df, counts_gap_threshold = 0.25) {
  df |>
    dplyr::mutate(
      inferred_distribution = dplyr::case_when(
        
        # 0) When counts gap between 0 and next min value >= threshold
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          (counts_gap_min_mean >= counts_gap_threshold) ~ "double_log1p",
        
        # 1) Small counts gap
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          (counts_gap_min_mean < counts_gap_threshold) ~ "log1p",
        
        # 2) No negatives, has large values
        !has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw",
        
        # 3) Large values, integer counts
        !has_negative & max_gt_20 & all_integer & !has_floating ~ "raw",
        
        # 4) Has negatives, compressed range
        has_negative & !max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 5) Has negatives and large values
        has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw",
        
        # fallback
        TRUE ~ NA_character_
      )
    )
}

x = samples_to_plot1 |> impute_x_approximate_distribution(0.25)
x = x |> dplyr::rename(inferred_distribution_0.25 = inferred_distribution)
# of those 25 samples, what transform method will be updated after increasing the counts_gap to 0.35, and whats theri distribution?
y =  x |> impute_x_approximate_distribution(0.35) |> dplyr::rename(inferred_distribution_0.35 = inferred_distribution) |>
  impute_x_approximate_distribution(0.45) |> dplyr::rename(inferred_distribution_0.45 = inferred_distribution) |>
  impute_x_approximate_distribution(0.65) |> dplyr::rename(inferred_distribution_0.65 = inferred_distribution) 
y |> dplyr::count(inferred_distribution_0.25, inferred_distribution_0.35, inferred_distribution_0.45, inferred_distribution_0.65)



# use counts_gap_min_mean 0.35
samples_to_plot2 = sample_target_tbl |> left_join(test_sample_summary_df, by = c("sample_id" = "sample_2")) |> 
  filter(counts_gap_min_mean >= 0.35) |>
  group_by(dataset_id) |> 
  slice_head(n=1) |>
  ungroup()

samples_to_plot2 = samples_to_plot2 |>left_join(sample_tbl, c("sample_id" = "sample_2"))

plot_raw_counts_hist_pdf(samples_to_plot2, "~/census_sample_counts_gap_min_mean_0.35.pdf")

# plot samples classified as log1p when increase counts_gap_min_mean from 0.25 to 0.35
samples_to_plot3 = y |> filter(inferred_distribution_0.25 == "double_log1p", inferred_distribution_0.35 == "log1p")
plot_raw_counts_hist_pdf(samples_to_plot3, "~/census_sample_reclassified_log1p_counts_gap_0.25-0.35.pdf")

# plot samples classified as log1p when increase counts_gap_min_mean from 0.25 to 0.45
samples_to_plot4 = y |> filter(inferred_distribution_0.25 == "double_log1p", inferred_distribution_0.35 == "double_log1p", inferred_distribution_0.45 == "log1p")
plot_raw_counts_hist_pdf(samples_to_plot4, "~/census_sample_reclassified_log1p_counts_gap_0.25-0.45.pdf")

# plot samples classified as log1p when increase counts_gap_min_mean from 0.25 to 0.55
samples_to_plot5 = y |> filter(inferred_distribution_0.25 == "double_log1p", inferred_distribution_0.35 == "double_log1p", 
                               inferred_distribution_0.45 == "double_log1p", inferred_distribution_0.65 == "log1p")
plot_raw_counts_hist_pdf(samples_to_plot5, "~/census_sample_reclassified_log1p_counts_gap_0.25-0.65.pdf")


# plot counts_gap_min_mean 0.25 and mode > 1 as double log
impute_x_approximate_distribution2 <- function(df, counts_gap_threshold = 0.25) {
  df |>
    dplyr::mutate(
      inferred_distribution = dplyr::case_when(
        
        # 0) When counts gap between 0 and next min value >= threshold
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          (counts_gap_min_mean >= counts_gap_threshold) & (positive_mode > 1) ~ "double_log1p",
        
        # 1) Small counts gap
        !has_negative & !max_gt_20 & !all_integer & !has_floating &
          !(
            (counts_gap_min_mean >= counts_gap_threshold) &
              (positive_mode > 1)
            
          ) ~ "log1p",
        
        # 2) No negatives, has large values
        !has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw",
        
        # 3) Large values, integer counts
        !has_negative & max_gt_20 & all_integer & !has_floating ~ "raw",
        
        # 4) Has negatives, compressed range
        has_negative & !max_gt_20 & !all_integer & !has_floating ~ "raw_limit_max_to_10",
        
        # 5) Has negatives and large values
        has_negative & max_gt_20 & !all_integer & !has_floating ~ "raw",
        
        # fallback
        TRUE ~ NA_character_
      )
    )
}
x1 = samples_to_plot1 |> impute_x_approximate_distribution2(0.25)
x1 = x1 |> dplyr::rename(inferred_distribution_0.25 = inferred_distribution)
# of those 25 samples, what transform method will be updated after increasing the counts_gap to 0.35, and whats theri distribution?
y1 =  x1 |> impute_x_approximate_distribution(0.35) |> dplyr::rename(inferred_distribution_0.35 = inferred_distribution) |>
  impute_x_approximate_distribution(0.45) |> dplyr::rename(inferred_distribution_0.45 = inferred_distribution) |>
  impute_x_approximate_distribution(0.65) |> dplyr::rename(inferred_distribution_0.65 = inferred_distribution) 
y1 |> dplyr::count(inferred_distribution_0.25, inferred_distribution_0.35, inferred_distribution_0.45, inferred_distribution_0.65)

# columns to track
cols <- c(
  "inferred_distribution_0.25",
  "inferred_distribution_0.35",
  "inferred_distribution_0.45",
  "inferred_distribution_0.65"
)

# reshape and summarise
plot_df <- y1 |>
  dplyr::select(all_of(cols)) |>
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "threshold",
    values_to = "distribution"
  ) |>
  dplyr::mutate(
    threshold = gsub("inferred_distribution_", "", threshold),
    threshold = as.numeric(threshold)
  ) |>
  dplyr::count(threshold, distribution)

plot_df2 <- y1 |>
  dplyr::select(all_of(cols)) |>
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "threshold",
    values_to = "distribution"
  ) |>
  dplyr::mutate(
    threshold = gsub("inferred_distribution_", "", threshold),
    threshold = as.numeric(threshold)
  ) |>
  dplyr::filter(distribution %in% c("double_log1p", "log1p")) |>
  dplyr::count(threshold, distribution)

p <- ggplot(plot_df, aes(x = threshold, y = n, color = distribution)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(
    aes(label = n),
    vjust = -0.6,
    size = 3,
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = c(0.25, 0.35, 0.45, 0.65)) +
  labs(
    x = "counts_gap_threshold",
    y = "Number of samples",
    color = "Distribution",
    title = "Distribution changes across thresholds"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) 

p2 <- ggplot(plot_df2, aes(x = threshold, y = n, color = distribution)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(
    aes(label = n),
    vjust = -0.6,
    size = 3,
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = c(0.25, 0.35, 0.45, 0.65)) +
  labs(
    x = "counts_gap_threshold",
    y = "Number of samples",
    color = "Distribution",
    title = "Distribution changes across thresholds"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
p | p2

ggsave(
  "~/increase_threshold_census_sample_distribution_benchmark.pdf",
  p | p2,
  width = 12,
  height = 6
)

p2

# plot 0.25 double_log1p one sample per dataset
df1 = y1 |> mutate(sample_id = stringr::str_remove(sample_id, ".h5ad")) |>
  left_join(sample_tbl, by = c("sample_id" = "sample_2", "dataset_id")) |> 
  dplyr::filter(inferred_distribution_0.25 == "double_log1p") |>
  group_by(dataset_id) |> 
  slice_head(n=1) |>
  ungroup() |> 
  head(40)
plot_raw_counts_hist_pdf(df1, "~/census_sample_counts_gap_min_mean_gt0.25_mode_gt1.pdf")


df2 = y1 |> mutate(sample_id = stringr::str_remove(sample_id, ".h5ad")) |>
  left_join(sample_tbl, by = c("sample_id" = "sample_2", "dataset_id")) |> 
  dplyr::filter(inferred_distribution_0.25 == "double_log1p", inferred_distribution_0.35 == "double_log1p") |>
  group_by(dataset_id) |> 
  slice_head(n=1) |>
  ungroup() |> 
  head(40)

plot_raw_counts_hist_pdf(df2, "~/census_sample_counts_gap_min_mean_gt0.35_mode_gt1.pdf")


df3 = y1 |> mutate(sample_id = stringr::str_remove(sample_id, ".h5ad")) |>
  left_join(sample_tbl, by = c("sample_id" = "sample_2", "dataset_id")) |> 
  dplyr::filter(inferred_distribution_0.25 == "double_log1p", inferred_distribution_0.35 == "log1p") |>
  group_by(dataset_id) |> 
  slice_head(n=1) |>
  ungroup() |> 
  head(40)
plot_raw_counts_hist_pdf(df3, "~/census_sample_reclassified_log1p_counts_gap_0.25-0.35_mode_gt1.pdf")

# Check those already failed sct samples distribution after applying min_mean 0.25, mode 1

sample_target_tbl <- readRDS("/vast/scratch/users/shen.m/sample_target_tbl.rds")
y1 =  y1 |> mutate(sample_id = stringr::str_remove(sample_id, ".h5ad"))
sample_target_tbl = sample_target_tbl |> left_join(y1, by = c("sample_id"))
sample_target_tbl |> dplyr::count(inferred_distribution_0.25)

