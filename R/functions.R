# Input: seurat, output: nested tibble of variable features
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat VariableFeatures
#' @importFrom glue glue
seurat_to_variable_features_by_cell_type = function(counts, assay, .cell_group = NULL, features_number_per_cell_group = 300){

  .cell_group = enquo(.cell_group)

  # Nest
  counts =
    counts |>
    nest(data = -!!.cell_group) |>

    # Filter more than 10 cells
    filter(map_int(data, ncol) > 100)

  # If I have enough information per cell type
  if(counts |> nrow() |> gt(1)){

    # Per cell type
    counts |>

      # Get feature within each cluster/cell-type
      mutate(feature = map(
        data,
        ~ .x |>
          FindVariableFeatures(nfeatures = features_number_per_cell_group, assay=assay) |>
          VariableFeatures(assay=assay)
      )) |>
      select(-data) |>
      unnest(feature) |>

      # Rename
      rename(group := !!.cell_group)

  }
  else {

    warning(glue("jascap says: you have only one distinct `{quo_name(.cell_group)}`, the per-cell-group variable gene detection will be skipped as it would olverlap with the global detection."))

    tibble()
  }
}

#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat VariableFeatures
seurat_to_variable_features_overall = function(counts, assay, features_number = 300){

  counts |>
    FindVariableFeatures(nfeatures = features_number, assay=assay) |>
    VariableFeatures(assay=assay) |>
    as_tibble() |>
    rename("feature" = "value") |>
    mutate(group= "variable_overall")

}

#' @importFrom rlang enquo
#' @importFrom rlang is_symbolic
#' @import tidyseurat
#' @importFrom Seurat NormalizeData
#' @importFrom stringr str_subset
#'
#' @export
#'
#'
#'
seurat_to_variable_features = function(
    counts,
    assay,
    .sample,
    .cell_group = NULL,
    features_number_independent_of_cell_groups = 300,
    features_number_per_cell_group = 300
){

  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  # If more than one sample balance the size
  if(counts |> distinct(!!.sample) |> nrow() |> gt(1))
    counts =
      counts |>

      # Sample up to a plateau to avoid extreme cell_type bias
      nest(data = -!!.sample) %>%
      mutate(n = map_int(data, ~ ncol(.x))) %>%
      mutate(upper_quantile = quantile(n, 0.75) %>% as.integer()) %>%
      mutate(data = map2(
        data, upper_quantile,
        ~ sample_n(.x, min(ncol(.x), .y), replace = FALSE)
      )) %>%
      filter(n>1) %>%
      unnest(data)

  # Normalise before - https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#normalizing-the-data-1
  counts  = counts |> NormalizeData(assay=assay)

  # Drop TCR, MT and RPL, RPS
  counts = counts[rownames(counts) |> str_subset("^MT|^RPL|^RPS|TRAV|TRBV|TRDV|TRGV", negate = TRUE), ]

  # Variable overall
  variable_df_overall = seurat_to_variable_features_overall(
    counts,
    assay,
    features_number = features_number_independent_of_cell_groups
  )

  # If cell_type_column_for_subsetting == null calculate clusters\
  if(!is_symbolic(.cell_group)){

    counts =
      counts |>
      FindVariableFeatures(nfeatures = number_features_overall, assay=assay)  |>
      NormalizeData(assay=assay) |>
      ScaleData(assay=assay) |>
      RunPCA(assay=assay) |>
      FindNeighbors(dims = 1:20) |>
      FindClusters(resolution = 0.5)

    .cell_group = as.symbol("seurat_clusters")
  }

  variable_df_by_cell_type = seurat_to_variable_features_by_cell_type(
    counts,
    assay,
    .cell_group = !!.cell_group,
    features_number_per_cell_group = features_number_per_cell_group
  )

  variable_df_overall  |>
    bind_rows(variable_df_by_cell_type)

}

#' @export
subset_top_rank_variable_genes_across_batches = function(
    table_across_cell_groups,
    table_within_cell_groups,
    .cell_group,
    .batch,
    features_number_independent_of_cell_groups = 2000,
    features_number_per_cell_group = 300
  ){

  .cell_group = enquo(.cell_group)
  .batch = enquo(.batch)

  batches_with_more_than_10_cell_types =
    table_within_cell_groups |>
    nest(data = -!!.batch) |>
    filter(map_int(data, ~ .x |> distinct(!!.cell_group) |> nrow()) > 10) |>
    unnest(data) |>
    pull(!!.batch)

  # Across cell types
  variable_across_cell_types =
    table_across_cell_groups |>

  # Filter files that have more than 10 cell types
  # because the genes will be overlapped with the cell type specific
  filter(!!.batch %in% batches_with_more_than_10_cell_types)  |>

  count(feature, !!.cell_group) |>
  with_groups(!!.cell_group, ~ .x |> arrange(desc(n)) |> slice(1:features_number_independent_of_cell_groups))


  # Within cell type
  variable_within_cell_types =
    table_within_cell_groups |>
    count(feature, !!.cell_group) |>
    with_groups(!!.cell_group, ~ .x |> arrange(desc(n)) |> slice(1:features_number_per_cell_group))

  # Unique
  bind_rows(

    # Cell type specific
    variable_within_cell_types,

    # Overall
    variable_across_cell_types

  ) |>
    pull(feature) |>
    unique()
}



#' @importFrom CellChat createCellChat setIdent subsetDB subsetData identifyOverExpressedGenes identifyOverExpressedInteractions projectData filterCommunication aggregateNet
#' @importFrom rlang quo_name enquo
#' @importFrom tibble tibble
#' @importFrom purrr map2 map map2_dbl
#'
#' @export
seurat_to_ligand_receptor_count = function(counts, .cell_group, assay, sample_for_plotting = ""){

  .cell_group = enquo(.cell_group)

  # If only one cell, return empty
  if((counts |> distinct(!!.cell_group) |> nrow()) < 2) return(tibble)

  counts_cellchat =
    counts |>

    # Filter
    filter(!is.na(!!.cell_group)) |>

    # Filter cell types with > 10 cells
    add_count(cell_type_harmonised, name = "n_cells") |>
    filter(n_cells>=10) |>


    # Convert from seurat MUST BE LOG-NORMALISED
    createCellChat(group.by = quo_name(.cell_group), assay = assay) |>
    setIdent( ident.use = quo_name(.cell_group))


  communication_results =
    tibble(DB = c("Secreted Signaling", "ECM-Receptor" , "Cell-Cell Contact" )) |>
    mutate(data = list(counts_cellchat)) |>
    mutate(data = map2(
      data, DB,
      ~ {
        print(.y)
        .x@DB <- subsetDB(CellChat::CellChatDB.human, search = .y)

        x = .x |>
          subsetData() |>
          identifyOverExpressedGenes() |>
          identifyOverExpressedInteractions() |>
          projectData(CellChat::PPI.human)

        if(nrow(x@LR$LRsig)==0) return(NA)

        x |>
          computeCommunProb() |>
          filterCommunication() |>
          computeCommunProbPathway() |>
          aggregateNet()

      }
    )) |>

    # Record sample
    mutate(sample = sample_for_plotting) |>

    # Add histogram
    mutate(tot_interactions = map2_dbl(
      data, DB,
      ~  .x |> when(
        !is.na(.) ~ sum(.x@net$count),
        ~ 0
      )
    )) |>

    # Add histogram
    mutate(cell_cell_count = map2_dbl(
      data, DB,
      ~  {
        my_data = .x

        # Return empty if no results
        if(is.na(my_data)) return(tibble(cell_from = character(), cell_to = character(),  weight = numeric()))

        my_data@net$count |>
          as_tibble(rownames = "cell_from") |>
          pivot_longer(-cell_from, names_to = "cell_to", values_to = "count")


      }
    )) |>

    # values_df_for_heatmap
    # Scores for each cell types across all others. How communicative is each cell type
    mutate(cell_vs_all_cells_per_pathway = map2(
      data  , sample,
      ~ when(
        .x,
        !is.na(.x) && length(.x@netP$pathways) > 0 ~
          netAnalysis_computeCentrality(., slot.name = "netP") |>
          cellchat_process_sample_signal(
            pattern = "all", signaling = .x@netP$pathways,
            title = .y, width = 5, height = 6, color.heatmap = "OrRd"
          ),
        ~ tibble(gene = character(), cell_type = character(), value = double())
      )
    ))

  genes = communication_results |> select(cell_vs_all_cells_per_pathway) |> unnest(cell_vs_all_cells_per_pathway) |> distinct(gene) |> pull(gene)

  # Hugh resolution
  communication_results |>
    mutate(cell_vs_cell_per_pathway = map(
      data,
      ~ {
        my_data = .x

        # Return empty if no results
        if(is.na(my_data)) return(tibble(gene = character(),  result = list()))

        tibble(gene = genes) |>
          mutate(result = map(gene, ~ {

            unparsed_result = cellchat_matrix_for_circle(my_data,  layout = "circle", signaling = .x)

            if(!is.null(unparsed_result))
              unparsed_result |>
              as_tibble(rownames = "cell_type_from") |>
              pivot_longer(-cell_type_from, names_to = "cell_type_to", values_to = "score")

            unparsed_result

          }))
      }
    )) |>

    mutate(cell_cell_weight = map(
      data,
      ~ {
        my_data = .x

        # Return empty if no results
        if(is.na(my_data)) return(tibble(cell_from = character(), cell_to = character(),  weight = numeric()))

        my_data@net$weight |>
          as_tibble(rownames = "cell_from") |>
          pivot_longer(-cell_from, names_to = "cell_to", values_to = "weight")


      }
    ))








}

#' @export
map_add_dispersion_to_se = function(se_df){
  
  se_df |> 
    mutate(se = map(
      se,
      ~ {
        counts = .x |> assay("counts")
        
        .x |> 
          left_join(
            
            # Dispersion data frame
            estimateDisp(counts)$tagwise.dispersion |> 
              setNames(rownames(counts)) |>
              enframe(name = ".feature", value = "dispersion")
          )
      }
    ))
  
}

#' @export
map_split_se_by_gene = function(se_df, num_chunks = 10){
  se_df |> 
    mutate(se = map(
      se,
      ~ {
        chunks =
          tibble(.feature = rownames(.x)) |>
          mutate(chunk___ = sample(seq_len(num_chunks), n(), replace = TRUE))
        
        .x |> 
          left_join(chunks) |>
          nest(se_chunk = -chunk___)
      }
    )) |>
    unnest(se) |> 
    select(-chunk___) |>
    mutate(se_md5 = map_chr(se_chunk, digest))
}

splitColData <- function(x, f) {
  # This is by @jma1991 
  # at https://github.com/drisso/SingleCellExperiment/issues/55
  
  i <- split(seq_along(f), f)
  
  v <- vector(mode = "list", length = length(i))
  
  names(v) <- names(i)
  
  for (n in names(i)) { v[[n]] <- x[, i[[n]]] }
  
  return(v)
  
}

splitRowData <- function(x, f) {
  
  i <- split(seq_along(f), f)
  
  v <- vector(mode = "list", length = length(i))
  
  names(v) <- names(i)
  
  for (n in names(i)) { v[[n]] <- x[i[[n]], ] }
  
  return(v)
  
}

#' @importFrom digest digest
#' @importFrom rlang enquo
#' 
#' @export
#' 
map_split_sce_by_gene = function(sce_df, .col, how_many_chunks_base = 10, max_cells_before_split = 4763){
  
  .col = enquo(.col)
  
  sce_df |> 
    mutate(!!.col := map(
      !!.col,
      ~ {
       
        how_many_splits = ceiling(ncol(.x)/max_cells_before_split)*how_many_chunks_base
        
        grouping_factor = sample(seq_len(how_many_splits), size = nrow(.x), replace = TRUE) |> as.factor()
        
        .x |> splitRowData(f = grouping_factor)

      }
    )) |>
    unnest(!!.col) |> 
    mutate(sce_md5 = map_chr(!!.col, digest))
}

#' @export
map_test_differential_abundance = function(
    se, max_rows_for_matrix_multiplication = NULL, 
    cores = 1
){
  
  se |> mutate(se_chunk = map(
    se_chunk,
    ~ .x |>
      
      # Test
      test_differential_abundance(
        ~ dex + (1 | cell),
        method = "glmmSeq_lme4",
        cores = cores, 
        max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication,
        .dispersion = dispersion
      )
    
  ))
  
  
  
  
  
  
}

#' @importFrom readr write_lines
#' @importFrom targets tar_config_get
#' 
#' @export
tar_script_append = function(code, script = targets::tar_config_get("script")){
  substitute(code) |> 
    deparse() |> 
    head(-1) |>
    tail(-1) |> 
    write_lines(script, append = TRUE)
}