# empty_droplet_id
#' @importFrom AnnotationDbi mapIds
#' @importFrom stringr str_subset
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr replace_na
#' @importFrom DropletUtils emptyDrops
#' @importFrom S4Vectors metadata
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom dplyr select
#' @export
empty_droplet_id <- function(input_read_RNA_assay,
                             filter_input){
  significance_threshold = 0.001
  # Genes to exclude
  location <- mapIds(
    EnsDb.Hsapiens.v86,
    keys=rownames(input_read_RNA_assay),
    column="SEQNAME",
    keytype="SYMBOL"
  )
  mitochondrial_genes = which(location=="MT") |> names()
  ribosome_genes = rownames(input_read_RNA_assay) |> str_subset("^RPS|^RPL")
  
  # if ("originalexp" %in% names(input_file@assays)) {
  #   barcode_ranks <- barcodeRanks(input_file@assays$originalexp@counts[!rownames(input_file@assays$originalexp@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE])
  # } else if ("RNA" %in% names(input_file@assays)) {
  #   barcode_ranks <- barcodeRanks(input_file@assays$RNA@counts[!rownames(input_file@assays$RNA@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE])
  # }
  
  # Calculate bar-codes ranks
  barcode_ranks = barcodeRanks(input_read_RNA_assay@assays$RNA@counts[!rownames(input_read_RNA_assay@assays$RNA@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE])
  
  # Set the minimum total RNA per cell for ambient RNA
  if(min(barcode_ranks$total) < 100) { lower = 100 } else {
    lower = quantile(barcode_ranks$total, 0.05)
    
    # write_lines(
    #   glue("{input_path} has supposely empty droplets with a lot of RNAm maybe a lot of ambient RNA? Please investigate"),
    #   file = glue("{dirname(output_path_result)}/warnings_emptyDrops.txt"),
    #   append = T
    # )
  }
  
  # Remove genes from input
  if (
    # If filter_input
    filter_input == "TRUE") {
    barcode_table <- input_read_RNA_assay@assays$RNA@counts[!rownames(input_read_RNA_assay@assays$RNA@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE] |>
      emptyDrops( test.ambient = TRUE, lower=lower) |>
      as_tibble(rownames = ".cell") |>
      mutate(empty_droplet = FDR >= significance_threshold) |>
      replace_na(list(empty_droplet = TRUE))
  }
  else {
    barcode_table <- select(., .cell) |>
      as_tibble() |>
      mutate( empty_droplet = FALSE)
  } 
  
  # barcode ranks
  barcode_table <- barcode_table |>
    left_join(
      barcode_ranks |>
        as_tibble(rownames = ".cell") |>
        mutate(
          knee =  metadata(barcode_ranks)$knee,
          inflection =  metadata(barcode_ranks)$inflection
        )
    )
  
  
  # barcode_table |>  saveRDS(output_path_result)
  
  # # Plot bar-codes ranks
  # plot_barcode_ranks =
  #   barcode_table %>%
  #   ggplot2::ggplot(aes(rank, total)) +
  #   geom_point(aes(color = empty_droplet, size = empty_droplet )) +
  #   geom_line(aes(rank, fitted), color="purple") +
  #   geom_hline(aes(yintercept = knee), color="dodgerblue") +
  #   geom_hline(aes(yintercept = inflection), color="forestgreen") +
  #   scale_x_log10() +
  #   scale_y_log10() +
  #   scale_color_manual(values = c("black", "#e11f28")) +
  #   scale_size_discrete(range = c(0, 2)) +
  #   theme_bw()
  
  # plot_barcode_ranks |> saveRDS(output_path_plot_rds)
  
  # ggsave(
  #   output_path_plot_pdf,
  #   plot = plot_barcode_ranks,
  #   useDingbats=FALSE,
  #   units = c("mm"),
  #   width = 183/2 ,
  #   height = 183/2,
  #   limitsize = FALSE
  # )
  
  barcode_table
  # return(list(barcode_table, plot_barcode_ranks))
}

# annotation_label_transfer
#' @importFrom celldex BlueprintEncodeData
#' @importFrom celldex MonacoImmuneData
#' @importFrom SingleR SingleR
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr select rename
#' @importFrom BiocGenerics ncol nrow
#' @importFrom scuttle logNormCounts
#' @importFrom dplyr left_join filter select
#' @importFrom Seurat CreateSeuratObject CreateAssayObject
#' @export
annotation_label_transfer <- function(input_read_RNA_assay,
                                      reference_azimuth,
                                      empty_droplets_tbl
){
  
  # SingleR
  sce =
    input_read_RNA_assay |>
    
    # Filter empty
    left_join(empty_droplets_tbl, by = ".cell") |>
    filter(!empty_droplet) |>
    as.SingleCellExperiment() |>
    logNormCounts()
  
  if(ncol(sce)==1){
    sce = S4Vectors::cbind(sce, sce)
    colnames(sce)[2]= "dummy___"
  }
  blueprint <- celldex::BlueprintEncodeData()
  
  data_singler =
    
    sce |>
    SingleR(
      ref = blueprint,
      assay.type.test= 1,
      labels = blueprint$label.fine
    )  |>
    as_tibble(rownames=".cell") |>
    nest(blueprint_scores_fine = starts_with("score")) |>
    select(-one_of("delta.next"),- pruned.labels) |>
    rename(blueprint_first.labels.fine = labels) |>
    
    left_join(
      
      sce |>
        SingleR(
          ref = blueprint,
          assay.type.test= 1,
          labels = blueprint$label.main
        )  |>
        as_tibble(rownames=".cell") |>
        nest(blueprint_scores_coarse = starts_with("score")) |>
        select(-one_of("delta.next"),- pruned.labels) |>
        rename( blueprint_first.labels.coarse = labels)
    )
  
  rm(blueprint)
  gc()
  
  
  MonacoImmuneData = celldex::MonacoImmuneData()
  
  data_singler =
    data_singler |>
    
    left_join(
      sce |>
        SingleR(
          ref = MonacoImmuneData,
          assay.type.test= 1,
          labels = MonacoImmuneData$label.fine
        )  |>
        as_tibble(rownames=".cell") |>
        nest(monaco_scores_fine = starts_with("score")) |>
        select(-delta.next,-pruned.labels) |>
        rename(monaco_first.labels.fine = labels)
      
    ) |>
    
    left_join(
      sce |>
        SingleR(
          ref = MonacoImmuneData,
          assay.type.test= 1,
          labels = MonacoImmuneData$label.main
        )  |>
        as_tibble(rownames=".cell") |>
        
        nest(monaco_scores_coarse = starts_with("score")) |>
        select(-delta.next,- pruned.labels) |>
        rename( monaco_first.labels.coarse = labels)
    )  |>
    filter(.cell!="dummy___")
  
  rm(MonacoImmuneData)
  gc()
  
  
  
  rm(sce)
  gc()
  
  # If not immune cells
  if(nrow(data_singler) == 0){
    
    tibble(.cell = character()) 
    #saveRDS(output_path)
    # output_path
    
    
  } else if(nrow(data_singler) <= 30){
    
    # If too little immune cells
    data_singler 
    #saveRDS(output_path)
    #output_path
    
  } else{
    
    print("Start Seurat")
    
    # Load reference PBMC
    # reference_azimuth <- LoadH5Seurat("data//pbmc_multimodal.h5seurat")
    # reference_azimuth |> saveRDS("analysis/annotation_label_transfer/reference_azimuth.rds")
    
    #reference_azimuth = readRDS(reference_azimuth_path)
    
    
    # Reading input
    input_read_RNA_assay =
      
      input_read_RNA_assay |>
      
      # Filter empty
      left_join(empty_droplets_tbl, by = ".cell") |>
      filter(!empty_droplet)
    
    
    # Subset
    RNA_assay = input_read_RNA_assay[["RNA"]][rownames(input_read_RNA_assay[["RNA"]]) %in% rownames(reference_azimuth[["SCT"]]),]
    #RNA_assay <- input_read_RNA_assay@assays$RNA[["counts"]][rownames(input_read_RNA_assay@assays$RNA[["counts"]])%in% rownames(reference_azimuth[["SCT"]]),]
    #ADT_assay = input_read_RNA_assay[["ADT"]][rownames(input_read_RNA_assay[["ADT"]]) %in% rownames(reference_azimuth[["ADT"]]),]
    input_read_RNA_assay <- CreateSeuratObject( counts = RNA_assay)
    
    if("ADT" %in% names(input_read_RNA_assay@assays) ) {
      ADT_assay = input_read_RNA_assay[["ADT"]][rownames(input_read_RNA_assay[["ADT"]]) %in% rownames(reference_azimuth[["ADT"]]),]
      if("ADT" %in% names(input_read_RNA_assay@assays) ) 
        input_read_RNA_assay[["ADT"]] = ADT_assay |> CreateAssayObject()
    }
    
    
    # Normalise RNA
    input_read_RNA_assay =
      input_read_RNA_assay |>
      
      # Normalise RNA - not informed by smartly selected variable genes
      SCTransform(assay="RNA") |>
      ScaleData(assay = "SCT") |>
      RunPCA(assay = "SCT")
    
    if("ADT" %in% names(input_read_RNA_assay@assays) ){
      VariableFeatures(input_read_RNA_assay, assay="ADT") <- rownames(input_read_RNA_assay[["ADT"]])
      input_read_RNA_assay =
        input_read_RNA_assay |>
        NormalizeData(normalization.method = 'CLR', margin = 2, assay="ADT") |>
        ScaleData(assay="ADT") |>
        RunPCA(assay = "ADT", reduction.name = 'apca')
    }
    
    
    # input_file =
    #   input_file |>
    #   FindMultiModalNeighbors(
    #     reduction.list = list("pca", "apca"),
    #     dims.list = list(1:30, 1:18),
    #     modality.weight.name = "RNA.weight"
    #   ) |>
    #   RunUMAP(
    #     nn.name = "weighted.nn",
    #     reduction.name = "wnn.umap",
    #     reduction.key = "wnnUMAP_"
    #   )
    
    # Define common anchors
    anchors <- FindTransferAnchors(
      reference = reference_azimuth,
      query = input_read_RNA_assay,
      normalization.method = "SCT",
      reference.reduction = "spca",
      dims = 1:50
    )
    
    # Mapping
    
    azimuth_annotation =
      tryCatch(
        expr = {
          MapQuery(
            anchorset = anchors,
            query = input_read_RNA_assay,
            reference = reference_azimuth ,
            refdata = list(
              celltype.l1 = "celltype.l1",
              celltype.l2 = "celltype.l2",
              predicted_ADT = "ADT"
            ),
            reference.reduction = "spca",
            reduction.model = "wnn.umap",
            query.dims = 1:2
          )
        },
        error = function(e){
          print(e)
          input_read_RNA_assay |> as_tibble() |> select(.cell)
        }
      ) |>
      as_tibble() |>
      select(.cell, any_of(c("predicted.celltype.l1", "predicted.celltype.l2")), contains("refUMAP"))
    
    # Save
    modified_data <- azimuth_annotation |>
      left_join(	data_singler	, by = join_by(.cell)	) 
    
    modified_data
  }
}
#alive_identification 
#' @importFrom scuttle perCellQCMetrics
#' @importFrom 
#' @importFrom 
#' @importFrom dplyr left_join filter
#' @export
#' 
alive_identification <- function(input_read_RNA_assay,
                                 empty_droplets_tbl,
                                 ann_lbl_trs) {
  input_read_RNA_assay =
    input_read_RNA_assay |>
    left_join(empty_droplets_tbl, by=".cell") |>
    filter(!empty_droplet)
  
  # Returns a named vector of IDs
  # Matches the gene id’s row by row and inserts NA when it can’t find gene names
  location <- mapIds(
    EnsDb.Hsapiens.v86,
    keys=rownames(input_read_RNA_assay),
    column="SEQNAME",
    keytype="SYMBOL"
  )
  
  which_mito = rownames(input_read_RNA_assay) |> str_which("^MT")
  
  mitochondrion =
    input_read_RNA_assay |>
    GetAssayData( slot = "counts", assay="RNA") |>
    
    # Join mitochondrion statistics
    # Compute per-cell quality control metrics for a count matrix or a SingleCellExperiment
    perCellQCMetrics(subsets=list(Mito=which_mito)) |>
    as_tibble(rownames = ".cell") |>
    select(-sum, -detected) |>
    
    # Join cell types
    left_join(ann_lbl_trs, by = ".cell") |>
    
    # Label cells
    nest(data = -blueprint_first.labels.fine) |>
    mutate(data = map(
      data,
      ~ .x |>
        mutate(high_mitochondrion = isOutlier(subsets_Mito_percent, type="higher")) |>
        
        # For compatibility
        mutate(high_mitochondrion = as.logical(high_mitochondrion))
    )) |>
    unnest(data)
  
  
  ribosome =
    
    input_read_RNA_assay |>
    
    select(.cell) |>
    
    # Join mitochondrion statistics
    
    
    mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = "RNA")[,1]  ) |>
    
    # Join cell types
    left_join(ann_lbl_trs, by = ".cell") |>
    
    # Label cells
    nest(data = -blueprint_first.labels.fine) |>
    mutate(data = map(
      data,
      ~ .x |>
        # Label cells
        
        
        mutate(high_ribosome = isOutlier(subsets_Ribo_percent, type="higher")) |>
        
        # For compatibility
        mutate(high_ribosome = as.logical(high_ribosome)) |>
        
        as_tibble() |>
        select(.cell, subsets_Ribo_percent, high_ribosome)
    )) |>
    unnest(data)
  
  modified_data <- mitochondrion |>
    left_join(ribosome, by=".cell") |>
    mutate(alive = !high_mitochondrion) # & !high_ribosome ) |>
  
  modified_data
}

#Doublet identification
#' @importFrom dplyr left_join filter
#' @export
#' 
doublet_identification <- function(input_read_RNA_assay, 
                                   empty_droplets_tbl, 
                                   alive_id, 
                                   ann_lbl_trs, 
                                   reference_label){
  
  input_read_RNA_assay |>
    
    # Filtering empty
    left_join(empty_droplets_tbl |> select(.cell, empty_droplet), by = ".cell") |>
    filter(!empty_droplet) |>
    
    # Filter dead
    left_join(alive_id |> select(.cell, high_mitochondrion, high_ribosome), by = ".cell") |>
    filter(!high_mitochondrion & !high_ribosome) |>
    
    # Annotate
    left_join(ann_lbl_trs, by = ".cell") |>
    
    # Convert
    as.SingleCellExperiment() |>
    
    # Double identification. If no label provided calculate clusters
    scDblFinder(clusters = ifelse(reference_label=="none", TRUE, reference_label)) |>
    
    # Parse
    as_tibble() |>
    select(.cell, contains("scDblFinder")) 
}

#Cell cylce scoring 
#' @importFrom dplyr left_join filter
#' @export
#' 
cell_cycle_scoring <- function(input_read_RNA_assay, 
                               empty_droplets_tbl){
  
  counts =
    input_read_RNA_assay |>
    left_join(empty_droplets_tbl, by = ".cell") |>
    filter(!empty_droplet) |>
    
    # Normalise needed
    NormalizeData() |>
    # Assign cell cycle scores of each cell 
    # Based on its expression of G2/M and S phase markers
    #Stores S and G2/M scores in object meta data along with predicted classification of each cell in either G2M, S or G1 phase
    CellCycleScoring(  
      s.features = cc.genes$s.genes,
      g2m.features = cc.genes$g2m.genes,
      set.ident = FALSE 
    ) |>
    
    as_tibble() |>
    select(.cell,  S.Score, G2M.Score, Phase) 
  
  return(counts)
  
}
#Non_batch_variation_removal
#' @importFrom dplyr left_join filter
#' @export
#' 
non_batch_variation_removal <- function(input_path_demultiplexed, 
                                        input_path_empty_droplets, 
                                        input_path_alive, 
                                        input_cell_cycle_scoring){
  
  counts =
    input_path_demultiplexed |>
    left_join(input_path_empty_droplets, by = ".cell") |>
    filter(!empty_droplet) |>
    
    left_join(
      input_path_alive |>
        select(.cell, subsets_Ribo_percent, subsets_Mito_percent),
      by=".cell"
    ) |>
    
    left_join(
      input_cell_cycle_scoring |>
        select(.cell, G2M.Score),
      by=".cell"
    )
  
  # filter(!high_mitochondrion | !high_ribosome)
  
  # variable_features = readRDS(input_path_merged_variable_genes)
  # 
  # # Set variable features
  # VariableFeatures(counts) = variable_features
  
  # Normalise RNA
  normalized_rna <- SCTransform(
    counts, 
    assay="RNA",
    return.only.var.genes=FALSE,
    residual.features = NULL,
    vars.to.regress = c("subsets_Mito_percent", "subsets_Ribo_percent", "G2M.Score"),
    vst.flavor = "v2",
    scale_factor=2186
  )
  
  # Normalise antibodies
  if ( "ADT" %in% names(normalized_rna@assays)) {
    normalized_data <- normalized_rna %>%
      NormalizeData(normalization.method = 'CLR', margin = 2, assay="ADT") %>%
      select(-subsets_Ribo_percent, -subsets_Mito_percent, -G2M.Score)
  } else { 
    normalized_data <- normalized_rna %>%
      # Drop alive columns
      select(-subsets_Ribo_percent, -subsets_Mito_percent, -G2M.Score)
  }
  
}
# Preprocessing_output
#' @importFrom dplyr left_join 
#' @importFrom dplyr filter
#' @export
#' 
preprocessing_output <- function(tissue, 
                                 input_path_non_batch_variation_removal, 
                                 input_path_alive, 
                                 input_cell_cycle_scoring, 
                                 input_path_annotation_label_transfer, 
                                 input_path_doublet_identification){
  processed_data <- input_path_non_batch_variation_removal |>
    # Filter dead cells
    left_join(
      input_path_alive |>
        select(.cell, alive, subsets_Mito_percent, subsets_Ribo_percent, high_mitochondrion, high_ribosome),
      by = ".cell"
    ) |>
    filter(alive) |>
    
    # Add cell cycle
    left_join(
      input_cell_cycle_scoring,
      by=".cell"
    ) |>
    
    # Filter doublets
    left_join(input_path_doublet_identification |> select(.cell, scDblFinder.class), by = ".cell") |>
    filter(scDblFinder.class=="singlet") |>
    
    # Filter Red blood cells and platelets
    left_join(input_path_annotation_label_transfer, by = ".cell") 
  
  if (tolower(tissue) == "pbmc") {
    filtered_data <- filter(processed_data, !predicted.celltype.l2 %in% c("Eryth", "Platelet"))
  } else {
    filtered_data <- processed_data
  }
}


# Pseudobulk_preprocessing
#' @importFrom dplyr left_join filter
#' @importFrom tidyseurat aggregate_cells
#' @importFrom tidybulk as_SummarizedExperiment
#' @importFrom dplyr select
#' @importFrom S4Vectors cbind
#' @export
#' 
pseudobulk_preprocessing <- function(reference_label_fine, 
                                     input_path_preprocessing_output, sample_column){
  assays = input_path_preprocessing_output[[1]]@assays |> names() |> intersect(c("RNA", "ADT"))
  
  sample_column = enquo(sample_column)
  
  pseudobulk =
    input_path_preprocessing_output |>
    
    # Aggregate
    map(~ { 
      library(rlang)
      .x |> 
        tidyseurat::aggregate_cells(c(!!sample_column, !!as.symbol(reference_label_fine)), slot = "data", assays=assays) |>
        tidybulk::as_SummarizedExperiment(.sample, .feature, any_of(c("RNA", "ADT"))) |>
        #tidybulk::as_SummarizedExperiment(.sample, .feature, c(RNA)) |>
        
        # Reshape to make RNA and ADT both features
        pivot_longer(
          cols = assays,
          names_to = "data_source",
          values_to = "count"
        ) |>
        filter(!count |> is.na()) |>
        
        # Some manipulation to get unique feature because RNA and ADT
        # both can have sma name genes
        rename(symbol = .feature) |>
        mutate(data_source = str_remove(data_source, "abundance_")) |>
        unite( ".feature", c(symbol, data_source), remove = FALSE) |>
        
        # Covert
        tidybulk::as_SummarizedExperiment(
          .sample = c( !!sample_column,  !!as.symbol(reference_label_fine)),
          .transcript = .feature,
          .abundance = count
        )
      
    })
  
  # pseudobulk |> saveRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/PPCG_tumour_microenvironment/PPCG_deconvolution_signatures_single_cell/PPCG_deconvolution_signatures_single_cell_PROCESSED_v3/preprocessing_results/pseudobulk_preprocessing/temp_pseudobulk_preprocessing_sample_output.rds")
  
  # This should not be needed if I create count files weel form the beginning
  # Select only common column
  common_columns =
    pseudobulk |>
    map(~ .x |> as_tibble() |> colnames()) |>
    unlist() |>
    table() %>%
    .[.==max(.)] |>
    names()
  
  all_genes =
    pseudobulk |>
    map(~ .x |> rownames()) |>
    unlist() |>
    unique() |>
    as.character()
  
  # Select and save
  output_path_sample_cell_type <- pseudobulk |>
    
    # Add missing genes
    map(~{
      
      missing_genes = all_genes |> setdiff(rownames(.x))
      
      missing_matrix = matrix(rep(0, length(missing_genes) * ncol(.x)), ncol = ncol(.x))
      
      rownames(missing_matrix) = missing_genes
      colnames(missing_matrix) = colnames(.x)
      
      new_se = SummarizedExperiment(assays = SimpleList(count = missing_matrix))
      colData(new_se) = colData(.x)
      #rowData(new_se) =  DataFrame(symbol = missing_genes, row.names = missing_genes)
      rowData(.x) = NULL
      .x = .x |> rbind(new_se)
      
      .x[all_genes,]
      
    }) |>
    
    map(~ .x |> select(any_of(common_columns)))   %>%
    
    do.call(S4Vectors::cbind, .)
  
  gc()
  
  # ONLY SAMPLE
  pseudobulk =
    input_path_preprocessing_output |>
    
    # Aggregate
    map(~
          .x |>
          tidyseurat::aggregate_cells(c(!!sample_column), slot = "counts", assays=assays) |>
          #tidybulk::as_SummarizedExperiment(!!sample_column, !!as.symbol(reference_label_fine), c(RNA, ADT)) |>
          tidybulk::as_SummarizedExperiment(!!sample_column, .feature, any_of(c("RNA", "ADT"))) |>
          # Reshape to make RNA and ADT both features
          pivot_longer(
            cols = assays,
            names_to = "data_source",
            values_to = "count"
          ) |>
          filter(!count |> is.na()) |>
          
          # Some manipulation to get unique feature because RNa and ADT both can have sma name genes
          rename(symbol = .feature) |>
          mutate(data_source = str_remove(data_source, "abundance_")) |>
          unite( ".feature", c(symbol, data_source), remove = FALSE) |>
          
          # Covert
          tidybulk::as_SummarizedExperiment(
            .sample = c( !!sample_column),
            .transcript = .feature,
            .abundance = count
          )
    )
  
  # This should not be needed if I create count files weel form the beginning
  # Select only common column
  common_columns =
    pseudobulk |>
    map(~ .x |> as_tibble() |> colnames()) |>
    unlist() |>
    table() %>%
    .[.==max(.)] |>
    names()
  
  # Select and save
  output_path_sample <- pseudobulk |>
    # Add missing genes
    map(~{
      
      missing_genes = all_genes |> setdiff(rownames(.x))
      
      missing_matrix = matrix(rep(0, length(missing_genes) * ncol(.x)), ncol = ncol(.x))
      
      rownames(missing_matrix) = missing_genes
      colnames(missing_matrix) = colnames(.x)
      
      new_se = SummarizedExperiment(assay = list(count = missing_matrix))
      colData(new_se) = colData(.x)
      #rowData(new_se) =  DataFrame(symbol = missing_genes, row.names = missing_genes)
      rowData(.x) = NULL
      .x = .x |> rbind(new_se)
      
      .x[all_genes,]
      
    }) |>
    
    map(~ .x |> select(any_of(common_columns)))   %>%
    
    do.call(S4Vectors::cbind, .)
  
  return(list(
    pseudobulk_by_sample = output_path_sample,
    pseudobulk_by_sample_and_cell_type = output_path_sample_cell_type
  ))
}
# # Reference_label_fine
#' @export
#' 
reference_label_fine_id <- function(tissue) {
  return(
    ifelse(tissue == "pbmc", "monaco_first.labels.fine",
    ifelse(tissue == "solid", "blueprint_first.labels.fine",
    ifelse(tissue == "atypical", "none",
    ifelse(tissue == "none", "monaco_first.labels.fine", NA)))))
}
# Reference_label_coarse
#' @export
#' 
reference_label_coarse_id <- function(tissue) {
  return(
    ifelse(tissue == "pbmc", "monaco_first.labels.coarse",
    ifelse(tissue == "solid", "blueprint_first.labels.coarse",
    ifelse(tissue == "atypical", "none",
    ifelse(tissue == "none", "monaco_first.labels.coarse", NA)))))
}

# Add_RNA_assay
#' @export
#'
add_RNA_assay <- function(input_read, RNA_assay_name){
  
  if (RNA_assay_name != "RNA"){
  input_read[["RNA"]] = input_read[[RNA_assay_name]]
  DefaultAssay(object = input_read) <- "RNA"
  input_read[[RNA_assay_name]] = NULL
  }

  # names(input_read@assays)<- names(input_read@assays) |> sapply(function(x) if(x == RNA_assay_name) "RNA" else x)
  input_read
}


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
map_add_dispersion_to_se = function(se_df, .col){
  
  .col = enquo(.col)
  
  se_df |>
    mutate(!!.col := map(
      !!.col,
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

#' @importFrom dplyr n
#' @importFrom purrr map_chr
#' @importFrom purrr map2
#' @importFrom dplyr left_join
#' @importFrom tidyr nest
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#'
#'
#' @export
map_split_se_by_gene = function(se_df, .col, .number_of_chunks){
  
  .col = enquo(.col)
  .number_of_chunks = enquo(.number_of_chunks)
  
  se_df |>
    mutate(!!.col := map2(
      !!.col, !!.number_of_chunks,
      ~ {
        chunks =
          tibble(.feature = rownames(.x)) |>
          mutate(chunk___ = min(1, .y):.y |> sample() |> rep(ceiling(nrow(.x)/max(1, .y))) |> head(nrow(.x)))

        # Join chunks
        grouping_factor = chunks |> pull(chunk___) |> as.factor()

         .x |> splitRowData(f = grouping_factor)
      }
    )) |>
    unnest(!!.col) |>
    mutate(se_md5 = ids::random_id(n()))
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

# library(parallel)
#
# splitRowDataParallel <- function(x, f, numCores = detectCores() - 1) {
#   i <- split(seq_along(f), f)
#
#   v <- mclapply(names(i), function(n) {
#     x[i[[n]], ]
#   }, mc.cores=numCores)
#
#   names(v) <- names(i)
#   return(v)
# }

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
    se, .col, formula, max_rows_for_matrix_multiplication = NULL,
    cores = 1
){
  
  .col = enquo(.col)
  
  se |> mutate(!!.col := map(
    !!.col,
    ~ .x |>
      
      # Test
      test_differential_abundance(
        formula,
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
#' 
tar_script_append = function(code, script = targets::tar_config_get("script")){
  substitute(code) |>
    deparse() |>
    head(-1) |>
    tail(-1) |>
    write_lines(script, append = TRUE)
}
