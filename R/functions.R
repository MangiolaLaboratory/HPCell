

#' Cell Type Annotation Transfer
#'
#' @description
#' `annotation_label_transfer` utilizes SingleR for cell-type identification using reference datasets
#' (Blueprint and Monaco Immune data). It can also perform cell type labeling using Azimuth when a reference
#' is provided.
#'
#' @param input_read_RNA_assay SingleCellExperiment object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param reference_azimuth Optional reference data for Azimuth.
#'
#' @return A tibble with cell type annotation data.
#'
#' @importFrom celldex BlueprintEncodeData
#' @importFrom celldex MonacoImmuneData
#' @importFrom scuttle logNormCounts
#' @importFrom SingleR SingleR
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat CreateAssayObject
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom magrittr extract2
#' 
#' @export
annotation_label_transfer <- function(input_read_RNA_assay,
                                      empty_droplets_tbl, 
                                      reference_azimuth = NULL,
                                      assay = NULL
){
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  # SingleR
  sce =
    input_read_RNA_assay |>
    
    # Filter empty
    left_join(empty_droplets_tbl, by = ".cell") |>
    filter(!empty_droplet) |>
    as.SingleCellExperiment() |>    # Add class to the tbl
    logNormCounts()
  
  if(ncol(sce)==1){
    sce = S4Vectors::cbind(sce, sce)
    colnames(sce)[2]= "dummy___"
  }
  blueprint <- celldex::BlueprintEncodeData()
  
  data_annotated =
    
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
  
  data_annotated =
    data_annotated |>
    
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
  if(nrow(data_annotated) == 0){
    
    tibble(.cell = character()) 
    #saveRDS(output_path)
    # output_path
    
  } 
  
  if(nrow(data_annotated) <= 30 | is.null(reference_azimuth)){
    
    # If too little immune cells
    return(data_annotated)
    #saveRDS(output_path)
    #output_path
    
  } else if (!is.null(reference_azimuth)) {
    
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
    RNA_assay = input_read_RNA_assay[[assay]][rownames(input_read_RNA_assay[[assay]]) %in% rownames(reference_azimuth[["SCT"]]),]
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
      SCTransform(assay=assay) |>
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
    modified_data <- data_annotated  |>
      left_join(azimuth_annotation, by = join_by(.cell)	) 
    
    return(modified_data)
  }
}

#' Alive Cell Identification
#'
#' @description
#' `alive_identification` filters out dead cells by analyzing mitochondrial and ribosomal gene expression percentages.
#'
#' @param input_read_RNA_assay SingleCellExperiment object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param annotation_label_transfer_tbl A tibble with annotation label transfer data.
#'
#' @return A tibble identifying alive cells.
#'
#' @importFrom scuttle perCellQCMetrics
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom tidyr unnest
#' @importFrom stringr str_which
#' @importFrom Seurat GetAssayData
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom scater isOutlier
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom purrr map
#' 
#' @export
alive_identification <- function(input_read_RNA_assay,
                                 empty_droplets_tbl,
                                 annotation_label_transfer_tbl,
                                 assay = NULL) {
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  input_read_RNA_assay =
    input_read_RNA_assay |>
    tidySummarizedExperiment::left_join(empty_droplets_tbl, by=".cell") |>
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
  
  # mitochondrion =
  #   input_read_RNA_assay |>
  #   GetAssayData( slot = "counts", assay=assay) |>
  # 
  #   # Join mitochondrion statistics
  #   # Compute per-cell quality control metrics for a count matrix or a SingleCellExperiment
  #   perCellQCMetrics(subsets=list(Mito=which_mito)) |>
  #   as_tibble(rownames = ".cell") |>
  #   select(-sum, -detected) |>
  # 
  #   # Join cell types if annotation_label_transfer_tbl provided
  #   {\(x)
  #     if (inherits(annotation_label_transfer_tbl, "tbl_df")) {
  #       left_join(x, annotation_label_transfer_tbl, by = ".cell") |>
  # 
  #       # Label cells
  #       nest(data = -blueprint_first.labels.fine) |>
  #         mutate(data = map(
  #           data,
  #           ~ .x |>
  #             mutate(high_mitochondrion = isOutlier(subsets_Mito_percent, type="higher")) |>
  # 
  #             # For compatibility
  #             mutate(high_mitochondrion = as.logical(high_mitochondrion))
  #         )) |>
  #         unnest(data)
  #     } else {
  #       x
  #     }
  #   }()
  
  #Extract counts for RNA assay
  rna_counts <- GetAssayData(input_read_RNA_assay, layer = "counts", assay=assay)
  
  # Compute per-cell QC metrics
  qc_metrics <- perCellQCMetrics(rna_counts, subsets=list(Mito=which_mito)) %>%
    as_tibble(rownames = ".cell") %>%
    dplyr::select(-sum, -detected)
  
  # Add cell type labels and determine high mitochondrion content, if annotation_label_transfer_tbl is provided
  if (inherits(annotation_label_transfer_tbl, "tbl_df")) {
    mitochondrion <- qc_metrics %>%
      left_join(annotation_label_transfer_tbl, by = ".cell") %>%
      nest(data = -blueprint_first.labels.fine) %>%
      mutate(data = map(data, ~ .x %>%
                          mutate(high_mitochondrion = isOutlier(subsets_Mito_percent, type="higher"),
                                 high_mitochondrion = as.logical(high_mitochondrion)))) %>%
      unnest(cols = data)
  } else {
    # Determing high mitochondrion content 
    mitochondrion <- qc_metrics %>%
      nest(data = everything()) %>%
      mutate(data = map(data, ~ .x %>%
                          mutate(high_mitochondrion = isOutlier(subsets_Mito_percent, type="higher"),
                                 high_mitochondrion = as.logical(high_mitochondrion)))) %>%
      unnest(cols = data)
  }
  
  
  if (inherits(annotation_label_transfer_tbl, "tbl_df")) {
    
    ribosome =
      input_read_RNA_assay |>
      select(.cell) |>
      #mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay)[,1]) |>
      
      # I HAVE TO DROP UNIQUE, AS SOON AS THE BUG IN SEURAT IS RESOLVED. UNIQUE IS BUG PRONE HERE.
      mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay)) |>
      left_join(annotation_label_transfer_tbl, by = ".cell") |>
      nest(data = -blueprint_first.labels.fine) |>
      mutate(data = map(
        data,
        ~ .x |>
          mutate(high_ribosome = isOutlier(subsets_Ribo_percent, type="higher")) |>
          mutate(high_ribosome = as.logical(high_ribosome)) |>
          as_tibble() |>
          select(.cell, subsets_Ribo_percent, high_ribosome)
      )) |>
      unnest(data)
    
  } else {
    
    # ribosome =
    #   input_read_RNA_assay |>
    #   select(.cell) |>
    #   mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay)[,1])
    ribosome =
      input_read_RNA_assay |>
      dplyr::select(.cell) |>
      #mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay)[,1]) |>
      mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay))|>
      nest(data = everything()) |>
      mutate(data = map(
        data,
        ~ .x |>
          mutate(high_ribosome = isOutlier(subsets_Ribo_percent, type="higher")) |>
          mutate(high_ribosome = as.logical(high_ribosome)) |>
          as_tibble() |>
          dplyr::select(.cell, subsets_Ribo_percent, high_ribosome)
      )) |>
      unnest(data)
  }
  #ribosome_meta <- as.data.frame(ribosome@meta.data)
  
  modified_data <- mitochondrion |>
    left_join(ribosome, by=".cell") |>
    mutate(alive = !high_mitochondrion) # & !high_ribosome ) |>
  
  modified_data
  
}


#' Doublet Identification
#'
#' @description
#' `doublet_identification` applies the scDblFinder algorithm to the filter_empty_droplets dataset. It supports integrating with
#' SingleR annotations if provided and outputs a tibble containing cells with their associated scDblFinder scores.
#'
#' @param input_read_RNA_assay SingleCellExperiment object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param alive_identification_tbl A tibble identifying alive cells.
#' @param annotation_label_transfer_tbl A tibble with annotation label transfer data.
#' @param reference_label_fine Optional reference label for fine-tuning.
#'
#' @return A tibble containing cells with their scDblFinder scores.
#'
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @import scDblFinder
#' @export
doublet_identification <- function(input_read_RNA_assay, 
                                   empty_droplets_tbl, 
                                   alive_identification_tbl, 
                                   annotation_label_transfer_tbl, 
                                   reference_label_fine,
                                   assay = NULL){
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  
  filter_empty_droplets <- input_read_RNA_assay |>
    # Filtering empty
    left_join(empty_droplets_tbl |> select(.cell, empty_droplet), by = ".cell") |>
    filter(!empty_droplet) |>
    
    # Filter dead
    left_join(alive_identification_tbl |> select(.cell, high_mitochondrion, high_ribosome), by = ".cell") |>
    filter(!high_mitochondrion & !high_ribosome) 
  
  # Annotate
  filter_empty_droplets <- filter_empty_droplets |> 
    left_join(annotation_label_transfer_tbl, by = ".cell")|>
    #scDblFinder(clusters = ifelse(reference_label_fine=="none", TRUE, reference_label_fine)) |>
    scDblFinder(clusters = NULL) 
  
  as_tibble(colData(filter_empty_droplets), rownames = ".cell")|> select(.cell, contains("scDblFinder")) 
  
}


#' Cell Cycle Scoring
#'
#' @description
#' Applies cell cycle scoring based on the expression of G2/M and S phase markers. 
#' Returns a tibble containing cell identifiers with their predicted classification 
#' into cell cycle phases: G2M, S, or G1 phase.
#'
#' @param input_read_RNA_assay SingleCellExperiment object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#'
#' @return A tibble with cell identifiers and their cell cycle phase classifications.
#'
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom Seurat CellCycleScoring
#' @export
cell_cycle_scoring <- function(input_read_RNA_assay, 
                               empty_droplets_tbl,
                               assay = NULL){
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
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
      s.features = Seurat::cc.genes$s.genes,
      g2m.features = Seurat::cc.genes$g2m.genes,
      set.ident = FALSE 
    ) |>
    
    as_tibble() |>
    select(.cell,  S.Score, G2M.Score, Phase) 
  
  return(counts)
  
}


#' Non-Batch Variation Removal
#'
#' @description
#' Regresses out variations due to mitochondrial content, ribosomal content, and 
#' cell cycle effects.
#'
#' @param input_path_demultiplexed Path to demultiplexed data.
#' @param input_path_empty_droplets Path to empty droplets data.
#' @param alive_identification_tbl A tibble from alive cell identification.
#' @param cell_cycle_score_tbl A tibble from cell cycle scoring.
#'
#' @return Normalized and adjusted data.
#'
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom Seurat NormalizeData
#' @import sctransform
#' @export
non_batch_variation_removal <- function(input_path_demultiplexed, 
                                        input_path_empty_droplets, 
                                        alive_identification_tbl, 
                                        cell_cycle_score_tbl,
                                        assay = NULL){
  
  # Get assay
  if(is.null(assay)) assay = input_path_demultiplexed@assays |> names() |> extract2(1)
  
  
  counts =
    input_path_demultiplexed |>
    left_join(input_path_empty_droplets, by = ".cell") |>
    filter(!empty_droplet) |>
    
    left_join(
      alive_identification_tbl |>
        select(.cell, subsets_Ribo_percent, subsets_Mito_percent),
      by=".cell"
    ) |>
    
    left_join(
      cell_cycle_score_tbl |>
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
    assay=assay,
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

#' Preprocessing Output
#'
#' @description
#' Incorporates outputs from doublets and dead cell filtering, cell cycle scoring, 
#' and optionally includes annotation label transfer information to generate a 
#' processed dataset ready for downstream analysis.
#'
#' @param tissue Type of tissue.
#' @param non_batch_variation_removal_S Result from non-batch variation removal.
#' @param alive_identification_tbl A tibble from alive cell identification.
#' @param cell_cycle_score_tbl A tibble from cell cycle scoring.
#' @param annotation_label_transfer_tbl A tibble from annotation label transfer.
#' @param doublet_identification_tbl A tibble from doublet identification.
#'
#' @return Processed and filter_empty_droplets dataset.
#'
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @export
preprocessing_output <- function(tissue, 
                                 non_batch_variation_removal_S, 
                                 alive_identification_tbl, 
                                 cell_cycle_score_tbl, 
                                 annotation_label_transfer_tbl, 
                                 doublet_identification_tbl){
  processed_data <- non_batch_variation_removal_S |>
    # Filter dead cells
    left_join(
      alive_identification_tbl |>
        select(.cell, alive, subsets_Mito_percent, subsets_Ribo_percent, high_mitochondrion, high_ribosome),
      by = ".cell"
    ) |>
    filter(alive) |>
    
    # Add cell cycle
    left_join(
      cell_cycle_score_tbl,
      by=".cell"
    ) |>
    
    # Filter doublets
    left_join(doublet_identification_tbl |> select(.cell, scDblFinder.class), by = ".cell") |>
    filter(scDblFinder.class=="singlet") 
  
  if (inherits(annotation_label_transfer_tbl, "tbl_df")){
    processed_data <- processed_data |>
      left_join(annotation_label_transfer_tbl, by = ".cell")
  }
  
  # Filter Red blood cells and platelets
  if (tolower(tissue) == "pbmc" & "predicted.celltype.l2" %in% c(rownames(annotation_label_transfer_tbl), colnames(annotation_label_transfer_tbl))) {
    filtered_data <- filter(processed_data, !predicted.celltype.l2 %in% c("Eryth", "Platelet"))
  } else {
    filtered_data <- processed_data
  }
}


#' Pseudobulk Preprocessing
#'
#' @description
#' Aggregates cells based on sample and cell type annotations, creating pseudobulk samples 
#' for each combination. Handles RNA and ADT assays, ensuring that missing genes are accounted 
#' for and aligns data across multiple samples.
#'
#' @param reference_label_fine Reference label for fine categorization.
#' @param preprocessing_output_S Processed dataset from preprocessing.
#' @param sample_column Column name indicating sample identifiers.
#'
#' @return List containing pseudobulk data aggregated by sample and by both sample and cell type.
#'
#' @import tidySingleCellExperiment
#' @import tidySummarizedExperiment
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom stringr str_remove
#' @importFrom tidyr unite
#' @importFrom tidyr pivot_longer
#' @importFrom tidyseurat aggregate_cells
#' @importFrom tidybulk as_SummarizedExperiment
#' @importFrom S4Vectors cbind
#' @importFrom purrr map
#' @importFrom scater isOutlier
#' @importFrom SummarizedExperiment rowData
#' @export
pseudobulk_preprocessing <- function(reference_label_fine, 
                                     preprocessing_output_S, 
                                     sample_column,
                                     assay = NULL){
  

  
  if (reference_label_fine %in% colnames(preprocessing_output_S[[1]]@meta.data)) {
    
    #sample_column = enquo(sample_column)
    #sample_symbol <- rlang::sym(rlang::quo_get_expr(sample_column))
    pseudobulk =
      preprocessing_output_S |>
      
      # Aggregate
      map(~ { 

        # Get assay
        if(is.null(assay)) {
          assay = .x@assays |> names() |> extract2(1)
          .x = add_RNA_assay(.x, assay)
          
        }
        
        assays = .x@assays |> names() |> intersect(c("RNA", "ADT"))
        
        
        .x |> 
          tidyseurat::aggregate_cells(c(!!as.symbol(sample_column), !!as.symbol(reference_label_fine)), slot = "data", assays=assays) |>
          tidybulk::as_SummarizedExperiment(.sample, .feature, any_of(c("RNA", "ADT"))) |>
          #tidybulk::as_SummarizedExperiment(.sample, .feature, c(RNA)) |>
          
          # Reshape to make RNA and ADT both features
          tidyr::pivot_longer(
            cols = assays,
            names_to = "data_source",
            values_to = "count"
          ) |>
          dplyr::filter(!count |> is.na()) |>
          
          # Some manipulation to get unique feature because RNA and ADT
          # both can have sma name genes
          rename(symbol = .feature) |>
          mutate(data_source = str_remove(data_source, "abundance_")) |>
          unite( ".feature", c(symbol, data_source), remove = FALSE) |>
          
          # Covert
          tidybulk::as_SummarizedExperiment(
            .sample = c( !!as.symbol(sample_column), !!as.symbol(reference_label_fine)),
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
      
      map(~ .x |> dplyr::select(any_of(common_columns)))   %>%
      
      do.call(S4Vectors::cbind, .)
    
    gc()
    
    # Pseuobul aggregated by sample ONLY SAMPLE
    pseudobulk =
      preprocessing_output_S |>
      
      # Aggregate
      map(~ {
        
        # Get assay
        if(is.null(assay)) {
          assay = .x@assays |> names() |> extract2(1)
          .x = add_RNA_assay(.x, assay)
          
        }
        
        assays = .x@assays |> names() |> intersect(c("RNA", "ADT"))
        
        
        .x |> 
          tidyseurat::aggregate_cells(c(!!as.symbol(sample_column)), slot = "data", assays=assays) |>
          tidybulk::as_SummarizedExperiment(.sample, .feature, any_of(c("RNA", "ADT"))) |>

            # Reshape to make RNA and ADT both features
            tidyr::pivot_longer(
              cols = assays,
              names_to = "data_source",
              values_to = "count"
            ) |>
            dplyr::filter(!count |> is.na()) |>
            
            # Some manipulation to get unique feature because RNa and ADT both can have sma name genes
            rename(symbol = .feature) |>
            mutate(data_source = str_remove(data_source, "abundance_")) |>
            unite( ".feature", c(symbol, data_source), remove = FALSE) |>
            
            # Covert
            tidybulk::as_SummarizedExperiment(
              .sample = c( !!as.symbol(sample_column)),
              .transcript = .feature,
              .abundance = count
            )
      }
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
      
      map(~ .x |> dplyr::select(any_of(common_columns)))   %>%
      
      do.call(S4Vectors::cbind, .)
    
    return(list(
      pseudobulk_by_sample = output_path_sample,
      pseudobulk_by_sample_and_cell_type = output_path_sample_cell_type
    ))
    }
  else {return(NULL)}
}






#' Ligand-Receptor Count from Seurat Data
#'
#' @description
#' Calculates ligand-receptor interactions for each cell type in a Seurat object using CellChat.
#'
#' @param counts Seurat object.
#' @param .cell_group Cell group variable.
#' @param assay Name of the assay to use.
#' @param sample_for_plotting Sample name for plotting.
#'
#' @return A list of communication results including interactions and signaling pathways.
#'
#' @importFrom CellChat createCellChat
#' @importFrom CellChat setIdent
#' @importFrom CellChat subsetDB
#' @importFrom CellChat subsetData
#' @importFrom CellChat identifyOverExpressedGenes
#' @importFrom CellChat identifyOverExpressedInteractions
#' @importFrom CellChat projectData
#' @importFrom CellChat filterCommunication
#' @importFrom CellChat aggregateNet
#' @importFrom rlang quo_name
#' @importFrom rlang enquo
#' @importFrom tibble tibble
#' @importFrom purrr map2
#' @importFrom purrr map
#' @importFrom purrr map2_dbl
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


#' Add Dispersion Estimates to SingleCellExperiment Object
#'
#' @description
#' `map_add_dispersion_to_se` function adds dispersion estimates to each feature (gene) in a 
#' SingleCellExperiment object. Dispersion estimates are added based on the abundance measure specified.
#'
#' @param se_df A data frame or list containing SingleCellExperiment objects.
#' @param .col A symbol indicating the column in `se_df` that contains SingleCellExperiment objects.
#' @param abundance (Optional) A character vector specifying the name of the assay to be used 
#'                  for dispersion estimation. If NULL or not provided, the first assay is used.
#'
#' @return The input data frame or list (`se_df`) with the specified `.col` modified to include 
#'         dispersion estimates in each SingleCellExperiment object.
#'
#' @details
#' The function iterates over each SingleCellExperiment object in the specified column of the input data frame 
#' or list. It calculates dispersion estimates for the features (genes) based on the specified abundance assay.
#' The results are joined back to each SingleCellExperiment object. If no abundance assay is specified, 
#' the function defaults to the first assay in each SingleCellExperiment object.
#'
#' @examples
#' # Assuming `se_list` is a list of SingleCellExperiment objects
#' result <- map_add_dispersion_to_se(se_list, .col = se_objects, abundance = "counts")
#'
#' @importFrom magrittr extract2
#' @importFrom edgeR estimateDisp
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom tibble enframe
#' @importFrom purrr map2
#' @export
map_add_dispersion_to_se = function(se_df, .col, abundance = NULL){
  
  .col = enquo(.col)
  
  if(abundance |> length() > 1) stop("HPCell says: for now only one feature abundance measure can be selected")
  
  se_df |>
    
    # This is adding a new column 
    mutate(assay_name = abundance) |> 
    mutate(!!.col := map2(
      !!.col, assay_name,
      ~ {
        
        # If not defined take the first assay
        if(is.null(.y) || .y == "NULL") .y = .x |> assays() |> extract2(1)
        
        counts = .x |> assay(.y)
        
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


#' Test Differential Abundance in SummarizedExperiment Object
#'
#' @description
#' Applies differential abundance testing to each SummarizedExperiment object in a data frame.
#'
#' @param se Data frame containing SummarizedExperiment objects.
#' @param .col Column in the data frame containing the SummarizedExperiment objects.
#' @param formula Formula for the differential abundance test.
#' @param max_rows_for_matrix_multiplication Maximum number of rows for matrix multiplication.
#' @param cores Number of cores to use for computation.
#' @param ... Additional parameters.
#'
#' @return Data frame with test results.
#'
#' @importFrom rlang enquo
#' @export
map_test_differential_abundance = function(
    
  se, .col, .formula, .abundance = NULL, max_rows_for_matrix_multiplication = NULL,
  cores = 1, ...
  
){
  
  .col = enquo(.col)
  .formula = enquo(.formula)
  
  se |> mutate(!!.col := map2(
    !!.col, !!.formula,
    ~ .x |>
      
      # Test
      test_differential_abundance(
        .y, 
        .abundance = !!as.symbol(.abundance),
        method = "glmmSeq_lme4",
        cores = cores,
        max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication,
        .dispersion = dispersion,
        ...
      ),
    ...
    
  ))
  
  
}


#' Split SummarizedExperiment Object by Gene
#'
#' @description
#' Splits each SummarizedExperiment object in a data frame into chunks by gene.
#'
#' @param se_df Data frame containing SummarizedExperiment objects.
#' @param .col Column in the data frame containing the SummarizedExperiment objects.
#' @param .number_of_chunks Number of chunks to split into.
#'
#' @return Data frame with SummarizedExperiment objects split into chunks.
#'
#' @importFrom dplyr n
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest
#' @importFrom dplyr select
#' @importFrom purrr map2
#' @importFrom tidyr nest
#' @importFrom rlang enquo
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

#' @export
map_split_se_by_number_of_genes = function(se_df, .col, chunk_size = 100){
  
  .col = enquo(.col)
  
  se_df |>
    mutate(!!.col := map(
      !!.col,
      ~ {
        total_rows = nrow(.x)
        num_chunks = ceiling(total_rows / chunk_size)
        
        chunks =
          tibble(.feature = rownames(.x)) |>
          mutate(chunk___ = rep(1:num_chunks, each = chunk_size, length.out = nrow(.x)))
        
        # Join chunks
        grouping_factor = chunks |> pull(chunk___) |> as.factor()
        
        .x |> splitRowData(f = grouping_factor)
      }
    )) |>
    unnest(!!.col) |>
    mutate(se_md5 = ids::random_id(n()))
}


#' @importFrom digest digest
#' @importFrom rlang enquo
#'
#' @export
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


