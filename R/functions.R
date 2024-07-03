## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Cell Type Annotation Transfer
#'
#' @description
#' `annotation_label_transfer` utilizes SingleR for cell-type identification using reference datasets
#' (Blueprint and Monaco Immune data). It can also perform cell type labeling using Azimuth when a reference
#' is provided.
#'
#' @param assay The assay to be used for analysis, specified as a character string.
#' @param input_read_RNA_assay A `SingleCellExperiment` or `Seurat` object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param reference_azimuth Optional reference data for Azimuth.
#' @param assay assay used, default = "RNA" 
#'
#' @return A tibble with cell type annotation data.
#'
#' @importFrom celldex BlueprintEncodeData MonacoImmuneData
#' @importFrom Seurat CreateAssayObject SCTransform CreateSeuratObject SCTransform
#' @importFrom Seurat VariableFeatures FindTransferAnchors MapQuery as.SingleCellExperiment
#' @importFrom scuttle logNormCounts
#' @importFrom SingleR SingleR
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr select join_by rename left_join filter
#' @importFrom magrittr extract2
#' 
#' @export
annotation_label_transfer <- function(input_read_RNA_assay,
                                      empty_droplets_tbl, 
                                      reference_azimuth = NULL,
                                      assay = NULL
){
  # Fix github checks 
  empty_droplet = NULL 
  pruned.labels = NULL 
  delta.next = NULL 
  .cell = NULL 
  
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  # SingleR
  if (inherits(input_read_RNA_assay, "Seurat")) {
    sce =
      input_read_RNA_assay |>
      # Filter empty
      left_join(empty_droplets_tbl, by = ".cell") |>
      dplyr::filter(!empty_droplet) |>
      as.SingleCellExperiment() |>
      logNormCounts()
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")){
    sce =
      input_read_RNA_assay |>
      # Filter empty
      left_join(empty_droplets_tbl, by = ".cell") |>
      dplyr::filter(!empty_droplet) |>
      logNormCounts()
  }

  
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
  
  # Convert SCE to SE to calculate SCT
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    assay(input_read_RNA_assay, assay) <- assay(input_read_RNA_assay, assay) |> as("dgCMatrix")
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL) |> 
      RenameAssays(originalexp = assay)
  }
  
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
    
    #print("Start Seurat")
    
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
    RNA_assay = input_read_RNA_assay[rownames(input_read_RNA_assay[[assay]]) %in% rownames(reference_azimuth[["SCT"]]),][[assay]]
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
      Seurat::VariableFeatures(input_read_RNA_assay, assay="ADT") <- rownames(input_read_RNA_assay[["ADT"]])
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
    anchors <- Seurat::FindTransferAnchors(
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
          Seurat::MapQuery(
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
      left_join(azimuth_annotation, by = dplyr::join_by(.cell)	) 
    
    return(modified_data)
  }
}

#' Alive Cell Identification
#'
#' @description
#' `alive_identification` filters out dead cells by analyzing mitochondrial and ribosomal gene expression percentages.
#'
#' @param assay The assay to be used for analysis, specified as a character string. 
#' @param input_read_RNA_assay A `SingleCellExperiment` or `Seurat` object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param annotation_label_transfer_tbl A tibble with annotation label transfer data.
#' @param assay assay used, default = "RNA" 
#'
#' @return A tibble identifying alive cells.
#'
#' @importFrom scuttle perCellQCMetrics
#' @importFrom AnnotationDbi mapIds
#' @importFrom dplyr left_join filter mutate select
#' @importFrom tidyr unnest
#' @importFrom stringr str_which
#' @importFrom Seurat GetAssayData PercentageFeatureSet
#' @importFrom scater isOutlier
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom purrr map
#' @importFrom magrittr not
#' @importFrom Matrix colSums
#' @importFrom magrittr extract2 not
#' @importFrom SummarizedExperiment assay colData
#' 
#' @export
alive_identification <- function(input_read_RNA_assay,
                                 empty_droplets_tbl,
                                 annotation_label_transfer_tbl = NULL,
                                 annotation_column = NULL,
                                 assay = NULL) {
  
  # Fix GCHECK notes
  empty_droplet = NULL
  detected = NULL
  .cell = NULL 
  high_mitochondrion = NULL 
  blueprint_first.labels.fine = NULL
  
  if(
    !is.null(annotation_column) && 
    !annotation_column %in% colnames(as_tibble(input_read_RNA_assay[1,1]))
  )
    stop("HPCell says: Your `group_by` columns are not present in your data. Please run annotate_cell_type_hpc() to get the cell type annotation that you can use as grouping for the cell-type-specific quality control and removal of dead cells.")
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  input_read_RNA_assay =
    input_read_RNA_assay |>
    left_join(empty_droplets_tbl, by=".cell") |>
    dplyr::filter(!empty_droplet)
  
  # Calculate nFeature_RNA and nCount_RNA if not exist in the data
  nFeature_name <- paste0("nFeature_", assay)
  nCount_name <- paste0("nCount_", assay)
  
  
  if (inherits(input_read_RNA_assay, "Seurat")) {
    counts <- GetAssayData(input_read_RNA_assay, assay = assay, slot = "counts")
    if (!any(str_which(colnames(input_read_RNA_assay[[]]), nFeature_name)) ||
        !any(str_which(colnames(input_read_RNA_assay[[]]), nCount_name))) {
      input_read_RNA_assay[[nFeature_name]] <-
        Matrix::colSums(counts > 0)
      input_read_RNA_assay[[nCount_name]] <-
        Matrix::colSums(counts)
    } else {
      input_read_RNA_assay
    }
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    counts <- assay(input_read_RNA_assay, assay = assay)
    if (!any(str_which(colnames(colData(input_read_RNA_assay)), nFeature_name)) ||
        !any(str_which(colnames(colData(input_read_RNA_assay)), nCount_name))) {
      colData(input_read_RNA_assay)[[nFeature_name]] <- Matrix::colSums(counts > 0)
      colData(input_read_RNA_assay)[[nCount_name]] <- Matrix::colSums(counts)
    }
  }
  
  
  
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
  #       nest(data = -all_of(annotation_column)) |>
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
  if (inherits(input_read_RNA_assay, "Seurat")){
    rna_counts <- GetAssayData(input_read_RNA_assay, layer = "counts", assay=assay)
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    rna_counts <- assay(input_read_RNA_assay, assay=assay)
    assay(input_read_RNA_assay, assay) <- assay(input_read_RNA_assay, assay) |> as("dgCMatrix")
    
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL) |>
      # avoid auto renaming assay name to originalexp after converting
      RenameAssays(originalexp = assay)
  }
  
  # Compute per-cell QC metrics
  qc_metrics <- perCellQCMetrics(rna_counts, subsets=list(Mito=which_mito)) %>%
    as_tibble(rownames = ".cell") %>%
    dplyr::select(-sum, -detected)
  
  # Compute ribosome statistics
  ribosome =
    input_read_RNA_assay |>
    select(.cell) |>
    #mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay)[,1]) |>
    
    # I HAVE TO DROP UNIQUE, AS SOON AS THE BUG IN SEURAT IS RESOLVED. UNIQUE IS BUG PRONE HERE.
    mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay))
  
  # Add cell type labels and determine high mitochondrion content, if annotation_label_transfer_tbl is provided
  if(annotation_column |> is.null() |> not()) {
    
   if (
     inherits(annotation_label_transfer_tbl, "tbl_df") &&
     annotation_column %in% colnames(annotation_label_transfer_tbl)
   ) {
     
     mitochondrion <- qc_metrics %>%
       left_join(annotation_label_transfer_tbl, by = ".cell") 
     
     ribosome =
       ribosome |>
       left_join(annotation_label_transfer_tbl, by = ".cell") 
   }
  
  
    else if (annotation_column %in% colnames(as_tibble(input_read_RNA_assay[1,1]))) {
      
      mitochondrion <- 
        qc_metrics %>%
        left_join(input_read_RNA_assay |> select(.cell, all_of(annotation_column)), by = ".cell") 
      
      
      ribosome = 
        ribosome |>
        left_join(input_read_RNA_assay |> select(.cell, all_of(annotation_column)), by = ".cell") 
    }

  
    mitochondrion = 
      mitochondrion %>%
      nest(data = -all_of(annotation_column)) 
    
    ribosome =
      ribosome |>
      nest(data = -all_of(annotation_column)) 
    
  } else {
    # Determing high mitochondrion content 
    mitochondrion <- qc_metrics %>%
      nest(data = everything()) 
    
    ribosome <- ribosome %>%
      nest(data = everything()) 
  }
  
  mitochondrion = mitochondrion %>%
    mutate(data = map(data, ~ .x %>%
                        mutate(high_mitochondrion = isOutlier(subsets_Mito_percent, type="higher"),
                               high_mitochondrion = as.logical(high_mitochondrion)))) %>%
    unnest(cols = data)

  ribosome = 
    ribosome |> 
    mutate(data = map(
      data,
      ~ .x |>
        mutate(high_ribosome = isOutlier(subsets_Ribo_percent, type="higher")) |>
        mutate(high_ribosome = as.logical(high_ribosome)) |>
        as_tibble() |>
        select(.cell, subsets_Ribo_percent, high_ribosome)
    )) |>
    unnest(data)
  
  # Merge
  mitochondrion |>
    left_join(ribosome, by=".cell") |>
    mutate(alive = !high_mitochondrion) # & !high_ribosome ) |>

}


#' Doublet Identification
#'
#' @description
#' `doublet_identification` applies the scDblFinder algorithm to the filter_empty_droplets dataset. It supports integrating with
#' SingleR annotations if provided and outputs a tibble containing cells with their associated scDblFinder scores.
#' 
#' @param assay The assay to be used for analysis, specified as a character string.
#' @param input_read_RNA_assay A `SingleCellExperiment` or `Seurat` object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param alive_identification_tbl A tibble identifying alive cells.
#' @param annotation_label_transfer_tbl A tibble with annotation label transfer data.
#' @param reference_label_fine Optional reference label for fine-tuning.
#' @param assay Name of the assay to use.
#'
#' @return A tibble containing cells with their scDblFinder scores.
#'
#' @importFrom dplyr left_join filter
#' @importFrom Matrix Matrix 
#' @importFrom SummarizedExperiment colData
#' @importFrom Seurat as.SingleCellExperiment
#' @import scDblFinder
#' @export
doublet_identification <- function(input_read_RNA_assay, 
                                   empty_droplets_tbl, 
                                   alive_identification_tbl, 
                                   #annotation_label_transfer_tbl, 
                                   #reference_label_fine,
                                   assay = NULL){
  
  # Fix GChecks 
  .cell = NULL 
  empty_droplet = NULL 
  high_mitochondrion = NULL 
  high_ribosome = NULL 
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  

  if (inherits(input_read_RNA_assay, "Seurat")) {
    filter_empty_droplets <- input_read_RNA_assay |>
      # Filtering empty
      Seurat::as.SingleCellExperiment() |>
      left_join(empty_droplets_tbl |> select(.cell, empty_droplet), by = ".cell") |>
      filter(!empty_droplet) |>
      
      # Filter dead
      left_join(alive_identification_tbl |> select(.cell, alive), by = ".cell") |>
      filter(alive) 
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    filter_empty_droplets <- input_read_RNA_assay |>
      # Filtering empty
      left_join(empty_droplets_tbl |> select(.cell, empty_droplet), by = ".cell") |>
      filter(!empty_droplet) |>
      
      # Filter dead
      left_join(alive_identification_tbl |> select(.cell, high_mitochondrion, high_ribosome), by = ".cell") |>
      filter(alive) 
  }
  
  # Annotate
  filter_empty_droplets <- filter_empty_droplets |> 
    #left_join(annotation_label_transfer_tbl, by = ".cell")|>
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
#' @param assay The assay to be used for analysis, specified as a character string.
#' @param input_read_RNA_assay A `SingleCellExperiment` or `Seurat` object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param assay Name of the assay to use.
#'
#' @return A tibble with cell identifiers and their cell cycle phase classifications.
#'
#' @importFrom dplyr left_join filter select
#' @importFrom tibble as_tibble
#' @importFrom Seurat CellCycleScoring as.Seurat NormalizeData
#' @export
cell_cycle_scoring <- function(input_read_RNA_assay, 
                               empty_droplets_tbl,
                               assay = NULL){
  #Fix GCHECK
  empty_droplet = NULL 
  .cell = NULL 
  S.Score = NULL 
  G2M.Score = NULL 
  Phase = NULL 
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  # Convert to Seurat in order to perform cell cycle scoring
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    assay(input_read_RNA_assay, assay) <- assay(input_read_RNA_assay, assay) |> as("dgCMatrix")
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL) |>
      RenameAssays(originalexp = assay)
  }
  
  counts <-
    input_read_RNA_assay |>
    left_join(empty_droplets_tbl, by = ".cell") |>
    dplyr::filter(!empty_droplet) |>
    
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
  
  counts
  
}


#' Non-Batch Variation Removal
#'
#' @description
#' Regresses out variations due to mitochondrial content, ribosomal content, and 
#' cell cycle effects.
#'
#' @param input_read_RNA_assay A `SingleCellExperiment` or `Seurat` object containing RNA assay data.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param alive_identification_tbl A tibble from alive cell identification.
#' @param cell_cycle_score_tbl A tibble from cell cycle scoring.
#' @param assay assay used, default = "RNA" 
#'
#' @return Normalized and adjusted data.
#'
#' @importFrom dplyr left_join filter
#' @importFrom Seurat NormalizeData VariableFeatures SCTransform
#' @export
non_batch_variation_removal <- function(input_read_RNA_assay, 
                                        empty_droplets_tbl, 
                                        alive_identification_tbl, 
                                        cell_cycle_score_tbl,
                                        assay = NULL,
                                        factors_to_regress = NULL){
  #Fix GChecks 
  empty_droplet = NULL 
  .cell <- NULL 
  subsets_Ribo_percent <- NULL  
  subsets_Mito_percent <- NULL  
  G2M.Score = NULL 
  
  # Your code for non_batch_variation_removal function here
  
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    assay(input_read_RNA_assay, assay) <- assay(input_read_RNA_assay, assay) |> as("dgCMatrix")
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL) |>
      RenameAssays(originalexp = assay)
  }
  
  counts =
    input_read_RNA_assay |>
    left_join(empty_droplets_tbl, by = ".cell") |>
    filter(!empty_droplet) |>
    
    left_join(
      alive_identification_tbl |>
        select(.cell, any_of(factors_to_regress)),
      by=".cell"
    ) 
  
  if(!is.null(cell_cycle_score_tbl)) 
    counts = counts |>
    
    left_join(
      cell_cycle_score_tbl |>
        select(.cell, any_of(factors_to_regress)),
      by=".cell"
    )
  
  # filter(!high_mitochondrion | !high_ribosome)
  
  # variable_features = readRDS(input_path_merged_variable_genes)
  # 
  # # Set variable features
  # VariableFeatures(counts) = variable_features
  
  # Normalise RNA
  normalized_rna <- Seurat::SCTransform(
    counts, 
    assay=assay,
    return.only.var.genes=FALSE,
    residual.features = NULL,
    vars.to.regress = factors_to_regress,
    vst.flavor = "v2",
    scale_factor=2186
  )
  
  my_assays = "SCT"
  
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    normalized_rna <- normalized_rna |> as.SingleCellExperiment(assay = assay)
  } else if (inherits(input_read_RNA_assay, "Seurat")) {
    normalized_rna
  }

  
  # Normalise antibodies
  if ( "ADT" %in% names(normalized_rna@assays)) {
    normalized_data <- normalized_rna %>%
      NormalizeData(normalization.method = 'CLR', margin = 2, assay="ADT") %>%
      select(-subsets_Ribo_percent, -subsets_Mito_percent, -G2M.Score)
    
    my_assays = my_assays |> c("CLR")
    
  } else { 
    normalized_data <- normalized_rna %>%
      # Drop alive columns
      select(-subsets_Ribo_percent, -subsets_Mito_percent, -G2M.Score)
  }
  
  normalized_data[[my_assays]] 
  
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
#' @import SeuratObject
#' @importFrom SummarizedExperiment left_join
#' @import tidySingleCellExperiment 
#' @import tidyseurat
#' @importFrom magrittr not
#' @export
preprocessing_output <- function(input_read_RNA_assay,
                                 empty_droplets_tbl,
                                 non_batch_variation_removal_S, 
                                 alive_identification_tbl, 
                                 cell_cycle_score_tbl, 
                                 annotation_label_transfer_tbl, 
                                 doublet_identification_tbl){
  #Fix GCHECKS 
  .cell <- NULL
  alive <- NULL
  subsets_Mito_percent <- NULL
  subsets_Ribo_percent <- NULL
  high_mitochondrion <- NULL
  high_ribosome <- NULL
  scDblFinder.class <- NULL
  predicted.celltype.l2 <- NULL
  
  if(empty_droplets_tbl |> is.null() |> not())
    input_read_RNA_assay =
      input_read_RNA_assay |>
      left_join(empty_droplets_tbl, by = ".cell") |>
      filter(!empty_droplet) 
  
  # Add normalisation
  if(!is.null(non_batch_variation_removal_S)){
    if(input_read_RNA_assay |> is("Seurat"))
      input_read_RNA_assay[["SCT"]] = non_batch_variation_removal_S
    else if(input_read_RNA_assay |> is("SingleCellExperiment"))
      assay(spe, "SCT") <- non_batch_variation_removal_S
  }
  
  
  
  input_read_RNA_assay <- input_read_RNA_assay |>
    
    # Filter dead cells
    left_join(
      alive_identification_tbl |>
        select(.cell, any_of(c("alive", "subsets_Mito_percent", "subsets_Ribo_percent", "high_mitochondrion", "high_ribosome"))),
      by = ".cell"
    ) |>
    filter(alive) |>
    
    # Filter doublets
    left_join(doublet_identification_tbl |> select(.cell, scDblFinder.class), by = ".cell") |>
    filter(scDblFinder.class=="singlet") 
  
  # Add cell cycle
  if(cell_cycle_score_tbl |> is.null() |> not())
    input_read_RNA_assay <- input_read_RNA_assay |> 
      left_join(
      cell_cycle_score_tbl,
      by=".cell"
    ) 
  
  # Attach annotation
  if (inherits(annotation_label_transfer_tbl, "tbl_df")){
    input_read_RNA_assay <- input_read_RNA_assay |>
      left_join(annotation_label_transfer_tbl, by = ".cell")
  }
  

  input_read_RNA_assay
  # # Filter Red blood cells and platelets
  # if (tolower(tissue) == "pbmc" & "predicted.celltype.l2" %in% c(rownames(annotation_label_transfer_tbl), colnames(annotation_label_transfer_tbl))) {
  #   filtered_data <- filter(processed_data, !predicted.celltype.l2 %in% c("Eryth", "Platelet"))
  # } else {
  #   filtered_data <- processed_data
  # }
}


#' Create pseudobulk
#'
#' @description
#' Aggregates cells based on sample and cell type annotations, creating pseudobulk samples 
#' for each combination. Handles RNA and ADT assays
#'
#' @param preprocessing_output_S Processed dataset from preprocessing.
#' @param assays A character vector specifying the assays to be included in the 
#' pseudobulk creation process, such as c("RNA", "ADT").
#' @param x A grouping variable used to aggregate cells into pseudobulk samples.
#' This variable should be present in the `preprocessing_output_S` object and 
#' typically represents a factor such as sample ID or condition.
#' @param ... Additional arguments passed to internal functions used within 
#' `create_pseudobulk`. This includes parameters for customization of 
#' aggregation, data transformation, or any other process involved in the 
#' creation of pseudobulk samples.
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

# Create pseudobulk for each sample 
create_pseudobulk <- function(preprocessing_output_S, sample_names ,x ,...) {
  #Fix GChecks 
  .sample = NULL 
  .feature = NULL 
  data_source = NULL 
  symbol = NULL 
  
  #browser()
  # x = enquo(x)
  
  if(preprocessing_output_S |> is("Seurat"))
    assays = Seurat::Assays(preprocessing_output_S)
  else if(preprocessing_output_S |> is("SingleCellExperiment"))
    assays = preprocessing_output_S@assays |> names()
  
  # Aggregate cells
  preprocessing_output_S |> 
    
    # Add sample
    mutate(sample_hpc = sample_names) |> 
    
    # Aggregate
    aggregate_cells(c(sample_hpc, !!x), slot = "data") |>
    as_SummarizedExperiment(.sample, .feature, any_of(c("RNA", "ADT"))) |>
    pivot_longer(cols = assays, names_to = "data_source", values_to = "count") |>
    filter(!count |> is.na()) |>
    
    # Some manipulation to get unique feature because RNA and ADT
    # both can have same name genes
    rename(symbol = .feature) |>
    mutate(data_source = stringr::str_remove(data_source, "abundance_")) |>
    unite(".feature", c(symbol, data_source), remove = FALSE) |>
    
    # Covert
    as_SummarizedExperiment(
      .sample = c(!!x),
      .transcript = .feature,
      .abundance = count
    ) 
}

#' Merge pseudobulk from all samples 
#'
#' @description
#' Merge pseudobulk from all samples. Ensures that missing genes are accounted 
#' for and aligns data across multiple samples.
#' @param create_pseudobulk_sample A list pseudobulk samples generated by `create_pseudobulk`
#' @param assays Default is set of `RNA` 
#' @param x User specified character vector for the column from which we subset the samples for pseudobulk analysis 
#' @param ... Additional arguments 
#' 
#' @importFrom purrr map
#' @importFrom dplyr select
#' @importFrom S4Vectors cbind
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment rowData<-
#' 
#' 
#' 
#' @export
#' 
pseudobulk_merge <- function(create_pseudobulk_sample, ...) {
  # Fix GCHECKS 
  . = NULL 

  # Select only common columns
  common_columns =
    create_pseudobulk_sample |>
    purrr::map(~ .x |> as_tibble() |> colnames()) |>
    unlist() |>
    table() %>%
    .[.==max(.)] |>
    names()
  
  # All genes 
  all_genes =
    create_pseudobulk_sample |>
    purrr::map(~ .x |> rownames()) |>
    unlist() |>
    unique() |>
    as.character()
  
  output_path_sample <- create_pseudobulk_sample |>
    # Add missing genes
    purrr::map(~{
      #browser()
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
    
    purrr::map(~ .x |> dplyr::select(any_of(common_columns)))   %>%
    
    do.call(S4Vectors::cbind, .)
  
  # Return the pseudobulk data for this single sample
  return(output_path_sample)
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
#' @importFrom dplyr distinct add_count
#' @export
seurat_to_ligand_receptor_count = function(counts, .cell_group, assay, sample_for_plotting = ""){
  
  #Fix GChecks 
  cell_type_harmonised <- NULL
  n_cells <- NULL
  DB <- NULL
  cell_vs_all_cells_per_pathway <- NULL
  gene <- NULL
  
  # Your code for seurat_to_ligand_receptor_count function here
  
  
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
# @examples
# # Assuming `se_list` is a list of SingleCellExperiment objects
# result <- map_add_dispersion_to_se(se_list, .col = se_objects, abundance = "counts")
#'
#' @importFrom magrittr extract2
#' @importFrom edgeR estimateDisp
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom tibble enframe
#' @importFrom purrr map2
#' @import dplyr 
#' @importFrom data.table :=
#' @export
map_add_dispersion_to_se = function(se_df, .col, abundance = NULL){
  
  # Fix GitChecks 
  assay_name = NULL 
  
  
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
#' @param .formula Formula for the differential abundance test.
#' @param .abundance Optional; a symbol indicating the column of the datasets that specifies the abundance measure. If NULL, a default is used.
#' @param max_rows_for_matrix_multiplication Maximum number of rows for matrix multiplication.
#' @param cores Number of cores to use for computation.
#' @param ... Additional parameters.
#'
#' @return Data frame with test results.
#'
#' @importFrom rlang enquo
#' @import dplyr 
#' @export
map_test_differential_abundance = function(
    
  se, .col, .formula, .abundance = NULL, max_rows_for_matrix_multiplication = NULL,
  cores = 1, ...
  
){
  
  .col = enquo(.col)
  .formula = enquo(.formula)
  
  se |> mutate(!!.col := map2(
    !!.col, !!.formula,
    ~ {
      
      if(ncol(.x) > 2000) method = "glmmseq_glmmTMB"
      else method = "glmmSeq_lme4"
      

      # Test
      test_differential_abundance(
        .x,
        .y, 
        .abundance = !!as.symbol(.abundance),
        method = method,
        cores = cores,
        max_rows_for_matrix_multiplication = max_rows_for_matrix_multiplication,
        .dispersion = dispersion,
        ...
      )
      },
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
#' @import dplyr 
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest
#' @importFrom dplyr select
#' @importFrom purrr map2
#' @importFrom tidyr nest
#' @importFrom rlang enquo
#' @importFrom ids random_id
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

#' map_split_se_by_number_of_genes
#' 
#' @description
#' Splits the singlecellexperiment dataframe into multiple chunks based on gene count
#' 
#' 
#' @importFrom rlang enquo 
#' @importFrom dplyr tibble
#' @importFrom purrr map
#' @importFrom tibble tibble
#' @importFrom ids random_id
#' @param se_df A `SingleCellExperiment` DataFrame that contains gene expression or other related data.
#' @param .col A symbol or string indicating the column in `se_df` which should be dynamically split into multiple chunks based on gene count.
#' @param chunk_size The maximum number of genes each chunk should contain; defaults to 100.
#'
#' @return Returns the input `SingleCellExperiment` DataFrame with an additional column `se_md5` containing unique random IDs for each chunk, and with the data split into chunks according to the specified number of genes.
#' 
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



#' map_split_sce_by_gene Split SingleCellExperiment by Gene 
#' 
#' @importFrom digest digest
#' @importFrom rlang enquo
#' @importFrom purrr map_chr
#' @description
#' Splits a SingleCellExperiment object into multiple chunks based on the number of cells.
#' This function dynamically partitions a SingleCellExperiment object into multiple chunks based on the number of cells per gene across the specified column. It computes the number of splits by dividing the total number of cells by a maximum threshold and multiplying the result by a base number of chunks. This approach allows handling of large datasets by reducing the complexity in each chunk, making it feasible to perform detailed analyses or computational tasks on subsets of data efficiently.
#' The function also assigns a unique MD5 hash to each chunk as an identifier, facilitating tracking and referencing of data subsets in subsequent analyses.
#' 
#' 
#' 
#' @param sce_df A dataframe (or tibble) where one of the columns contains SingleCellExperiment objects 
#' @param .col A symbol or string indicating the column in `sce_df` which should be dynamically split into multiple chunks
#' @param how_many_chunks_base A base number of chunks to divide the data into, adjusted by the actual size of the data in each group.
#' @param max_cells_before_split The maximum number of cells a single chunk can have before it is split into another chunk.
#' 
#' @return Returns the input SingleCellExperiment DataFrame with an additional column `sce_md5` containing MD5 hashes of the chunks, and with the data split according to the specified parameters.
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
    mutate(sce_md5 = purrr::map_chr(!!.col, digest))
}

#' Find variable genes 
#' 
#' @param input_seurat Single Seurat object (Input data)
#' @param empty_droplet Single dataframe containing empty droplet filtering information 
#' @return A vector of variable gene names
#' 
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat FindVariableFeatures
#' 
#' @export
# Find set of variable genes 
find_variable_genes <- function(input_seurat, empty_droplet){
  
  # Set the assay of choice
  assay_of_choice = input_seurat@assays |> names() |> extract2(1)
  
  # Ensure "HTO" and "ADT" assays are removed if present
  if("HTO" %in% names(input_seurat@assays)) input_seurat[["HTO"]] = NULL
  if("ADT" %in% names(input_seurat@assays)) input_seurat[["ADT"]] = NULL
  
  # Filter out empty droplets
  seu<- dplyr::left_join(input_seurat, empty_droplet) |>
    dplyr::filter(!empty_droplet)
  
  # Update Seurat object meta.data after filtering
  # input_seurat@meta.data <- seu
  
  # Scale data
  input_seurat <- ScaleData(seu, assay=assay_of_choice, return.only.var.genes=FALSE)
  
  # Find and retrieve variable features
  input_seurat <- Seurat::FindVariableFeatures(input_seurat, assay=assay_of_choice, nfeatures = 500)
  my_variable_genes <- Seurat::VariableFeatures(input_seurat, assay=assay_of_choice)
  
  return(my_variable_genes)
}

#' Harmonize cell type annotations based on consensus
#'
#' This function harmonizes cell type annotations by matching them with a reference annotation
#' and applying specific rules for non-immune cell types.
#'
#' @param single_cell_data A data frame containing single-cell data with cell type annotations.
#' @param .sample_column The column name specifying sample information.
#' @param .cell_type The column name for the cell type annotations.
#' @param .azimuth The column name for Azimuth annotations.
#' @param .blueprint The column name for Blueprint annotations.
#' @param .monaco The column name for Monaco annotations.
#'
#' @return A data frame with harmonized cell type annotations.
#'
#'
#' @importFrom dplyr across
#' @importFrom readr read_csv
#' @importFrom dplyr bind_rows
#' @importFrom dplyr join_by
#' @importFrom data.table :=
#'
#'
#' @export
annotation_consensus = function(single_cell_data, .sample_column, .cell_type, .azimuth, .blueprint, .monaco){
  # Fix GITCHECK notes
  .sample = NULL 
  cell_type = NULL 
  cell_annotation_azimuth_l2 = NULL 
  cell_annotation_blueprint_singler = NULL 
  cell_annotation_monaco_singler = NULL 
  .cell = NULL 
  cell_type_harmonised = NULL 
  confidence_class = NULL
  
  
  # Fix GCHECK notes
  .cell = NULL
  cell_type_harmonised = NULL
  confidence_class = NULL 
  cell_annotation_azimuth_l2 = NULL 
  cell_annotation_blueprint_singler = NULL
  cell_annotation_monaco_singler = NULL
  cell_type = NULL 
  .sample = NULL 
  .sample_column = enquo(.sample_column)
  .azimuth = enquo(.azimuth)
  .blueprint = enquo(.blueprint)
  .monaco = enquo(.monaco)
  .cell_type = enquo(.cell_type)
  
  # reference_annotation =
  #   CuratedAtlasQueryR::get_metadata() |> 
  #   filter(cell_type_harmonised!="immune_unclassified" | is.na(cell_type_harmonised)) |> 
  #   select(cell_type,
  #          cell_type_harmonised, 
  #          cell_annotation_azimuth_l2, 
  #          cell_annotation_blueprint_singler, 
  #          cell_annotation_monaco_singler,
  #          confidence_class
  #   ) |> 
  #   as_tibble() |> 
  #   mutate(cell_type_clean = cell_type |> clean_cell_types()) |> 
  #   HPCell::clean_cell_types_deeper() |> 
  #   select(-cell_type) |> 
  #   
  #   count(cell_type_harmonised, cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler, confidence_class, cell_type_clean) |>
  #   with_groups(c(cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler, cell_type_clean), ~ .x |> arrange(desc(n)) |> slice(1) )
  # 
  # reference_annotation |> saveRDS("reference_annotation_16_jan_2024.rds")
  
  reference_annotation = readRDS("reference_annotation_16_jan_2024.rds")
  
  annotation=
    single_cell_data |>
    rename(
      .sample := !!.sample_column,
      cell_type := !!.cell_type,
      cell_annotation_azimuth_l2 := !!.azimuth,
      cell_annotation_blueprint_singler := !!.blueprint,
      cell_annotation_monaco_singler := !!.monaco
    ) |>
    select(.cell, .sample, cell_type, cell_annotation_azimuth_l2,cell_annotation_blueprint_singler, cell_annotation_monaco_singler) |> 
    mutate(across(c(cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler),	tolower	)) |>
    mutate(across(c(cell_annotation_azimuth_l2, cell_annotation_blueprint_singler, cell_annotation_monaco_singler),	clean_cell_types	)) |>
    
    is_strong_evidence(cell_annotation_azimuth_l2, cell_annotation_blueprint_singler) |> 
    
    # Clean cell types
    mutate(cell_type_clean = cell_type |> clean_cell_types()) |> 
    left_join(read_csv("~/PostDoc/CuratedAtlasQueryR/dev/metadata_cell_type.csv"),  by = "cell_type") |> 
    clean_cell_types_deeper() |>
    
    # Reference annotation link
    left_join(reference_annotation ) 
  
  annotation_connie_non_immune = 
    annotation |> 
    filter(cell_type_harmonised |> is.na()) |> 
    
    harmonise_names_non_immune() |> 
    
    # Fix some gaps in the original code
    mutate(cell_type_harmonised = case_when(
      cell_type |> tolower() |> str_detect("endothelial") ~ "endothelial_cell",
      cell_type |> tolower() |> str_detect("enodothelial") ~ "endothelial_cell",
      cell_type |> tolower() |> str_detect("epithelial") ~ "epithelial_cell",
      cell_type |> tolower() |> str_detect("fibroblast") ~ "fibroblast",
      TRUE ~ cell_type
    )) |>
    
    mutate(confidence_class = 1)
  
  single_cell_data |> 
    left_join(
      annotation |> 
        filter(!cell_type_harmonised |> is.na()) |> 
        bind_rows(annotation_connie_non_immune) |>
        select(.cell, .sample, cell_type_harmonised, confidence_class),
      by = join_by(.cell, !!.sample_column == .sample)
    )
  
}



