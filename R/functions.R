## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Identify Empty Droplets in Single-Cell RNA-seq Data
#'
#' @description
#' `empty_droplet_id` distinguishes between empty and non-empty droplets using the DropletUtils package.
#' It excludes mitochondrial and ribosomal genes, calculates barcode ranks, and optionally filters input data
#' based on these criteria. The function returns a tibble containing log probabilities, FDR, and a classification
#' indicating whether cells are empty droplets.
#'
#' @param input_read_RNA_assay SingleCellExperiment or Seurat object containing RNA assay data.
#' @param filter_empty_droplets Logical value indicating whether to filter the input data.
#'
#' @return A tibble with columns: logProb, FDR, empty_droplet (classification of droplets).
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stringr str_subset
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr replace_na
#' @importFrom DropletUtils emptyDrops barcodeRanks
#' @importFrom S4Vectors metadata
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom biomaRt useMart getBM
#' 
#' @export
empty_droplet_id <- function(input_read_RNA_assay,
                             total_RNA_count_check  = -Inf,
                             assay = NULL,
                             feature_nomenclature,
                             sample_name){
  
  if(ncol(input_read_RNA_assay) == 0) {
    warning("HPCell says: the sample ", sample_name, " has no cells and returned NULL")
    return(NULL)
  }
  
  if(input_read_RNA_assay |> is.null()) {
    warning("HPCell says: the sample ", sample_name, " is NULL and returned NULL")
    return(NULL)
  }
  

  
  #Fix GChecks 
  FDR = NULL 
  .cell = NULL 
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  # Get counts
  if (inherits(input_read_RNA_assay, "Seurat")) {
    counts <- GetAssayData(input_read_RNA_assay, assay, slot = "counts")
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    counts <- assay(input_read_RNA_assay, assay)
  }
  
  
  significance_threshold = 0.001
  
  # Genes to exclude
  if (feature_nomenclature == "symbol") {
    location <- mapIds(
      EnsDb.Hsapiens.v86,
      keys=rownames(input_read_RNA_assay),
      column="SEQNAME",
      keytype="SYMBOL"
    )
    mitochondrial_genes = which(location=="MT") |> names()
    ribosome_genes = rownames(input_read_RNA_assay) |> str_subset("^RPS|^RPL")
    
  } else if (feature_nomenclature == "ensembl") {
    # all_genes are saved in data/all_genes.rda to avoid recursively accessing biomaRt backend for potential timeout error
    data(ensembl_genes_biomart)
    all_mitochondrial_genes <- ensembl_genes_biomart[grep("MT", ensembl_genes_biomart$chromosome_name), ] 
    all_ribosome_genes <- ensembl_genes_biomart[grep("^(RPL|RPS)", ensembl_genes_biomart$external_gene_name), ]
    mitochondrial_genes <- all_mitochondrial_genes |> 
      filter(ensembl_gene_id %in% rownames(input_read_RNA_assay)) |> pull(ensembl_gene_id)
    ribosome_genes <- all_ribosome_genes |> 
      filter(ensembl_gene_id %in% rownames(input_read_RNA_assay)) |> pull(ensembl_gene_id)
  }
  
  
  # if ("originalexp" %in% names(input_file@assays)) {
  #   barcode_ranks <- barcodeRanks(input_file@assays$originalexp@counts[!rownames(input_file@assays$originalexp@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE])
  # } else if ("RNA" %in% names(input_file@assays)) {
  #   barcode_ranks <- barcodeRanks(input_file@assays$RNA@counts[!rownames(input_file@assays$RNA@counts) %in% c(mitochondrial_genes, ribosome_genes),, drop=FALSE])
  # }
  
  
  filtered_counts <- counts[!(rownames(counts) %in% c(mitochondrial_genes, ribosome_genes)),, drop=FALSE ]
  
  n_expressed_genes_non_zero = (filtered_counts > 0) |> colSums()
  
  filter_empty_droplets = n_expressed_genes_non_zero |> min() < 200
  
  if(!filter_empty_droplets)
    return(  input_read_RNA_assay |> 
               as_tibble() |> 
               select(.cell) |>
               mutate( empty_droplet = FALSE))
  
  quantile_expressed_genes = n_expressed_genes_non_zero |> quantile(0.05)
  
  # Check if empty droplets have been identified
  # nFeature_name <- paste0("nFeature_", assay)
  
  #if (any(input_read_RNA_assay[[nFeature_name]] < total_RNA_count_check)) {
  # filter_empty_droplets <- "TRUE"
  # }
  # else {
  #   filter_empty_droplets <- "FALSE"
  # }
  
  # Attempt to run emptyDrops() and handle potential errors
  tryCatch({
    # WE CANNOT USE AMBIENT PARAMETER BECAUSE WITH PERCULIAR DATASETS WITH 
    # cells with more zeros have also more total RNA counts dysfunction stalls 
    # for example for this sample
    #.cell                                                   dataset_id                           sample_id 
    #AAACCCAAGCTAATCC___eec804b9-2ae5-44f0-a1b5-d721e21257de eec804b9-2ae5-44f0-a1b5-d721e21257de 485c0dac47c6bd0b91fd3ae9d7de7385
    
    emptyDrops(filtered_counts) |> 
      as_tibble(rownames = ".cell") |>
      mutate(empty_droplet = FDR >= significance_threshold) |>
      replace_na(list(empty_droplet = TRUE)) |> 
      mutate(filter_empty_method = "emptyDrops")
    
  }, error = function(e) {
    # Check if the error message matches the specific error
    if (grepl("no counts available to estimate the ambient profile", e$message)) {
      # You can also print a message if you like
      message("Error encountered: ", e$message)
      message("Setting do_filter to FALSE.")
      
      # Return NULL or an empty object as appropriate
      input_read_RNA_assay |> 
        as_tibble() |> 
        select(.cell) |>
        mutate( empty_droplet = n_expressed_genes_non_zero < 200)  |> 
        mutate(filter_empty_method = "expressed_genes_more_than_200")
    } else {
      # For other errors, re-throw the error
      stop(e)
    }
  })
  
  
  # barcode ranks
  # Calculate bar-codes ranks
  # barcode_ranks <- barcodeRanks(filtered_counts)
  # 
  # barcode_table <- barcode_table |>
  #   left_join(
  #     barcode_ranks |>
  #       as_tibble(rownames = ".cell") |>
  #       mutate(
  #         knee =  metadata(barcode_ranks)$knee,
  #         inflection =  metadata(barcode_ranks)$inflection
  #       )
  #   )
  
  
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
  
}

#' Identify Empty Droplets in Single-Cell RNA-seq Data
#'
#' @description
#' `empty_droplet_threshold` distinguishes between empty and non-empty droplets by threshold. 
#' It excludes mitochondrial and ribosomal genes, and filters input data
#' based on defined values of `nCount_RNA` and `nFeature_RNA`
#' The function returns a tibble containing RNA count, RNA feature count indicating whether cells are empty droplets.
#'
#' @param input_read_RNA_assay SingleCellExperiment or Seurat object containing RNA assay data.
#' @param filter_empty_droplets Logical value indicating whether to filter the input data.
#' @param RNA_feature_threshold An optional integer for the number of feature count. Default is 200
#'
#' @return A tibble with columns: Cell, nFeature_RNA, empty_droplet (classification of droplets).
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom stringr str_subset
#' @importFrom dplyr left_join mutate
#' @importFrom tidyr replace_na
#' @importFrom DropletUtils emptyDrops barcodeRanks
#' @importFrom S4Vectors metadata
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom biomaRt useMart getBM
#' 
#' @export
empty_droplet_threshold<- function(input_read_RNA_assay,
                                   total_RNA_count_check  = -Inf,
                                   assay = NULL,
                                   feature_nomenclature,
                                   RNA_feature_threshold = 200,
                                   sample_name){
  
  
  if(ncol(input_read_RNA_assay) == 0) {
    warning("HPCell says: the sample ", sample_name, " has no cells and returned NULL")
    return(NULL)
  }
  
  if(input_read_RNA_assay |> is.null()) {
    warning("HPCell says: the sample ", sample_name, " is NULL and returned NULL")
    return(NULL)
  }
  
  #Fix GChecks 
  FDR = NULL 
  .cell = NULL 
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  # Check if empty droplets have been identified
  nFeature_name <- paste0("nFeature_", assay)
  
  filter_empty_droplets <- "TRUE"
  
  significance_threshold = 0.001
  # Genes to exclude
  if (feature_nomenclature == "symbol") {
    location <- mapIds(
      EnsDb.Hsapiens.v86,
      keys=rownames(input_read_RNA_assay),
      column="SEQNAME",
      keytype="SYMBOL"
    )
    mitochondrial_genes = which(location=="MT") |> names()
    ribosome_genes = rownames(input_read_RNA_assay) |> str_subset("^RPS|^RPL")
    
  } else if (feature_nomenclature == "ensembl") {
    # all_genes are saved in data/all_genes.rda to avoid recursively accessing biomaRt backend for potential timeout error
    data(ensembl_genes_biomart)
    all_mitochondrial_genes <- ensembl_genes_biomart[grep("MT", ensembl_genes_biomart$chromosome_name), ] 
    all_ribosome_genes <- ensembl_genes_biomart[grep("^(RPL|RPS)", ensembl_genes_biomart$external_gene_name), ]
    
    mitochondrial_genes <- all_mitochondrial_genes |> 
      filter(ensembl_gene_id %in% rownames(input_read_RNA_assay)) |> pull(ensembl_gene_id)
    ribosome_genes <- all_ribosome_genes |> 
      filter(ensembl_gene_id %in% rownames(input_read_RNA_assay)) |> pull(ensembl_gene_id)
  }
  
  # Get counts
  if (inherits(input_read_RNA_assay, "Seurat")) {
    counts <- GetAssayData(input_read_RNA_assay, assay, slot = "counts")
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    counts <- assay(input_read_RNA_assay, assay)
  }
  filtered_counts <- counts[!(rownames(counts) %in% c(mitochondrial_genes, ribosome_genes)),, drop=FALSE ]
  
  # filter based on nCount_RNA and nFeature_RNA
  result <- colSums(filtered_counts > 0 ) |> enframe(name = ".cell", value = "nFeature_RNA") |> 
    #left_join(colSums(filtered_counts) |> enframe(name = ".cell", value = "nCount_RNA"), by = ".cell") |>
    mutate(empty_droplet = nFeature_RNA < RNA_feature_threshold)
  
  # # Discard samples with nFeature_RNA density mode < threshold, avoid potential downstream error
  # density_est = result |> pull(nFeature_RNA) |> density()
  # density_value = density_est$x[which.max(density_est$y)]
  # if (density_value < RNA_feature_threshold) return(NULL)
  
  result
}

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
#' @importFrom celldex BlueprintEncodeData
#' @importFrom celldex MonacoImmuneData
#' 
#' # Seurat
#' @importFrom Seurat CreateAssayObject
#' @importFrom Seurat SCTransform
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat FindTransferAnchors
#' @importFrom Seurat MapQuery
#' @importFrom Seurat as.SingleCellExperiment
#' @import Seurat
#' 
#' @importFrom scuttle logNormCounts
#' @importFrom SingleR SingleR
#' @importFrom tibble as_tibble
#' @importFrom tibble tibble
#' @importFrom dplyr select
#' @importFrom dplyr join_by
#' @importFrom dplyr rename
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom magrittr extract2
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assay<-
#' @importFrom Azimuth RunAzimuth
#' @importFrom stringr str_detect
#' @importFrom tidyr nest
#' @importFrom S4Vectors cbind
#' @importFrom S4Vectors Assays
#' @importFrom S4Vectors RenameAssays
#' 
#' @export
annotation_label_transfer <- function(input_read_RNA_assay,
                                      empty_droplets_tbl = NULL, 
                                      reference_azimuth = NULL,
                                      assay = NULL,
                                      feature_nomenclature,
                                      sample_name
){
  # Fix github checks 
  empty_droplet = NULL 
  pruned.labels = NULL 
  delta.next = NULL 
  .cell = NULL 
  
  if(ncol(input_read_RNA_assay) == 0) {
    warning("HPCell says: the sample ", sample_name, " has no cells and returned NULL")
    return(NULL)
  }
  
  if(input_read_RNA_assay |> is.null()) {
    warning("HPCell says: the sample ", sample_name, " is NULL and returned NULL")
    return(NULL)
  }
  
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  # TEMPORARY FOR SOME REASON THE MIN COUNTS IS NOT 0 FOR SOME SAMPLES
  input_read_RNA_assay = check_if_assay_minimum_count_is_zero_and_correct_TEMPORARY(input_read_RNA_assay, assay)
  
  
  if (!is.null(empty_droplets_tbl)) {
    input_read_RNA_assay =
      input_read_RNA_assay |>
      left_join(empty_droplets_tbl, by=".cell") |>
      dplyr::filter(!empty_droplet)
  }
  
  # SingleR
  if (inherits(input_read_RNA_assay, "Seurat")) {
    sce =
      input_read_RNA_assay |>
      as.SingleCellExperiment() |>
      logNormCounts(assay.type = assay)
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")){
    sce =
      # Filter empty
      input_read_RNA_assay|>
      logNormCounts(assay.type = assay)
  }
  
  # This because an error is num cell = 1
  if(ncol(input_read_RNA_assay)==1){
    input_read_RNA_assay = S4Vectors::cbind(input_read_RNA_assay, input_read_RNA_assay)
    colnames(input_read_RNA_assay)[2]= "dummy___"
  }
  
  blueprint <- celldex::BlueprintEncodeData(
    ensembl = feature_nomenclature == "ensembl"
    #legacy = TRUE
  )
  
  data_annotated =
    
    input_read_RNA_assay |>
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
      
      input_read_RNA_assay |>
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
  
  MonacoImmuneData <- celldex::MonacoImmuneData(
    ensembl = feature_nomenclature == "ensembl"
    #legacy = TRUE
  )
  
  data_annotated =
    data_annotated |>
    
    left_join(
      input_read_RNA_assay |>
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
      input_read_RNA_assay |>
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
    
    library(Seurat) # !!! If this is not here gives error, but this has to go for Bioconductor
    
    # Convert SCE to SE to calculate SCT
    if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
      
      assay(input_read_RNA_assay, assay) <- 
        assay(input_read_RNA_assay, assay) |> 
        as("dgCMatrix")
      
      input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL,  counts = assay) 
      
      # Rename assay
      assay_name_old = input_read_RNA_assay |> Assays() |> _[[1]]
      input_read_RNA_assay = input_read_RNA_assay |>
        RenameAssays(
          assay.name = assay_name_old,
          new.assay.name = assay)
    } 
    
    options(future.globals.maxSize = 16 * 1024^3)
    
    azimuth_annotation = 
      tryCatch({
        
        if(ncol(input_read_RNA_assay)<200) k.weight = 25
        else k.weight = 50
        
        input_read_RNA_assay |> RenameAssays(assay.name = assay, new.assay.name = "RNA") |> 
          Azimuth::RunAzimuth(reference = reference_azimuth, assay = "RNA", umap.name = "refUMAP") |> 
          as_tibble() |>
          dplyr::select(.cell, any_of(
            c(
              "predicted.celltype.l1",
              "predicted.celltype.l2",
              "predicted.celltype.l3",
              "predicted.celltype.l1.score",
              "predicted.celltype.l2.score",
              "predicted.celltype.l3.score"
            )
          )) |> 
          nest(azimuth_scores_celltype = c(ends_with("score"))) |>
          dplyr::rename(azimuth_predicted.celltype.l1 = predicted.celltype.l1,
                        azimuth_predicted.celltype.l2 = predicted.celltype.l2,
                        azimuth_predicted.celltype.l3 = predicted.celltype.l3)
      },
      error = function(e) {
        if(!str_detect(e$message, "Please set k.weight to be smaller than the number of anchors|Number of anchor cells is less than k.weight|number of items to replace is not a multiple of replacement length")) 
          stop("HPCell says: Seurat Azimuth failed, probably for the small number of cells, which is .", ncol(input_read_RNA_assay), " Please investigate -> ", e$message)
        print(e)
        input_read_RNA_assay |> as_tibble() |> dplyr::select(.cell)
      })
    
    # Save
    data_annotated  |>
      left_join(azimuth_annotation, by = dplyr::join_by(.cell)	)
    

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
#' @importFrom magrittr not
#' @importFrom Matrix colSums
#' @importFrom magrittr extract2 not
#' @importFrom SummarizedExperiment assay colData assay<-
#' @importFrom tidyselect all_of
#' 
#' @export
alive_identification <- function(input_read_RNA_assay,
                                 empty_droplets_tbl = NULL,
                                 annotation_label_transfer_tbl = NULL,
                                 annotation_column = NULL,
                                 assay = NULL,
                                 feature_nomenclature) {
  
  # Fix GCHECK notes
  empty_droplet = NULL
  detected = NULL
  .cell = NULL 
  high_mitochondrion = NULL 
  
  if(ncol(input_read_RNA_assay) == 0) {
    warning("HPCell says: the sample ", sample_name, " has no cells and returned NULL")
    return(NULL)
  }
  
  if(input_read_RNA_assay |> is.null()) {
    warning("HPCell says: the sample ", sample_name, " is NULL and returned NULL")
    return(NULL)
  }
  
  if(
    !is.null(annotation_column) && 
    !annotation_column %in% colnames(as_tibble(input_read_RNA_assay[1,1]))
  )
    stop("HPCell says: Your `group_by` columns are not present in your data. Please run annotate_cell_type_hpc() to get the cell type annotation that you can use as grouping for the cell-type-specific quality control and removal of dead cells.")
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  if (!is.null(empty_droplets_tbl)) {
    input_read_RNA_assay =
      input_read_RNA_assay |>
      left_join(empty_droplets_tbl, by=".cell") |>
      dplyr::filter(!empty_droplet)
  } 
  
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
    counts <- SummarizedExperiment::assay(input_read_RNA_assay, assay = assay)
    if (!any(str_which(colnames(colData(input_read_RNA_assay)), nFeature_name)) ||
        !any(str_which(colnames(colData(input_read_RNA_assay)), nCount_name))) {
      colData(input_read_RNA_assay)[[nFeature_name]] <- Matrix::colSums(counts > 0)
      colData(input_read_RNA_assay)[[nCount_name]] <- Matrix::colSums(counts)
    }
  }
  
  
  
  # Returns a named vector of IDs
  # Matches the gene id’s row by row and inserts NA when it can’t find gene names
  if (feature_nomenclature == "symbol") {
    location <- mapIds(
      EnsDb.Hsapiens.v86,
      keys=rownames(input_read_RNA_assay),
      column="SEQNAME",
      keytype="SYMBOL"
    )
  }
  
  
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
    rna_counts <- SummarizedExperiment::assay(input_read_RNA_assay, assay=assay)
    SummarizedExperiment::assay(input_read_RNA_assay, assay) <- 
      SummarizedExperiment::assay(input_read_RNA_assay, assay) |> as("dgCMatrix")
    
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL, 
                                                              counts = assay) 
    
    # Rename assay
    assay_name_old = input_read_RNA_assay |> Assays() |> _[[1]]
    input_read_RNA_assay = input_read_RNA_assay |>
      RenameAssays(
        assay.name = assay_name_old,
        new.assay.name = assay)
  }
  
  # Compute per-cell QC metrics
  qc_metrics <- perCellQCMetrics(rna_counts, subsets=list(Mito=which_mito)) %>%
    as_tibble(rownames = ".cell") %>%
    dplyr::select(-sum, -detected)
  
  # I HAVE TO DROP UNIQUE, AS SOON AS THE BUG IN SEURAT IS RESOLVED. UNIQUE IS BUG PRONE HERE.
  percentage_output = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay)
  percentage_output = percentage_output[!duplicated(names(percentage_output))]
  # Compute ribosome statistics
  ribosome =
    input_read_RNA_assay |>
    select(.cell) |>
    #mutate(subsets_Ribo_percent = PercentageFeatureSet(input_read_RNA_assay,  pattern = "^RPS|^RPL", assay = assay)[,1]) |>
    mutate(subsets_Ribo_percent = percentage_output)
  
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
                                   empty_droplets_tbl = NULL, 
                                   alive_identification_tbl = NULL, 
                                   #annotation_label_transfer_tbl, 
                                   #reference_label_fine,
                                   assay = NULL){
  
  # Fix GChecks 
  .cell = NULL 
  empty_droplet = NULL 
  
  if(ncol(input_read_RNA_assay) == 0) {
    warning("HPCell says: the sample ", sample_name, " has no cells and returned NULL")
    return(NULL)
  }
  
  if(input_read_RNA_assay |> is.null()) {
    warning("HPCell says: the sample ", sample_name, " is NULL and returned NULL")
    return(NULL)
  }
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  
  if (inherits(input_read_RNA_assay, "Seurat")) {
    input_read_RNA_assay <- input_read_RNA_assay |>
      # Filtering empty
      Seurat::as.SingleCellExperiment() 
  } 
  
  if (!is.null(empty_droplets_tbl)) { 
    
    filter_empty_droplets <- input_read_RNA_assay |>
      # Filtering empty
      left_join(empty_droplets_tbl |> select(.cell, empty_droplet), by = ".cell") |>
      filter(!empty_droplet) 
    
    # Filtering dead
    if(alive_identification_tbl |> is.null() |> not())
      input_read_RNA_assay = input_read_RNA_assay |>
        left_join(alive_identification_tbl |> select(.cell, alive), by = ".cell") |>
        filter(alive)
  } 
  
  # Condition as scDblFinder only accept assay "counts"
  if (!"counts" %in% (SummarizedExperiment::assays(filter_empty_droplets) |> names())){
    SummarizedExperiment::assay(filter_empty_droplets, "counts") <-  
      SummarizedExperiment::assay(filter_empty_droplets, assay)
    SummarizedExperiment::assay(filter_empty_droplets, assay) <- NULL
  }
  # Annotate
  input_read_RNA_assay |> 
    #left_join(annotation_label_transfer_tbl, by = ".cell")|>
    #scDblFinder(clusters = ifelse(reference_label_fine=="none", TRUE, reference_label_fine)) |>
    scDblFinder(clusters = NULL) |> 
    colData() |> 
    as_tibble(rownames = ".cell") |> 
    select(.cell, contains("scDblFinder")) 
  
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
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom Seurat CellCycleScoring
#' @importFrom Seurat as.Seurat
#' @importFrom SeuratObject RenameAssays
#' @importFrom Seurat NormalizeData
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
cell_cycle_scoring <- function(input_read_RNA_assay, 
                               empty_droplets_tbl = NULL,
                               feature_nomenclature,
                               assay = NULL){
  #Fix GCHECK
  empty_droplet = NULL 
  .cell = NULL 
  S.Score = NULL 
  G2M.Score = NULL 
  Phase = NULL 
  
  if(ncol(input_read_RNA_assay) == 0) {
    warning("HPCell says: the sample ", sample_name, " has no cells and returned NULL")
    return(NULL)
  }
  
  if(input_read_RNA_assay |> is.null()) {
    warning("HPCell says: the sample ", sample_name, " is NULL and returned NULL")
    return(NULL)
  }
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  # Convert to Seurat in order to perform cell cycle scoring
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    assay(input_read_RNA_assay, assay) <- assay(input_read_RNA_assay, assay) |> 
      as("dgCMatrix")
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL, 
                                                              counts = assay) 
    
    # Rename assay
    assay_name_old = input_read_RNA_assay |> Assays() |> _[[1]]
    input_read_RNA_assay = input_read_RNA_assay |>
      RenameAssays(
        assay.name = assay_name_old,
        new.assay.name = assay)
  }
  
  if (feature_nomenclature == "ensembl") {
    s.features_tidy = Seurat::cc.genes$s.genes |> 
      convert_gene_names(current_nomenclature = "symbol") |>
      filter(stringr::str_detect(gene_id, "ENSG*")) |> dplyr::pull(gene_id)
    g2m.features_tidy = Seurat::cc.genes$g2m.genes |> 
      convert_gene_names(current_nomenclature = "symbol") |>
      filter(stringr::str_detect(gene_id, "ENSG*")) |> dplyr::pull(gene_id)
  } else if (feature_nomenclature == "symbol") {
    s.features_tidy = Seurat::cc.genes$s.genes
    g2m.features_tidy = Seurat::cc.genes$g2m.genes
  }
  
  # avoid small number of cells 
  if (!is.null(empty_droplets_tbl)) {
    filtered_counts <- input_read_RNA_assay |>
      left_join(empty_droplets_tbl, by = ".cell") |>
      dplyr::filter(!empty_droplet)
  } 
  
  counts <- filtered_counts |>
    # Normalise needed
    NormalizeData() |>
    
    # Assign cell cycle scores of each cell
    # Based on its expression of G2/M and S phase markers
    #Stores S and G2/M scores in object meta data along with predicted classification of each cell in either G2M, S or G1 phase
    CellCycleScoring(s.features = s.features_tidy,
                     g2m.features = g2m.features_tidy,
                     set.ident = FALSE) |>
    
    as_tibble() |>
    select(.cell, S.Score, G2M.Score, Phase) 
  
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
#' @importFrom SummarizedExperiment assay assay<-
#' @export
non_batch_variation_removal <- function(input_read_RNA_assay, 
                                        empty_droplets_tbl = NULL, 
                                        alive_identification_tbl = NULL, 
                                        cell_cycle_score_tbl = NULL,
                                        assay = NULL,
                                        factors_to_regress = NULL,
                                        external_path){
  #Fix GChecks 
  empty_droplet = NULL 
  .cell <- NULL 
  
  if(ncol(input_read_RNA_assay) == 0) {
    warning("HPCell says: the sample ", sample_name, " has no cells and returned NULL")
    return(NULL)
  }
  
  if(input_read_RNA_assay |> is.null()) {
    warning("HPCell says: the sample ", sample_name, " is NULL and returned NULL")
    return(NULL)
  }
  
  # Your code for non_batch_variation_removal function here
  class_input = input_read_RNA_assay |> class()
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    assay(input_read_RNA_assay, assay) <- assay(input_read_RNA_assay, assay) |> as("dgCMatrix")
    
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL, 
                                                              counts = assay) 
    
    # Rename assay
    assay_name_old = input_read_RNA_assay |> Assays() |> _[[1]]
    input_read_RNA_assay_transform = input_read_RNA_assay |>
      RenameAssays(
        assay.name = assay_name_old,
        new.assay.name = assay)
  }
  
  # avoid small number of cells 
  if (!is.null(empty_droplets_tbl)) {
    filtered_counts <- input_read_RNA_assay_transform |>
      left_join(empty_droplets_tbl, by = ".cell") |>
      dplyr::filter(!empty_droplet)
  } 
  
  counts =
    filtered_counts |>
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
  normalized_rna <- 
    input_read_RNA_assay |> 
    Seurat::SCTransform(
      counts, 
      assay=assay,
      return.only.var.genes=FALSE,
      residual.features = NULL,
      vars.to.regress = factors_to_regress,
      vst.flavor = "v2",
      scale_factor=2186,  
      conserve.memory=T, 
      min_cells=0,
    )  |> 
    GetAssayData(assay="SCT")
  
  
  if (class_input == "SingleCellExperiment") {
    
    write_HDF5_array_safe(normalized_rna, "SCT", external_path)
    
  } else if (class_input ==  "Seurat") {
    
    normalized_rna 
    
  }
  
  
  # # Normalise antibodies
  # if ( "ADT" %in% names(normalized_rna@assays)) {
  #   normalized_data <- normalized_rna %>%
  #     NormalizeData(normalization.method = 'CLR', margin = 2, assay="ADT") %>%
  #     select(-subsets_Ribo_percent, -subsets_Mito_percent, -G2M.Score)
  #   
  #   my_assays = my_assays |> c("CLR")
  #   
  # } else { 
  #   normalized_data <- normalized_rna %>%
  #     # Drop alive columns
  #     select(-subsets_Ribo_percent, -subsets_Mito_percent, -G2M.Score)
  # }
  
  
  
  
  
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
#' @importFrom dplyr left_join filter select
#' @import SeuratObject
#' @importFrom SummarizedExperiment left_join assay assay<-
#' @import tidySingleCellExperiment 
#' @import tidyseurat
#' @importFrom magrittr not
#' @importFrom SingleCellExperiment altExp
#' @importFrom SingleCellExperiment altExp<-
#' @export
preprocessing_output <- function(input_read_RNA_assay,
                                 empty_droplets_tbl = NULL,
                                 non_batch_variation_removal_S = NULL, 
                                 alive_identification_tbl = NULL, 
                                 cell_cycle_score_tbl = NULL, 
                                 annotation_label_transfer_tbl = NULL, 
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
  
  if(ncol(input_read_RNA_assay) == 0) {
    warning("HPCell says: the sample ", sample_name, " has no cells and returned NULL")
    return(NULL)
  }
  
  if(input_read_RNA_assay |> is.null()) {
    warning("HPCell says: the sample ", sample_name, " is NULL and returned NULL")
    return(NULL)
  }
  
  if (empty_droplets_tbl |> is.null() |> not()) {
    input_read_RNA_assay =
      input_read_RNA_assay |>
      left_join(empty_droplets_tbl, by = ".cell") |>
      filter(!empty_droplet)
  } 
  
  # Add normalisation
  if(!is.null(non_batch_variation_removal_S)){
    if(input_read_RNA_assay |> is("Seurat"))
      input_read_RNA_assay[["SCT"]] = non_batch_variation_removal_S
    else if(input_read_RNA_assay |> is("SingleCellExperiment")){
      message("HPCell says: in order to attach SCT assay to the SingleCellExperiment, SCT was added to external experiments slot")
      #input_read_RNA_assay = input_read_RNA_assay[rownames(non_batch_variation_removal_S), 
      # altExp(input_read_RNA_assay) = SingleCellExperiment(assay  = list(SCT = non_batch_variation_removal_S))
      assay(input_read_RNA_assay, "SCT") <- non_batch_variation_removal_S
      
    }
  }
  
  
  # Filtering dead
  if(alive_identification_tbl |> is.null() |> not())
    input_read_RNA_assay = input_read_RNA_assay |>
    left_join(alive_identification_tbl |> select(.cell, alive), by = ".cell") |>
    filter(alive) 
  
  
  
  # Filter doublets
  if(doublet_identification_tbl |> is.null() |> not())
    input_read_RNA_assay <- input_read_RNA_assay |>
    left_join(doublet_identification_tbl |> select(.cell, scDblFinder.class), by = ".cell") |>
    filter(scDblFinder.class=="singlet") 
  
  # attach cell cycle
  if(cell_cycle_score_tbl |> is.null() |> not())
    input_read_RNA_assay = 
    input_read_RNA_assay |>
    left_join(
      cell_cycle_score_tbl ,
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
#' @importFrom digest digest
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' 
#' @export

# Create pseudobulk for each sample 
create_pseudobulk <- function(input_read_RNA_assay, 
                              sample_names_vec, 
                              empty_droplets_tbl = NULL,
                              alive_identification_tbl = NULL,
                              cell_cycle_score_tbl = NULL,
                              annotation_label_transfer_tbl = NULL,
                              doublet_identification_tbl = NULL,  
                              x = c() , 
                              external_path, assays = NULL) {
  #Fix GChecks 
  .sample = NULL 
  .feature = NULL 
  data_source = NULL 
  symbol = NULL 
  
  dir.create(external_path, showWarnings = FALSE, recursive = TRUE)
  
  preprocessing_output_S = 
    preprocessing_output(
      input_read_RNA_assay,
      empty_droplets_tbl,
      non_batch_variation_removal_S = NULL, 
      alive_identification_tbl, 
      cell_cycle_score_tbl, 
      annotation_label_transfer_tbl, 
      doublet_identification_tbl
    )
  
  
  if(assays |> is.null()){
    if(preprocessing_output_S |> is("Seurat"))
      assays = Seurat::Assays(preprocessing_output_S)
    else if(preprocessing_output_S |> is("SingleCellExperiment"))
      assays = preprocessing_output_S@assays |> names()
    
  }
  
  # Aggregate cells
  pseudobulk = 
    preprocessing_output_S |> 
    
    # Add sample
    mutate(sample_hpc = sample_names_vec) |> 
    
    # Aggregate
    tidySingleCellExperiment::aggregate_cells(c(sample_hpc, !!sym(x)), 
                                              slot = "data", assays = assays) 
  
  # If I start from Seurat
  if(pseudobulk |> is("data.frame"))
    pseudobulk = pseudobulk |>
    as_SummarizedExperiment(.sample, .feature, any_of(assays)) 
  
  rowData(pseudobulk)$feature_name = rownames(pseudobulk)
  
  pseudobulk = pseudobulk |>
    pivot_longer(cols = assays, names_to = "data_source", values_to = "count") |>
    filter(!count |> is.na()) |>
    
    # Some manipulation to get unique feature because RNA and ADT
    # both can have same name genes
    rename(symbol = .feature) |>
    mutate(data_source = stringr::str_remove(data_source, "abundance_")) |>
    unite(".feature", c(symbol, data_source), remove = FALSE) |>
    
    # Covert
    as_SummarizedExperiment(
      .sample = .sample,
      .transcript = .feature,
      .abundance = count
    ) 
  
  file_name = glue("{external_path}/{digest(pseudobulk)}")
  
  pseudobulk |>
    
    # Conver to H5
    saveHDF5SummarizedExperiment(dir = file_name, replace=TRUE, as.sparse=TRUE)
  
}

#' Merge pseudobulk from all samples 
#'
#' @description
#' Merge pseudobulk from all samples. Ensures that missing genes are accounted 
#' for and aligns data across multiple samples.
#' @param pseudobulk_list A list pseudobulk samples generated by `create_pseudobulk`
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
pseudobulk_merge <- function(pseudobulk_list, external_path, ...) {
  
  dir.create(external_path, showWarnings = FALSE, recursive = TRUE)
  
  
  # Fix GCHECKS 
  . = NULL 
  
  # Select only common columns
  # investiagte common_columns, as data_source is not a common column in the pilot data
  common_columns =
    pseudobulk_list |>
    purrr::map(~ .x |> as_tibble() |> colnames()) |>
    unlist() |>
    table() %>%
    .[.==max(.)] |>
    names()
  
  # All genes 
  all_genes =
    pseudobulk_list |>
    purrr::map(~ .x |> rownames()) |>
    unlist() |>
    unique() |>
    as.character()
  
  
  se <- pseudobulk_list |>
    
    # Add missing genes
    purrr::map(~{
      missing_genes = all_genes |> setdiff(rownames(.x))
      
      if(missing_genes |> length() == 0) return(.x)
      else
        .x |> add_missingh_genes_to_se(all_genes, missing_genes)
      
    }) |>
    
    purrr::map(~ .x |> dplyr::select(any_of(common_columns)))   %>%
    
    do.call(S4Vectors::cbind, .) 
  
  
  file_name = glue("{external_path}/{digest(se)}")
  
  se = 
    se |> 
    
    saveHDF5SummarizedExperiment(dir = file_name, replace=TRUE, as.sparse=TRUE) 
  
  # Return the pseudobulk data for this single sample
  return(se)
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
#' @importFrom tidybulk test_differential_abundance
#' 
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


#' @export
is_target = function(x) {
  
  if(x |> is.null()) return(NULL)
  
  if(x |> is("character") |> not())
    stop("HPCell says: the input to `is_target` must be a character")
  
  as.name(x) 
} 
