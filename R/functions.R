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
                             feature_nomenclature){
  
  if(input_read_RNA_assay |> is.null()) return(NULL)
  if(ncol(input_read_RNA_assay) == 0) return(NULL)
  
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
  #   ggplot2::ggplot(aes(rank, total)) 
  #   geom_point(aes(color = empty_droplet, size = empty_droplet )) 
  #   geom_line(aes(rank, fitted), color="purple") 
  #   geom_hline(aes(yintercept = knee), color="dodgerblue") 
  #   geom_hline(aes(yintercept = inflection), color="forestgreen") 
  #   scale_x_log10() 
  #   scale_y_log10() 
  #   scale_color_manual(values = c("black", "#e11f28")) 
  #   scale_size_discrete(range = c(0, 2)) 
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
#' `empty_droplet_threshold` identifies empty droplets by applying a gene expression threshold per sample.
#'   It excludes mitochondrial and ribosomal genes, and classifies droplets as empty if 
#'   the number of expressed genes falls below the specified threshold.
#'   
#'   The function returns a tibble containing the number of expressed genes, 
#'   total RNA count for each cell, and a logical annotation indicating whether the droplet was classified as empty.
#'
#' @param input_read_RNA_assay SingleCellExperiment or Seurat object containing RNA assay data.
#' @param filter_empty_droplets Logical value indicating whether to filter the input data.
#' @param RNA_feature_threshold An optional integer for the number of feature expressed in a sample.
#'
#' @return A tibble with columns: Cell, nFeature_expressed_in_sample, nCount_RNA, empty_droplet (classification of droplets).
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
                                   RNA_feature_threshold){
  if(input_read_RNA_assay |> is.null()) return(NULL)
  if(ncol(input_read_RNA_assay) == 0) return(NULL)
  
  #Fix GChecks 
  FDR = NULL 
  .cell = NULL 
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  significance_threshold = 0.001
  
  # # Rule of thumb threshold
  # expressed_genes_threshold = 0.025
  # if (is.null(RNA_feature_threshold)) RNA_feature_threshold = min(floor(dim(input_read_RNA_assay)[1]*expressed_genes_threshold), 500)
  # 
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
  
  # Generate library size for each cell
  library_size <- colSums(filtered_counts) |> enframe(name = ".cell", value = "nCount_RNA")
  
  # filter based on number of expressed genes
  result <- colSums(filtered_counts > 0 ) |> enframe(name = ".cell", value = "nFeature_expressed_in_sample") |> 
    left_join(library_size, by = ".cell") |>
    mutate(empty_droplet = nFeature_expressed_in_sample < RNA_feature_threshold)
  
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
#' @importFrom Seurat CreateAssayObject
#' @importFrom Seurat SCTransform
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat FindTransferAnchors
#' @importFrom Seurat MapQuery
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom Seurat Assays
#' @importFrom SeuratObject RenameAssays
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
#' 
#' @export
annotation_label_transfer <- function(input_read_RNA_assay,
                                      empty_droplets_tbl = NULL, 
                                      reference_azimuth = NULL,
                                      assay = NULL,
                                      feature_nomenclature
){
  # Fix github checks 
  empty_droplet = NULL 
  pruned.labels = NULL 
  delta.next = NULL 
  .cell = NULL 
  
  if(input_read_RNA_assay |> is.null()) return(NULL)
  if(ncol(input_read_RNA_assay) == 0) return(NULL)
  
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
    input_read_RNA_assay =
      input_read_RNA_assay |>
      as.SingleCellExperiment() |>
      logNormCounts()
  } else if (inherits(input_read_RNA_assay, "SingleCellExperiment")){
    input_read_RNA_assay =
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
      
      assay = input_read_RNA_assay@assays |> names() |> extract2(1)
      
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
#' @param cell_type_ensembl_harmonised_tbl A tibble with annotated cell type label data.
#' @param cell_type_column A character vector indicating the cell type column used for grouping during quality control and dead cell removal.
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
#' @importFrom magrittr extract2
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment colData
#' @importFrom tidyselect all_of
#' 
#' @export
alive_identification <- function(input_read_RNA_assay,
                                 empty_droplets_tbl = NULL,
                                 cell_type_ensembl_harmonised_tbl = NULL,
                                 cell_type_column = NULL,
                                 assay = NULL,
                                 feature_nomenclature) {
  
  # Fix GCHECK notes
  empty_droplet = NULL
  detected = NULL
  .cell = NULL 
  high_mitochondrion = NULL 
  
  if (input_read_RNA_assay |> is.null()) return(NULL)
  
  if(
    is.null(cell_type_column) && 
    !cell_type_column %in% colnames(as_tibble(input_read_RNA_assay[1,1]))
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
  
  # In rare cases, all cells in a sample are empty droplets
  if (ncol(input_read_RNA_assay) == 0) return(NULL)
  
  # In rare cases, a cell in a sample is non empty droplet
  if (ncol(input_read_RNA_assay) == 1) input_read_RNA_assay = input_read_RNA_assay |> duplicate_single_column_assay()
  
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
  
  input_read_RNA_assay <- input_read_RNA_assay %>%
    AddMetaData(
      metadata = PercentageFeatureSet(input_read_RNA_assay, pattern = "^MT-", assay = assay), 
      col.name = "percent.mt"
    )
  
  # Returns a named vector of IDs
  # Matches the gene id's row by row and inserts NA when it can't find gene names
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
  #   # Join cell types if cell_type_ensembl_harmonised_tbl provided
  #   {\(x)
  #     if (inherits(cell_type_ensembl_harmonised_tbl, "tbl_df")) {
  #       left_join(x, cell_type_ensembl_harmonised_tbl, by = ".cell") |>
  # 
  #       # Label cells
  #       nest(data = -all_of(cell_type_column)) |>
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
    
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL) 
    
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
  
  # Add cell type labels and determine high mitochondrion content, if cell_type_ensembl_harmonised_tbl is provided
  if(cell_type_column |> is.null() |> not()) {
    
    if (
      inherits(cell_type_ensembl_harmonised_tbl, "tbl_df") &&
      cell_type_column %in% colnames(cell_type_ensembl_harmonised_tbl)
    ) {
      
      mitochondrion <- qc_metrics %>%
        left_join(cell_type_ensembl_harmonised_tbl, by = ".cell") 
      
      ribosome =
        ribosome |>
        # Only retrieve metadata so nesting in the next step won't break
        left_join(cell_type_ensembl_harmonised_tbl, by = ".cell") |> as_tibble()
    }
    
    
    else if (cell_type_column %in% colnames(as_tibble(input_read_RNA_assay[1,1]))) {
      
      mitochondrion <- 
        qc_metrics %>%
        left_join(input_read_RNA_assay |> as_tibble() |> select(.cell, all_of(cell_type_column)), by = ".cell") 
      
      
      ribosome = 
        ribosome |>
        # Only retrieve metadata so nesting in the next step won't break
        left_join(input_read_RNA_assay |> as_tibble() |> select(.cell, all_of(cell_type_column)), by = ".cell") |> as_tibble()
        
    }
    
    
    mitochondrion = 
      mitochondrion %>%
      nest(data = -all_of(cell_type_column)) 
    
    ribosome =
      ribosome |>
      nest(data = -all_of(cell_type_column)) 
    
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
    left_join(ribosome) |>
    mutate(alive = !high_mitochondrion) |> # & !high_ribosome ) |>
  # Select informative columns
    select(.cell, {{cell_type_column}}, contains("subsets"), contains("observation"),
           contains("high"), alive)
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
#' @param assay Name of the assay to use.
#'
#' @return A tibble containing cells with their scDblFinder scores.
#'
#' @importFrom dplyr left_join filter select
#' @importFrom Matrix Matrix 
#' @importFrom SummarizedExperiment colData assayNames assayNames<-
#' @importFrom Seurat as.SingleCellExperiment
#' @import scDblFinder
#' @export
doublet_identification <- function(input_read_RNA_assay, 
                                   empty_droplets_tbl = NULL, 
                                   # annotation_label_transfer_tbl,
                                   # reference_label_fine,
                                   assay = NULL){
  
  # Fix GChecks 
  .cell = NULL 
  empty_droplet = NULL 
  
  if (is.null(input_read_RNA_assay)) return(NULL)
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  
  if (inherits(input_read_RNA_assay, "Seurat")) {
    input_read_RNA_assay <- input_read_RNA_assay |>
      Seurat::as.SingleCellExperiment() 
    
  } 
  # Filtering empty
  if (!is.null(empty_droplets_tbl)) { 
    input_read_RNA_assay <- input_read_RNA_assay |>

      left_join(empty_droplets_tbl |> select(.cell, empty_droplet), by = ".cell") |>
      filter(!empty_droplet) 
  }
  
  # scDblFinder() can identify doublets from non-empty droplet cells, so no need to filter alive
  
  # In rare cases, all cells in a sample are empty droplets or dead
  if (ncol(input_read_RNA_assay) == 0) return(NULL)
  
  # scDblFinder can only handle counts assay, thus rename
  assayNames(input_read_RNA_assay)[assayNames(input_read_RNA_assay) == assay] <- "counts"
  
  # Annotate
  # Mark doublets to Unknown when scDblFinder fails and keep them in the downstream analysis.
  result <- tryCatch({
    # Run scDblFinder
    input_read_RNA_assay <- input_read_RNA_assay |>
      # By default, artificial doublets will be considered unidentifiable when score threshold below 0.2
      scDblFinder(clusters = NULL)
  }, error = function(e) {
    # Error handling
    message("Error in scDblFinder: ", e$message)
    input_read_RNA_assay <- input_read_RNA_assay |> mutate(scDblFinder.class = "Unknown")
    }
  )

  result |>
    colData() |>
    as_tibble(rownames = ".cell") |>
    select(.cell, scDblFinder.class)

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
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> extract2(1)
  
  # Convert to Seurat in order to perform cell cycle scoring
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    assay(input_read_RNA_assay, assay) <- assay(input_read_RNA_assay, assay) |> as("dgCMatrix")
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL) 
    
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
    input_read_RNA_assay = input_read_RNA_assay |>
      RenameAssays(
        assay.name = assay_name_old,
        new.assay.name = assay)
  }
  
  # avoid small number of cells 
  if (!is.null(empty_droplets_tbl)) {
    input_read_RNA_assay <- input_read_RNA_assay |>
      left_join(empty_droplets_tbl, by = ".cell") |>
      dplyr::filter(!empty_droplet)
  } 

  if (!is.null(alive_identification_tbl)) {
  input_read_RNA_assay =
    input_read_RNA_assay |>
    left_join(
      alive_identification_tbl |>
        select(.cell, any_of(factors_to_regress)),
      by=".cell"
    ) 
  }
  
  if(!is.null(cell_cycle_score_tbl)) 
    input_read_RNA_assay = input_read_RNA_assay |>
    
    left_join(
      cell_cycle_score_tbl |>
        select(.cell, any_of(factors_to_regress)),
      by=".cell"
    )
  
  # filter(!high_mitochondrion | !high_ribosome)
  
  # variable_features = readRDS(input_path_merged_variable_genes)
  # 
  # # Set variable features
  # VariableFeatures(input_read_RNA_assay) = variable_features
  
  # Normalise RNA
  input_read_RNA_assay <- 
    input_read_RNA_assay |> 
      Seurat::SCTransform(
      assay=assay,
      return.only.var.genes=FALSE,
      residual.features = NULL,
      vars.to.regress = factors_to_regress,
      vst.flavor = "v2",
      scale_factor=2186,  
      conserve.memory=T, 
      min_cells=0
    )  |> 
    GetAssayData(assay="SCT")

  if (class_input == "SingleCellExperiment") {

    if(input_read_RNA_assay[,1,drop=FALSE] |> is.nan() |> any())
             warning("HPCell says: some features might be all 0s, NaN are added by Seurat in the SCT assay, and kept in the assay because SingleCellExperiment requires same feature set for all assays.")
             
    write_HDF5_array_safe(input_read_RNA_assay, "SCT", external_path)
    
  } else if (class_input ==  "Seurat") {

    # Remove NaN features from SCT assay
    input_read_RNA_assay <- input_read_RNA_assay[!apply(input_read_RNA_assay, 1, function(row) all(is.nan(row))), ]
                                                      
    input_read_RNA_assay 
    
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

#' Preprocess metacells with the SuperCell approach
#' 
#' This function preprocesses a single-cell gene expression matrix for downstream simplification using PCA 
#' and k-nearest neighbor (kNN) graph construction. It includes options for scaling, feature selection, 
#' approximate sampling, and PCA computation methods.
#'
#' @param input_read_RNA_assay A `SingleCellExperiment` or `Seurat` object containing RNA assay data.
#' @param assay assay used, default = "RNA" 
#' @param genes.use a vector of genes used to compute PCA
#' @param genes.exclude a vector of genes to be excluded when computing PCA
#' @param n.var.genes if \code{"genes.use"} is not provided, \code{"n.var.genes"} genes with the largest variation are used
#' @param k.knn parameter to compute single-cell kNN network
#' @param do.scale whether to scale gene expression matrix when computing PCA
#' @param n.pc number of principal components to use for construction of single-cell kNN network
#' @param fast.pca use \link[irlba]{irlba} as a faster version of prcomp (one used in Seurat package)
#' @param do.approx compute approximate kNN in case of a large dataset (>50'000)
#' @param approx.N number of cells to subsample for an approximate approach. By default, 5000 cells are used 
#'   for approximation to capture biological meaningful result.
#' @param seed seed to use to subsample cells for an approximate approach
#' @param ... other parameters of \link{build_knn_graph} function
#' @return A list of variables to be passed to the `SuperCell::SCimplify` gamma involved function. 
#' @importFrom Matrix t
#' @importFrom stats var prcomp
#' @importFrom irlba irlba
#' @importFrom SuperCell build_knn_graph
#' @export
preprocess_SCimplify <- function(input_read_RNA_assay,
                                 assay = NULL,
                                 genes.use = NULL,
                                 genes.exclude = NULL,
                                 n.var.genes = min(1000, nrow(input_read_RNA_assay)),
                                 k.knn = 5,
                                 do.scale = TRUE,
                                 n.pc = 10,
                                 fast.pca = TRUE,
                                 do.approx = FALSE,
                                 approx.N = 5000,
                                 seed = 12345,
                                 ...){
  #Fix GChecks 
  empty_droplet = NULL 
  .cell <- NULL 
  
  # Your code for non_batch_variation_removal function here
  class_input = input_read_RNA_assay |> class()
  
  # For small number of cells
  if (ncol(input_read_RNA_assay) < 10) {
    k.knn = ncol(input_read_RNA_assay) - 1
    n.pc = ncol(input_read_RNA_assay) - 1
  }
  
  # Get assay
  if(is.null(assay)) assay = input_read_RNA_assay@assays |> names() |> magrittr::extract2(1)
  
  # Convert to SE if the input is SCE
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    assay(input_read_RNA_assay, assay) <- assay(input_read_RNA_assay, assay) |> as("dgCMatrix")
    
    input_read_RNA_assay <- input_read_RNA_assay |> as.Seurat(data = NULL, 
                                                              counts = assay) 
    
    # Rename assay
    assay_name_old = DefaultAssay(input_read_RNA_assay)
    input_read_RNA_assay_transform = input_read_RNA_assay |>
      RenameAssays(
        assay.name = assay_name_old,
        new.assay.name = assay)
  }
  
  # Get normalise and scale gene expression matrix with rows to be genes and cols to be cells
  normalized_rna <- 
    input_read_RNA_assay |> 
    NormalizeData(normalization.method = "LogNormalize") |> 
    FindVariableFeatures(nfeatures = 2000) |>
    ScaleData() |>
    RunPCA(npcs = min(50, ncol(input_read_RNA_assay) - 1), verbose = F) |> 
    RunUMAP(reduction = "pca", dims = c(1:min(30, ncol(input_read_RNA_assay) - 1)), 
            n.neighbors = min(30, ncol(input_read_RNA_assay) - 1), verbose = F) |> 
    Seurat::GetAssayData(slot = "data")
  
  N.c <- ncol(normalized_rna)
  
  # if(gamma > 100 & N.c < 100000){
  #   warning(paste0("Graining level (gamma = ", gamma, ") seems to be very large! Please, consider using smaller gamma, the suggested range is 10-50."))
  # }
  
  if(is.null(rownames(normalized_rna))){
    if(!(is.null(genes.use) | is.null(genes.exclude))){
      stop("rownames(normalized_rna) is Null \nGene expression matrix normalized_rna is expected to have genes as rownames")
    } else {
      warning("colnames(normalized_rna) is Null, \nGene expression matrix normalized_rna is expected to have genes as rownames! \ngenes will be created automatically in a form 'gene_i' ")
      rownames(normalized_rna) <- paste("gene", 1:nrow(normalized_rna), sep = "_")
    }
  }
  
  if(is.null(colnames(normalized_rna))){
    warning("colnames(normalized_rna) is Null, \nGene expression matrix normalized_rna is expected to have cellIDs as colnames! \nCellIDs will be created automatically in a form 'cell_i' ")
    colnames(normalized_rna) <- paste("cell", 1:N.c, sep = "_")
  }
  
  cell.ids <- colnames(normalized_rna)
  
  keep.genes    <- setdiff(rownames(normalized_rna), genes.exclude)
  normalized_rna             <- normalized_rna[keep.genes,]
  
  
  if(is.null(genes.use)){
    n.var.genes <- min(n.var.genes, nrow(normalized_rna))
    if(N.c > 50000){
      set.seed(seed)
      idx         <- sample(N.c, 50000)
      gene.var    <- apply(normalized_rna[,idx], 1, stats::var)
    } else {
      gene.var    <- apply(normalized_rna, 1, stats::var)
    }
    
    genes.use   <- names(sort(gene.var, decreasing = TRUE))[1:n.var.genes]
  }
  
  if(length(intersect(genes.use, genes.exclude)) > 0){
    stop("Sets of genes.use and genes.exclude have non-empty intersection")
  }
  
  genes.use <- genes.use[genes.use %in% rownames(normalized_rna)]
  normalized_rna <- normalized_rna[genes.use,]
  
  if(do.approx & approx.N >= N.c){
    do.approx <- FALSE
    warning("approx.N is larger or equal to the number of single cells, thus, an exact simplification will be performed")
  }
  
  # if(do.approx & (approx.N < round(N.c/gamma))){
  #   approx.N <- round(N.c/gamma)
  #   warning(paste("approx.N is set to N.SC", approx.N))
  # }
  # 
  # if(do.approx & ((N.c/gamma) > (approx.N/3))){
  #   warning("approx.N is not much larger than desired number of super-cells, so an approximate simplification may take londer than an exact one!")
  # }
  
  if(do.approx){
    set.seed(seed)
    approx.N            <- min(approx.N, N.c)
    presample           <- sample(1:N.c, size = approx.N, replace = FALSE)
    presampled.cell.ids <- cell.ids[sort(presample)]
    rest.cell.ids       <- setdiff(cell.ids, presampled.cell.ids)
  } else {
    presampled.cell.ids <- cell.ids
    rest.cell.ids       <- c()
  }
  
  normalized_rna.for.pca            <- Matrix::t(normalized_rna[genes.use, presampled.cell.ids])
  if(do.scale){ normalized_rna.for.pca            <- scale(normalized_rna.for.pca) }
  normalized_rna.for.pca[is.na(normalized_rna.for.pca)] <- 0
  
  if(is.null(n.pc[1]) | min(n.pc) < 1){stop("Please, provide a range or a number of components to use: n.pc")}
  if(length(n.pc)==1) n.pc <- 1:n.pc
  
  if(fast.pca & (N.c < 1000)){
    warning("Normal pca is computed because number of cell is low for irlba::irlba()")
    fast.pca <- FALSE
  }
  
  if(!fast.pca){
    PCA.presampled <- tryCatch({
      stats::prcomp(normalized_rna.for.pca, rank. = max(n.pc), scale. = FALSE, center = FALSE)
    }, error = function(e) {
      # Print error message
      cat("Error in PCA computation: ", e$message, "\nExcluding zero variance and retrying...\n")
      
      # Update normalized_rna.for.pca to exclude columns with zero variance
      normalized_rna.for.pca <- normalized_rna.for.pca[, apply(normalized_rna.for.pca, 2, var) != 0]
      
      # Rerun PCA on the updated dataset
      stats::prcomp(normalized_rna.for.pca, rank. = max(n.pc), scale. = FALSE, center = FALSE)
    })
   
  } else {
    set.seed(seed)
    PCA.presampled          <- irlba::irlba(normalized_rna.for.pca, nv = max(n.pc, 25))
    PCA.presampled$x        <- PCA.presampled$u %*% diag(PCA.presampled$d)
    PCA.presampled$rotation <- PCA.presampled$v
  }
  
  
  sc.nw <- SuperCell::build_knn_graph(
    X = PCA.presampled$x[,n.pc],
    k = k.knn, from = "coordinates",
    #use.nn2 = use.nn2,
    dist_method = "euclidean",
    #directed = directed,
    #DoSNN = DoSNN,
    #pruning = pruning,
    #which.snn = which.snn,
    #kmin = kmin,
    ...
  )
  
  list(sc.nw = sc.nw, PCA.presampled = PCA.presampled, 
       normalized_rna.for.pca = normalized_rna.for.pca, 
       presampled.cell.ids = presampled.cell.ids, 
       rest.cell.ids = rest.cell.ids, genes.use = genes.use, cell.ids = cell.ids,
       do.approx = do.approx, n.pc = n.pc, k.knn = k.knn)

}

#' Detection of metacells with the SuperCell approach
#'
#' This function detects metacells (former super-cells) from single-cell gene expression matrix
#'
#'
#' @param preprocessed A list returned by `preprocess_SCimplify` containing preprocessed single-cell data, 
#'        PCA results, and kNN graph.
#' @param cell.annotation a vector of cell type annotation, if provided, metacells that contain single cells of different cell type annotation will be split in multiple pure metacell (may result in slightly larger numbe of metacells than expected with a given gamma)
#' @param cell.split.condition a vector of cell conditions that must not be mixed in one metacell. If provided, metacells will be split in condition-pure metacell (may result in significantly(!) larger number of metacells than expected)
#' @param gamma graining level of data (proportion of number of single cells in the initial dataset to the number of metacells in the final dataset)
#' @param block.size number of cells to map to the nearest metacell at the time (for approx coarse-graining)
#' @param igraph.clustering clustering method to identify metacells (available methods "walktrap" (default) and "louvain" (not recommended, gamma is ignored)).
#' @param return.singlecell.NW whether return single-cell network (which consists of approx.N if \code{"do.approx"} or all cells otherwise)
#' @param return.hierarchical.structure whether return hierarchical structure of metacell
#' @param ... other parameters of \link{build_knn_graph} function
#'
#' @return A tibble with column 'cell' and 'membership' indicating which metacell cluster each cell belongs to.
#' @importFrom igraph cluster_walktrap cluster_louvain contract simplify E V
#' @importFrom Matrix t
#' @importFrom proxy dist
#' @export
postprocess_SCimplify <- function(preprocessed,
                      cell.annotation = NULL,
                      cell.split.condition = NULL,
                      gamma,
                      block.size = 10000,
                      igraph.clustering = c("walktrap", "louvain"),
                      return.singlecell.NW = TRUE,
                      return.hierarchical.structure = TRUE,
                      ...) {
  
  sc.nw = preprocessed$sc.nw
  PCA.presampled = preprocessed$PCA.presampled
  normalized_rna.for.pca = preprocessed$normalized_rna.for.pca
  presampled.cell.ids = preprocessed$presampled.cell.ids
  rest.cell.ids = preprocessed$rest.cell.ids
  genes.use = preprocessed$genes.use
  cell.ids = preprocessed$cell.ids
  do.approx = preprocessed$do.approx
  n.pc = preprocessed$n.pc
  k.knn = preprocessed$k.knn
  #normalized_rna = preprocessed$normalized_rna
  
  N.c <- length(preprocessed$cell.ids)
  
  k   <- round(N.c / gamma)
  
  if (igraph.clustering[1] == "walktrap") {
    g.s              <- igraph::cluster_walktrap(sc.nw$graph.knn)
    g.s$membership   <- igraph::cut_at(g.s, k)
    
  } else if (igraph.clustering[1] == "louvain") {
    warning(paste(
      "igraph.clustering =",
      igraph.clustering,
      ", gamma is ignored"
    ))
    g.s    <- igraph::cluster_louvain(sc.nw$graph.knn)
    
  } else {
    stop(
      paste(
        "Unknown clustering method (",
        igraph.clustering,
        "), please use louvain or walkrtap"
      )
    )
  }
  
  membership.presampled        <- g.s$membership
  names(membership.presampled) <- presampled.cell.ids
  
  ## Split super-cells containing cells from different annotations or conditions
  if (!is.null(cell.annotation) | !is.null(cell.split.condition)) {
    if (is.null(cell.annotation))
      cell.annotation <- rep("a", N.c)
    if (is.null(cell.split.condition))
      cell.split.condition <- rep("s", N.c)
    names(cell.annotation) <- names(cell.split.condition) <- cell.ids
    
    split.cells <- interaction(cell.annotation[presampled.cell.ids], cell.split.condition[presampled.cell.ids], drop = TRUE)
    
    membership.presampled.intr <- interaction(membership.presampled, split.cells, drop = TRUE)
    membership.presampled <- as.numeric(membership.presampled.intr)
    names(membership.presampled) <- presampled.cell.ids
  }
  
  
  
  SC.NW                        <- igraph::contract(sc.nw$graph.knn, membership.presampled)
  if (!do.approx) {
    SC.NW                        <- igraph::simplify(SC.NW,
                                                     remove.loops = T,
                                                     edge.attr.comb = "sum")
  }
  
  
  if (do.approx) {
    PCA.averaged.SC      <- as.matrix(Matrix::t(supercell_GE(t(
      PCA.presampled$x[, n.pc]
    ), groups = membership.presampled)))
    normalized_rna.for.roration       <- Matrix::t(normalized_rna[genes.use, rest.cell.ids])
    
    
    
    if (do.scale) {
      normalized_rna.for.roration <- scale(normalized_rna.for.roration)
    }
    normalized_rna.for.roration[is.na(normalized_rna.for.roration)] <- 0
    
    
    membership.omitted   <- c()
    if (is.null(block.size) | is.na(block.size))
      block.size <- 10000
    
    N.blocks <- length(rest.cell.ids) %/% block.size
    if (length(rest.cell.ids) %% block.size > 0)
      N.blocks <- N.blocks+1
    
    
    if (N.blocks > 0) {
      for (i in 1:N.blocks) {
        # compute knn by blocks
        idx.begin <- (i - 1) * block.size+1
        idx.end   <- min(i * block.size, length(rest.cell.ids))
        
        cur.rest.cell.ids    <- rest.cell.ids[idx.begin:idx.end]
        
        PCA.ommited          <- normalized_rna.for.roration[cur.rest.cell.ids, ] %*% PCA.presampled$rotation[, n.pc] ###
        
        D.omitted.subsampled <- proxy::dist(PCA.ommited, PCA.averaged.SC) ###
        
        membership.omitted.cur        <- apply(D.omitted.subsampled, 1, which.min) ###
        names(membership.omitted.cur) <- cur.rest.cell.ids ###
        
        membership.omitted   <- c(membership.omitted, membership.omitted.cur)
      }
    }
    
    membership.all_       <- c(membership.presampled, membership.omitted)
    membership.all        <- membership.all_
    
    
    names_membership.all <- names(membership.all_)
    ## again split super-cells containing cells from different annotation or split conditions
    if (!is.null(cell.annotation) | !is.null(cell.split.condition)) {
      split.cells <- interaction(cell.annotation[names_membership.all], cell.split.condition[names_membership.all], drop = TRUE)
      
      
      membership.all.intr <- interaction(membership.all_, split.cells, drop = TRUE)
      
      membership.all      <- as.numeric(membership.all.intr)
      
    }
    
    
    SC.NW                        <- igraph::simplify(SC.NW,
                                                     remove.loops = T,
                                                     edge.attr.comb = "sum")
    names(membership.all) <- names_membership.all
    membership.all <- membership.all[cell.ids]
    
  } else {
    membership.all       <- membership.presampled[cell.ids]
  }
  membership       <- membership.all
  
  supercell_size   <- as.vector(table(membership))
  
  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight / 10)
  
  if (igraph::vcount(SC.NW) == length(supercell_size)) {
    igraph::V(SC.NW)$size          <- supercell_size
    igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)
  } else {
    igraph::V(SC.NW)$size          <- as.vector(table(membership.all_))
    igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)
    warning("Supercell graph was not splitted")
  }
  
  res <- list(
    graph.supercells = SC.NW,
    gamma = gamma,
    N.SC = length(unique(membership)),
    membership = membership,
    supercell_size = supercell_size,
    genes.use = genes.use,
    simplification.algo = igraph.clustering[1],
    do.approx = do.approx,
    n.pc = n.pc,
    k.knn = k.knn,
    sc.cell.annotation. = cell.annotation,
    sc.cell.split.condition. = cell.split.condition
  )
  
  if (return.singlecell.NW) {
    res$graph.singlecell <- sc.nw$graph.knn
  }
  if (!is.null(cell.annotation) | !is.null(cell.split.condition)) {
    res$SC.cell.annotation. <- supercell_assign(cell.annotation, res$membership)
    res$SC.cell.split.condition. <- supercell_assign(cell.split.condition, res$membership)
  }
  
  if (igraph.clustering[1] == "walktrap" &
      return.hierarchical.structure)
    res$h_membership <- g.s
  
  metacell_classification <- tibble(cell = res$membership |> names(), 
                                    membership = res$membership)
  
  metacell_classification
}

#' Calculate Appropriate Gamma Values for Metacell Analysis
#'
#' This function determines viable gamma (γ) values to be used in metacell analysis. Gamma is a graining level 
#' parameter that controls the degree of cell aggregation when creating metacells. It represents the ratio 
#' between the original number of cells and the desired number of metacells.
#'
#' For example:
#' - γ = 2: combines cells to create metacells, aiming for half as many metacells as original cells
#' - γ = 4: aims for one-fourth as many metacells
#' - γ = 8: aims for one-eighth as many metacells
#' And so on, using powers of 2.
#'
#' The function starts with γ = 2 and doubles it repeatedly (2, 4, 8, 16...) until the ratio of 
#' cells/gamma would result in metacells that are smaller than the minimum allowed size. Higher gamma 
#' values mean more aggressive aggregation (fewer, larger metacells), while lower gamma values preserve 
#' more granularity (more, smaller metacells).
#'
#' @param cell_count Integer, the total number of cells.
#' @param min_cells_per_metacell Integer, the minimum number of cells allowed per metacell. Defaults to 30.
#' @return An Integer vector of viable gamma values. If no viable gamma values are found, returns 0.
calculate_gamma <- function(cell_count, min_cells_per_metacell = 1) {
  gamma = 2
  gamma_values <- integer()
  while (cell_count / gamma >= min_cells_per_metacell) {
    gamma_values <- c(gamma_values, gamma)
    gamma <- gamma * 2
  }
  if (length(gamma_values) == 0) {
    return(0)  # Return 0 if no viable gamma values
  }
  return(gamma_values)
}

#' Calculate Metacell Membership Scores Across Different Gamma Parameters
#'
#' This function processes single-cell data to identify metacell membership across various gamma settings. 
#' It preprocesses the single-cell data, calculates gamma values based on the number of columns (typically genes), 
#' and postprocesses each gamma setting to assign cells to metacells. It then aggregates these results and 
#' handles missing values by taking the maximum value in each group, ignoring NAs.
#'
#' @param sample_sce a SingleCellExperiment object containing pre-loaded single-cell RNA-seq data.
#' @param min_cells_per_metacell An integer of minimum cells in each metacell.
#' @return A tibble with metacells membership scores across computed gamma settings.
#' @importFrom purrr map
#' @importFrom dplyr rename group_by summarise group_split
#' @examples
#' # Assume 'sce' is a SingleCellExperiment object with a cell type
#' calculate_metacell(sce)
calculate_metacell_for_a_sample_per_cell_type <- function(sample_sce,
                                                          min_cells_per_metacell = 1) {
  # Preprocess the single-cell data
  preprocessed_sce = sample_sce |> preprocess_SCimplify()
  
  # Calculate the number of metacells can be produced
  gammas <- calculate_gamma(sample_sce |> colnames() |> length(),
                            min_cells_per_metacell)
  
  # Postprocess data for each gamma, rename columns, and aggregate results
  gammas |> map(~ postprocess_SCimplify(preprocessed_sce, gamma = .x) |> 
                  dplyr::rename(!!paste0("gamma", .x) :=  membership)) |> 
    bind_rows() |> 
    
    # Group by cell and summarise by taking the max value across all variables, removing NAs
    group_by(cell) |>
    summarise(across(everything(), max, na.rm = TRUE))
}

#' Calculate Metacell Membership for Each Cell Type
#'
#' This function processes a SingleCellExperiment object by grouping cells according
#' to their type, calculates metacell membership for each group, and combines the
#' results into a single tibble.
#'
#' @param sample_sce A SingleCellExperiment object containing single-cell data.
#' @param cell_type_tbl A tibble of cell type.
#' @param empty_droplets_tbl A tibble identifying empty droplets.
#' @param alive_identification_tbl A tibble from alive cell identification.
#' @param doublet_identification_tbl A tibble from doublet identification.
#' @param x A character vector of cell type aggregation column.
#' @param min_cells_per_metacell An integer of minimum cells in each metacell.
#' @return A tibble with metacell membership data for each cell type.
#' @export
split_sample_cell_type_calculate_metacell_membership <- function(sample_sce, 
                                                                 cell_type_tbl,
                                                                 empty_droplets_tbl = NULL,
                                                                 alive_identification_tbl = NULL,
                                                                 doublet_identification_tbl = NULL,
                                                                 x="cell_type",
                                                                 min_cells_per_metacell = NULL) {
  
  if (sample_sce |> is.null()) return(NULL)
  
  if (cell_type_tbl |> is.null()) return(NULL)
  
  # avoid small number of cells 
  if (!is.null(empty_droplets_tbl)) {
    sample_sce <- sample_sce |>
      left_join(empty_droplets_tbl, by = ".cell") |>
      dplyr::filter(!empty_droplet)
  } 
  
  # remove dead cells
  if (!is.null(alive_identification_tbl)) {
    sample_sce =
      sample_sce |>
      left_join(
        alive_identification_tbl ,
        by=".cell"
      ) |> dplyr::filter(alive) 
  }
  
  # remove doublets
  if (!is.null(doublet_identification_tbl)) {
    sample_sce =
      sample_sce |>
      left_join(
        doublet_identification_tbl ,
        by=".cell"
      ) |> dplyr::filter(scDblFinder.class != "doublet") 
  }
  
  # In rare cases, all cells in a sample are from empty droplets or dead or doublets
  if (ncol(sample_sce) == 0) return(NULL)
  
  # # In rare cases, only one cell in a sample is left
  # if (ncol(sample_sce) == 1) sample_sce = sample_sce |> duplicate_single_column_assay()
  
  metacell_gamma_membership_tibble <- sample_sce |> left_join(cell_type_tbl) |> 
    dplyr::group_split(!!sym(x)) |>
    # We need to include all good quality single cells in metacell. 
    # For those cells that cant be halved further, mark metacell_2 to 1
    purrr::map( ~ if (ncol(.x) <= 2) {
       .x |> 
        SummarizedExperiment::colData() |> as.data.frame() |> tibble::rownames_to_column("cell") |> 
        select(cell) |> mutate(gamma2 = 1) |> as_tibble()
      } else if (ncol(.x) >2) {calculate_metacell_for_a_sample_per_cell_type(.x)}) |>
  
      # calculate_metacell_for_a_sample_per_cell_type(.x, min_cells_per_metacell)} else return(NULL)) |> 
    bind_rows()
  
  metacell_gamma_membership_tibble
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
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr left_join  
#' @import SeuratObject
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assay<-
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
                                 cell_type_ensembl_harmonised_tbl = NULL,
                                 annotation_label_transfer_tbl = NULL, 
                                 doublet_identification_tbl = NULL){
  #Fix GCHECKS 
  .cell <- NULL
  alive <- NULL
  subsets_Mito_percent <- NULL
  subsets_Ribo_percent <- NULL
  high_mitochondrion <- NULL
  high_ribosome <- NULL
  scDblFinder.class <- NULL
  predicted.celltype.l2 <- NULL
  
  if (empty_droplets_tbl |> is.null() |> not()) {
    input_read_RNA_assay =
      input_read_RNA_assay |>
      left_join(empty_droplets_tbl, by = ".cell") |>
      filter(!empty_droplet)
  } 
  
  # Add normalisation
  if(!is.null(non_batch_variation_removal_S)){
    if(input_read_RNA_assay |> is("Seurat")) {
      non_batch_variation_removal_S_assay <- CreateAssay5Object(data = non_batch_variation_removal_S)
      input_read_RNA_assay[["SCT"]] <- non_batch_variation_removal_S_assay

    } else if(input_read_RNA_assay |> is("SingleCellExperiment")){
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
    filter(scDblFinder.class!="doublet") 
  
  # attach cell cycle
  if(cell_cycle_score_tbl |> is.null() |> not())
    input_read_RNA_assay = 
    input_read_RNA_assay |>
    left_join(
      cell_cycle_score_tbl ,
      by=".cell"
    )
  
  # Attach annotation
  try({
      if (inherits(annotation_label_transfer_tbl, "tbl_df") && nrow(annotation_label_transfer_tbl) > 0){
        input_read_RNA_assay <- input_read_RNA_assay |>
          left_join(annotation_label_transfer_tbl, by = ".cell") |>
          left_join(cell_type_ensembl_harmonised_tbl)
        
        # Replace NA annotation column with "other", as annotations are single-cell level, not related to pseudobulk
        annotation_columns <- c("blueprint_first.labels.fine",  "blueprint_first.labels.coarse",
                                "monaco_first.labels.fine", "monaco_first.labels.coarse", 
                                "blueprint_first_labels_fine", "monaco_first_labels_fine", 
                                "azimuth_predicted_celltype_l2", "azimuth", "blueprint", "monaco")
        
        cell_type_concensus_columns <- c("cell_type_unified_ensemble", "data_driven_ensemble")
        
        input_read_RNA_assay <- input_read_RNA_assay |> mutate(across(all_of(annotation_columns), 
                                                                      ~tidyr::replace_na(., "Other"))) |> 
          
          mutate(across(all_of(cell_type_concensus_columns), as.character)) |>
          
          # Replace NA with Unknown because non-immune cells are regarded as "Other" in cell_type_unified_ensemble
          mutate(across(all_of(cell_type_concensus_columns), 
                        ~tidyr::replace_na(., "Unknown")))
      }
    }, silent = TRUE)
  
  
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
#' @param container_type A character vector specifying the output file type. Ideally it should match to the input file type.
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
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<-
#' @importFrom digest digest
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' 
#' @export
# Create pseudobulk for each sample 
create_pseudobulk <- function(input_read_RNA_assay, 
                              sample_names_vec, 
                              empty_droplets_tbl = NULL,
                              alive_identification_tbl = NULL,
                              cell_cycle_score_tbl = NULL,
                              annotation_label_transfer_tbl = NULL,
                              cell_type_ensembl_harmonised_tbl = NULL,
                              doublet_identification_tbl = NULL,  
                              x = c() , 
                              external_path, assays = NULL,
                              container_type) {
  #Fix GChecks 
  .sample = NULL 
  .feature = NULL 
  data_source = NULL 
  symbol = NULL 
  
  dir.create(external_path, showWarnings = FALSE, recursive = TRUE)
  
  if (input_read_RNA_assay |> is.null()) return(NULL)
  
  preprocessing_output_S = 
    preprocessing_output(
      input_read_RNA_assay,
      empty_droplets_tbl,
      non_batch_variation_removal_S = NULL, 
      alive_identification_tbl,
      cell_cycle_score_tbl = NULL,
      cell_type_ensembl_harmonised_tbl,
      annotation_label_transfer_tbl, 
      doublet_identification_tbl
    )
  
  # In rare cases, empty droplets observed across cells in a sample
  if (ncol(preprocessing_output_S) == 0) return(NULL)
  
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
    #aggregate_cells(c(sample_hpc, any_of(x)), slot = "data", assays = assays) 
    tidySingleCellExperiment::aggregate_cells(c(sample_hpc, !!sym(x)), slot = "data", assays = assays)
  
  # If I start from Seurat
  if(pseudobulk |> is("data.frame"))
    pseudobulk = pseudobulk |>
    tidybulk::as_SummarizedExperiment(.sample, .feature, any_of(assays)) 
  
  rowData(pseudobulk)$feature_name = rownames(pseudobulk)
  colData(pseudobulk)$pseudobulk_sample = colnames(pseudobulk)
  
  pseudobulk = pseudobulk |>
    pivot_longer(cols = assays, names_to = "data_source", values_to = "count") |>
    filter(!count |> is.na()) |>
    
    # Some manipulation to get unique feature because RNA and ADT
    # both can have same name genes
    dplyr::rename(symbol = .feature) |>
    mutate(data_source = stringr::str_remove(data_source, "abundance_")) |>
    tidyr::unite(".feature", c(symbol, data_source), remove = FALSE) |>
    
    tidybulk::as_SummarizedExperiment(
      .sample = .sample,
      .transcript = .feature,
      .abundance = count
    ) 
  
  # Covert pseudobulk to SCE representation as zellkonverter::writeH5AD 
  #  does not support saving a SummarizedExperiment
  if (container_type == "anndata") {
    pseudobulk = SingleCellExperiment(
      assays = assays(pseudobulk),
      rowData = rowData(pseudobulk),
      colData = colData(pseudobulk)
    )
  }

  file_name = glue::glue("{external_path}/{digest(pseudobulk)}")
  
  # Maybe do not need to save Anndata externally for cellNexus because pseudobulk can be read by tar_read_raw
  pseudobulk |>
    
    # Convert to Anndata
    save_experiment_data(
      dir = file_name, 
      container_type = container_type
    )
  
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

#' Perform Human Cell-Cell Communication Analysis
#' @description This function performs cells communication analysis.
#'   It processes single-cell RNA sequencing data to identify and analyze intercellular communication networks.
#' @param input_read_RNA_assay A SingleCellExperiment or Seurat object containing gene expression data
#' @param empty_droplets_tbl Optional tibble identifying empty droplets to be filtered out
#' @param alive_identification_tbl Optional tibble identifying dead cells to be filtered out
#' @param doublet_identification_tbl Optional A tibble from doublet identification.
#' @param cell_type_tbl Optional A tibble containing cell, cell type, and sample_id information
#' @param assay Character string specifying which assay to use
#' @param cell_type_column Character string specifying the column name containing cell type annotations
#' @param feature_nomenclature Character vector specifying gene in Symbol or Ensemble format
#' @param reference_db The ligand-receptor interaction database curated in CellChat tool. Choose between human or mouse.
#' @param ... Additional arguments passed to \code{CellChat::subsetDB}
#' @return A CellChat tibble containing the inferred communication at the level of 
#'     ligands/receptors
#' @importFrom CellChat createCellChat subsetDB subsetData identifyOverExpressedGenes 
#'   identifyOverExpressedInteractions computeCommunProb filterCommunication subsetCommunication
#'   normalizeData smoothData aggregateNet setIdent
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter mutate
#' @export
cell_communication <- function(input_read_RNA_assay,
                               empty_droplets_tbl = NULL,
                               alive_identification_tbl = NULL,
                               doublet_identification_tbl = NULL,
                               cell_type_tbl = NULL,
                               assay = NULL,
                               cell_type_column = NULL,
                               feature_nomenclature,
                               reference_db = "human",
                               ...){
  
  # Input should not be NULL
  if (is.null(input_read_RNA_assay)) return(NULL)
  
  # Get assay
  if(is.null(assay)) my_assay = input_read_RNA_assay@assays |> names() |> magrittr::extract2(1)
  
  # Identify cell type column
  if(
    is.null(cell_type_column) && 
    !cell_type_column %in% colnames(as_tibble(input_read_RNA_assay[1,1]))
  ) stop("HPCell says: Your `cell_type_column` columns are not present in your data. Please run celltype_consensus_constructor() to get the cell type annotation that you can use as grouping.")
  
  # Avoid small number of cells 
  if (!is.null(empty_droplets_tbl)) {
    input_read_RNA_assay <- input_read_RNA_assay |>
      left_join(empty_droplets_tbl, by = ".cell") |>
      dplyr::filter(!empty_droplet)
  } 
  
  # Avoid dead cells
  if (!is.null(alive_identification_tbl)) {
    input_read_RNA_assay <- input_read_RNA_assay |>
      left_join(alive_identification_tbl |> select(.cell, alive), by = ".cell") |>
      dplyr::filter(alive)
  } 
  
  # Avoid doublet
  if (!is.null(doublet_identification_tbl)) {
    input_read_RNA_assay <- input_read_RNA_assay |> 
      left_join(doublet_identification_tbl |> select(.cell, scDblFinder.class), by = ".cell") |>
      filter(scDblFinder.class!="doublet") 
  }
  
  # Append cell type
  if (!is.null(cell_type_tbl) && 
      cell_type_column %in% colnames(cell_type_tbl)) {
    input_read_RNA_assay <- input_read_RNA_assay |>
      left_join(cell_type_tbl) |> 
      select(everything(), !!cell_type_column) |> 
      filter(!is.na(.data[[cell_type_column]]))
  } else if (!is.null(cell_type_tbl) && !cell_type_column %in% colnames(cell_type_tbl))
      stop ("HPCell says: Your `cell_type_column` does not present in `cell_type_tbl` data")
  
  # Note: CellChat only takes gene symbols as input, thus conversion step is required for ensemble IDs
  if (feature_nomenclature == "ensembl") {
    
    gene_map = rownames(input_read_RNA_assay) |> convert_gene_names(current_nomenclature = feature_nomenclature) |> 
      filter(!is.na(gene_name))
    
    input_read_RNA_assay = input_read_RNA_assay[gene_map$gene_id, ]
    rownames(input_read_RNA_assay) = gene_map$gene_name
  }
  
  if (inherits(input_read_RNA_assay, "SingleCellExperiment")) {
    counts = input_read_RNA_assay |> assay(my_assay)
    meta = input_read_RNA_assay |> colData() |> as.data.frame()
  } else if (inherits(input_read_RNA_assay, "Seurat")){
    counts = GetAssayData(input_read_RNA_assay, assay = my_assay) 
    meta = input_read_RNA_assay[[]]
  }
  
  if (meta |> nrow() == 0) return(NULL)
  
  # CellChat identifyOverExpressedGenes() would only support at least 2 groups
  if (meta |> distinct(.data[[cell_type_column]]) |> pull() |> length() <= 1) return(NULL)

  # Choose cellchat reference
  DB <- switch(reference_db, human = CellChat::CellChatDB.human, mouse = CellChat::CellChatDB.mouse)
  projectionDB <- switch(reference_db, human = CellChat::PPI.human, mouse = CellChat::PPI.mouse)
    
  CellChatDB <- DB
  
  CellChatDB.use <- subsetDB(CellChatDB, search =c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), key = c("annotation")) 
  
  # CellChat Only Takes log-Normalized data
  cellchat = counts |> 
    as("dgCMatrix") |> 
    normalizeData(do.log = TRUE) |>
    createCellChat(group.by = cell_type_column, meta = meta, assay = my_assay)
  
  cellchat@DB <- CellChatDB.use
  
  # Preprocessing
  cellchat  = cellchat |> 
    subsetData() |> 
    identifyOverExpressedGenes() |> 
    identifyOverExpressedInteractions() |> 
    smoothData(adj = projectionDB)
  
  # Return NULL when none of LR pairs are found
  if (nrow(cellchat@LR$LRsig) == 0) return(NULL)
    
  cellchat = cellchat |> 
    
    # Use projected data
    computeCommunProb(raw.use = FALSE) |> 
    
    # Filter the number of cells in each group are less than 10
    filterCommunication(min.cells = 10) |> 
    computeCommunProbPathway() |>
    aggregateNet()
  
  gc()
  
  # Extract the inferred cellular communication network as a data frame
  # By default, slot.name = "net" extracts the inferred communication at the level of ligands/receptors
  # Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
  # If all arguments are NULL, it returns a data frame consisting of all the inferred cell-cell communications
  lr_tbl <- tryCatch(
    subsetCommunication(cellchat, slot.name = "net", thresh = NULL) |> 
      dplyr::rename(lr_prob = prob,
                    lr_pval = pval),
    error = function(e) {
      message("Error in subsetCommunication(): ", e$message)
      return(NULL)
    }
  )
  
  pathway_tbl <- tryCatch(
    subsetCommunication(cellchat, slot.name = "netP", thresh = NULL) |> 
      dplyr::rename(pathway_prob = prob,
                    pathway_pval = pval),
    error = function(e) {
      message("Error in subsetCommunication(): ", e$message)
      return(NULL)
    }
  )
  
  if (nrow(lr_tbl) == 0 || is.null(lr_tbl)) return(NULL)
  
  cell_interaction_count = cellchat@net$count |> as_tibble(rownames = "source") |>
    pivot_longer(-source, names_to = "target", values_to = "interaction_count") 
  
  cell_interaction_weight =  cellchat@net$weight |> as_tibble(rownames = "source") |>
    pivot_longer(-source, names_to = "target", values_to = "interaction_weight")
  
  result = lr_tbl |> left_join(pathway_tbl,  by = c("source", "target", "pathway_name")) |>
      mutate(sample_id = unique(cellchat@meta$sample_id)) |>
      left_join(cell_interaction_count) |> left_join(cell_interaction_weight) |>
      as_tibble()
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
