#' Dummy HPC Dataset
#'
#' This dataset, named `dummy_hpc`, is a synthetic example dataset used for demonstrating
#' high-performance computing (HPC) data processing and analysis techniques.
#'
#' @format A data frame with multiple rows and columns representing synthetic HPC data. Each column represents a different variable in the dataset.
#'
#' @details
#' The `dummy_hpc` dataset is created for educational and demonstration purposes. It includes
#' simulated data points that resemble typical HPC workload metrics. The dataset can be used to
#' showcase various data processing, analysis, and visualisation techniques in the context of HPC.
#'
#' @usage
#' data(dummy_hpc)
#'
#' @examples
#' # Load the dataset
#' data(dummy_hpc)
#'
#' # Display the first few rows of the dataset
#' head(dummy_hpc)
#'
#' # Example analysis: summary statistics
#' summary(dummy_hpc)
#'
#' @keywords datasets
#' 
#' @noRd
#' 
"dummy_hpc"


#' A data frame of Ensembl genes retrieved from biomaRt package
#' 
#' This dataset contains Ensembl gene IDs, external gene names, and chromosome names
#' retrieved using the biomaRt package.
#' 
#' @format A data frame map of ensembl_gene_id, external_gene_name and chromosome_name
#' 
#' @usage 
#' data(ensembl_genes_biomart)
#' 
#' @source biomaRt::getBM()
#' 
#' @keywords datasets
#' @docType data
"ensembl_genes_biomart"