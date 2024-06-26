% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{create_pseudobulk}
\alias{create_pseudobulk}
\title{Create pseudobulk}
\usage{
create_pseudobulk(preprocessing_output_S, sample_names, x, ...)
}
\arguments{
\item{preprocessing_output_S}{Processed dataset from preprocessing.}

\item{x}{A grouping variable used to aggregate cells into pseudobulk samples.
This variable should be present in the \code{preprocessing_output_S} object and
typically represents a factor such as sample ID or condition.}

\item{...}{Additional arguments passed to internal functions used within
\code{create_pseudobulk}. This includes parameters for customization of
aggregation, data transformation, or any other process involved in the
creation of pseudobulk samples.}

\item{assays}{A character vector specifying the assays to be included in the
pseudobulk creation process, such as c("RNA", "ADT").}
}
\value{
List containing pseudobulk data aggregated by sample and by both sample and cell type.
}
\description{
Aggregates cells based on sample and cell type annotations, creating pseudobulk samples
for each combination. Handles RNA and ADT assays
}
