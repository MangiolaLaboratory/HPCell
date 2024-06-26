% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{is_strong_evidence}
\alias{is_strong_evidence}
\title{Check for Strong Evidence}
\usage{
is_strong_evidence(
  single_cell_data,
  cell_annotation_azimuth_l2,
  cell_annotation_blueprint_singler
)
}
\arguments{
\item{single_cell_data}{A data frame containing single-cell data.}

\item{cell_annotation_azimuth_l2}{A column representing Azimuth L2 cell annotation.}

\item{cell_annotation_blueprint_singler}{A column representing Blueprint Singler cell annotation.}
}
\value{
A data frame with a column indicating strong evidence.
}
\description{
This function checks for strong evidence in cell annotations.
}
