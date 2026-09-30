\name{evaluate_hpc}
\alias{evaluate_hpc}
\alias{evaluate_hpc.HPCell}
\title{Execute the HPCell Targets Pipeline}
\usage{
evaluate_hpc(input_hpc)

\method{evaluate_hpc}{HPCell}(input_hpc)
}
\arguments{
\item{input_hpc}{An \code{HPCell} object constructed by \code{initialise_hpc()} and
extended with one or more pipeline step functions.}
}
\value{
Either the assembled single-cell object (if \code{get_single_cell()} was
called) or a \code{tibble} with targets metadata.
}
\description{
Closes the pipeline target list and calls \code{targets::tar_make()} to execute
all queued steps. If \code{get_single_cell()} was the last step, returns the
assembled single-cell object; otherwise returns the \code{tar_meta()} table.
}
