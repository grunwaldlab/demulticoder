% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_taxonomy.R
\name{assign_taxonomyDada2}
\alias{assign_taxonomyDada2}
\title{Assign taxonomy}
\usage{
assign_taxonomyDada2(
  asv_abund_matrix,
  temp_directory_path,
  minBoot = 0,
  tryRC = FALSE,
  verbose = FALSE,
  multithread = TRUE,
  locus = barcode
)
}
\arguments{
\item{asv_abund_matrix}{The ASV abundance matrix}

\item{minBoot}{(Optional). Default 50. 
The minimum bootstrap confidence for assigning a taxonomic level.}

\item{tryRC}{(Optional). Default FALSE. 
If TRUE, the reverse-complement of each sequences will be used for classification if it is a better match to the reference
sequences than the forward sequence.}

\item{verbose}{(Optional). Default FALSE.
If TRUE, print status to standard output.}

\item{multithread}{(Optional). Default is FALSE.
If TRUE, multithreading is enabled and the number of available threads is automatically determined.   
If an integer is provided, the number of threads to use is set by passing the argument on to
\code{\link{setThreadOptions}}.}

\item{ref_database}{The reference database used for taxonomic inference steps}
}
\description{
Assign taxonomy
}
\keyword{internal}
