% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assign_taxonomy.R
\name{process_single_barcode}
\alias{process_single_barcode}
\title{Run 'DADA2' taxonomy functions for single metabarcode}
\usage{
process_single_barcode(
  data_tables,
  temp_directory_path,
  output_directory_path,
  asv_abund_matrix,
  metabarcode = metabarcode,
  barcode_params
)
}
\arguments{
\item{data_tables}{The data tables containing the paths to read files, metadata, and metabarcode information with associated primer sequences}

\item{asv_abund_matrix}{The final abundance matrix containing amplified sequence variants}
}
\description{
Run 'DADA2' taxonomy functions for single metabarcode
}
\keyword{internal}
