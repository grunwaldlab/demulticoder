% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_asv_abund_matrix.R
\name{get_fastq_paths}
\alias{get_fastq_paths}
\title{Retrieve the paths of the filtered and trimmed Fastq files}
\usage{
get_fastq_paths(data_tables, my_direction, my_primer_pair_id)
}
\arguments{
\item{data_tables}{The data tables containing the paths to read files, metadata, and metabarcode information with associated primer sequences}

\item{my_direction}{Whether primer is in forward or reverse direction}

\item{my_primer_pair_id}{The specific metabarcode ID}
}
\description{
Retrieve the paths of the filtered and trimmed Fastq files
}
\keyword{internal}
