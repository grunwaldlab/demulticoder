% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/format_databases.R
\name{format_db_16s}
\alias{format_db_16s}
\title{An 16s database that has modified headers and is output in the
reference_databases folder}
\usage{
format_db_16s(
  data_tables,
  data_path,
  output_directory_path,
  temp_directory_path,
  db_16s
)
}
\arguments{
\item{data_tables}{The data tables containing the paths to read files, metadata, primer sequences}

\item{data_path}{Path to the data directory}

\item{output_directory_path}{The path to the directory where resulting files
are output}

\item{temp_directory_path}{User-defined temporary directory to place reads
throughout the workflow metadata, and primer_info files}

\item{db_16s}{The name of the database}
}
\value{
An 16s database that has modified headers and is output in the
reference_databases folder
}
\description{
An 16s database that has modified headers and is output in the
reference_databases folder
}
\keyword{internal}
