% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cut_primers_trim_reads.R
\name{run_cutadapt}
\alias{run_cutadapt}
\title{Core function for running cutadapt}
\usage{
run_cutadapt(
  cutadapt_path,
  cutadapt_data_barcode,
  barcode_params,
  minCutadaptlength
)
}
\arguments{
\item{cutadapt_path}{A path to the cutadapt program.}

\item{minCutadaptlength}{Read lengths that are lower than this threshold will
be discarded. Default is 50.}

\item{cutadapt_data}{Directory_data folder with trimmed and filtered reads
for each sample.}
}
\value{
Trimmed read.
}
\description{
Core function for running cutadapt
}
\keyword{internal}
