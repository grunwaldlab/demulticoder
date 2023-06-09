% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Trim_functions.R
\name{run_cutadapt}
\alias{run_cutadapt}
\title{Core function for running cutadapt}
\usage{
run_cutadapt(cutadapt_path, cutadapt_data, minCutadaptlength = 20)
}
\arguments{
\item{cutadapt_path}{A path to the cutadapt program}

\item{cutadapt_data}{directory_data folder with trimmed and filtered reads for each sample}

\item{minCutadaptlength}{Read lengths that are lower than this threshold will be discarded. Default is 50.}
}
\description{
Core function for running cutadapt
}