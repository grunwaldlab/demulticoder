#TODO fix up documentation 

devtools::load_all("~/demulticoder")
devtools::document()

prepare_reads(
  data_directory = "~/demulticoder/inst/extdata", 
  output_directory = "~/testing_package11", 
  tempdir_id = "temp",
  tempdir_path = "~/",
  overwrite_existing = TRUE)

cut_trim(
  analysis_setup,
  cutadapt_path="/opt/homebrew/bin/cutadapt",
  overwrite_existing = TRUE)

make_asv_abund_matrix(
  analysis_setup,
  overwrite_existing = TRUE)

#still need to fix align file
assign_tax(
  analysis_setup,
  asv_abund_matrix,
  db_its = "fungidb.fasta",
  db_rps10 = "oomycetedb.fasta",
  retrieve_files=TRUE,
  overwrite_existing=TRUE)

convert_asv_matrix_to_objs(analysis_setup,save_outputs=TRUE,overwrite=TRUE)

load("~/demulticoder/inst/extdata/obj_dada_its.RData") 
obj_dada_its<-obj_dada

devtools::check()
library("pkgdown")


pkgdown::build_site("./")

install.packages("testthat")
library(testthat)

usethis::use_testthat(3)

library(testthat)
library(demulticoder)

testthat::use_test()
devtools::test()

test_check("demulticoder")
