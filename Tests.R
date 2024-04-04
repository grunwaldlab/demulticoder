devtools::load_all("~/demulticoder")
devtools::document()

prepare_reads(
  maxN = 0, 
  data_directory = "~/demulticoder/inst/extdata", 
  output_directory = "~/testing_package10", 
  tempdir_id = "tempdir",
  tempdir_path="~/",
  overwrite_existing = FALSE)

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
  db_16s = "bacteriadb.fasta",
  retrieve_files=TRUE,
  overwrite_existing=TRUE)

convert_asv_matrix_to_objs(save_outputs=TRUE, overwrite=TRUE)

load("~/demulticoder/inst/extdata/obj_dada_its.RData") 
obj_dada_its<-obj_dada


library("pkgdown")
library(pkgdown)


pkgdown::build_site("~/demulticoder")
