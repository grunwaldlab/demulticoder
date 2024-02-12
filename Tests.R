devtools::load_all("/Users/masudermann/demulticoder")
devtools::document()

prepare_reads(
  maxN = 0, 
  data_directory = "~/demulticoder/inst/extdata", 
  output_directory = "~/testing_package4", 
  tempdir_id = "test11",
  tempdir_path="~/",
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
  retrieve_files=TRUE,
  overwrite_existing=TRUE)

convert_asv_matrix_to_objs(save_outputs=TRUE, overwrite=TRUE)

#TODO-cleanup functions again and revise documentation
#fix align file
#IDtaxa implementation for silva or unite?
#BLAST addition? 

library(testthat)
check()

usethis::use_pkgdown()
