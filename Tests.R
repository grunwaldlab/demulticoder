devtools::load_all("/Users/masudermann/demulticoder")
devtools::document()

prepare_reads(
  maxN = 0, 
  data_directory = "~/demulticoder/inst/extdata", 
  output_directory = "~/testing_package3", 
  tempdir_id = "test10",
  tempdir_path="~/",
  overwrite_existing = FALSE)

cut_trim(
  analysis_setup,
  cutadapt_path="/opt/homebrew/bin/cutadapt",
  overwrite_existing = FALSE)

make_asv_abund_matrix(
  analysis_setup,
  overwrite_existing = FALSE)

assign_tax(
  analysis_setup,
  asv_abund_matrix,
  retrieve_files=TRUE,
  overwrite_existing=TRUE)

convert_asv_matrix_to_objs(save_outputs=TRUE, overwrite=TRUE)

library(testthat)
check()