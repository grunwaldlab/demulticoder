devtools::load_all("/Users/masudermann/rps10package")
devtools::document()
sessionInfo()


#devtools::install
#devtools::build_rmd("vignettes/Introduction.Rmd")
#browseVignettes('rps10package')

prepare_reads(
  maxN = 0, 
  data_directory = "inst/extdata", 
  output_directory = "~/results_test1", 
  tempdir_id = "run3",
  overwrite_existing = TRUE)

cut_trim(
  analysis_setup,
  cutadapt_path="/opt/homebrew/bin/cutadapt",
  maxEE=8, 
  truncQ=5,
  minCutadaptlength = 50,
  multithread=TRUE,
  overwrite_existing = TRUE)

make_asv_abund_matrix(
  analysis_setup,
  verbose = TRUE,
  maxMismatch = 2,
  minOverlap = 15,
  overwrite_existing = TRUE)

assignTax(
  analysis_setup,
  asv_abund_matrix,
  barcode = "rps10_its",
  retrieve_files=TRUE,
  overwrite_existing=TRUE)


asv_matrix_to_taxmap(save_taxmap = TRUE)
taxmap_to_phyloseq(save_phyloseq=TRUE)


summary_table <- process_single_barcode(analysis_setup$data_tables, analysis_setup$dir_paths$temp_directory, analysis_setup$dir_paths$output_directory, asv_abund_matrix, barcode = "rps10")
assign("summary_table", summary_table, envir = .GlobalEnv)

library(testthat)
check()

#R CMD build

#usethis::use_readme_rmd()
#devtools::build_readme()
#usethis::use_test("rps10package2")
#usethis::use_test("PrepareReads_CountPrimers")
#usethis::use_vignette("Package_Workflow")


#To do-16S functionality, and if so any permutation of two barcodes
#Would need to resolve issues with cutandtrim-being able to adjust the parameters for barcodes individually
#IDtaxa-requires specific db formats but this could solve some issues as well, and maybe make others? 