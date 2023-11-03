devtools::load_all("/Users/masudermann/rps10package")
devtools::document()
#devtools::install
#devtools::build_rmd("vignettes/Introduction.Rmd")
#browseVignettes('rps10package')

prepare_reads(
  maxN = 0, 
  data_directory = "inst/extdata", 
  output_directory = "~/output_package2", 
  tempdir_id = "run5",
  overwrite_existing = FALSE)

cut_trim(
  analysis_setup,
  cutadapt_path="/opt/homebrew/bin/cutadapt",
  verbose = TRUE,
  #maxEE = 2,
  #truncQ = 5,
  #minLen = 200,
  #maxLen = 297,
  #minCutadaptlength = 50, 
  overwrite_existing = FALSE)


make_asv_abund_matrix(
  analysis_setup,
  #minOverlap = 15,
  #maxMismatch = 2,
  verbose = TRUE,
  multithread = TRUE,
  overwrite_existing = FALSE)

assignTax(
  analysis_setup,
  asv_abund_matrix,
  multithread = TRUE,
  barcode = "rps10_its", 
  retrieve_files=TRUE,
  overwrite_existing=FALSE)


asv_matrix_to_taxmap(save_taxmap = TRUE)
taxmap_to_phyloseq(save_phyloseq=TRUE)


library(testthat)
check()

#R CMD build

#usethis::use_readme_rmd()
#devtools::build_readme()
#usethis::use_test("rps10package2")
#usethis::use_test("PrepareReads_CountPrimers")
#usethis::use_vignette("Package_Workflow")


