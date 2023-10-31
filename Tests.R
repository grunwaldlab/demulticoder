devtools::load_all("/Users/masudermann/rps10package")
devtools::document()
#devtools::install
#devtools::build_rmd("vignettes/Introduction.Rmd")
#browseVignettes('rps10package')

#TODO-make separate data directory vs. working directory-otherwise test folder gets messy
#Revise instruction and documentation so more consistent and clearer
#If everything comes together finally add functionality for other barcodes. How hard to integrate ITSx? 

prepare_reads(
  maxN = 0, 
  data_directory = "inst/extdata", 
  output_directory = "~/output_package2", 
  tempdir_id = "run1")

cut_trim(
  analysis_setup,
  cutadapt_path="/opt/homebrew/bin/cutadapt",
  verbose = TRUE,
  maxEE = 2,
  truncQ = 5,
  minLen = 200,
  maxLen = 297,
  minCutadaptlength = 50)

make_asv_abund_matrix(
  analysis_setup,
  minOverlap = 15,
  maxMismatch = 2,
  verbose = TRUE,
  multithread = TRUE)

assignTax(
  analysis_setup,
  asv_abund_matrix,
  multithread = TRUE,
  barcode = "rps10", 
  retrieve_files=TRUE
)


asv_matrix_to_taxmap(asv_abund_matrix,
                     min_read_depth = 10,
                     minimum_bootstrap = 0)

taxmap_to_phyloseq()


library(testthat)
check()

#R CMD build

#usethis::use_readme_rmd()
#devtools::build_readme()
#usethis::use_test("rps10package2")
#usethis::use_test("PrepareReads_CountPrimers")
#usethis::use_vignette("Package_Workflow")


