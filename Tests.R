devtools::load_all("/Users/masudermann/rps10package")
devtools::document()
#devtools::install
#devtools::build_rmd("vignettes/Introduction.Rmd")
#browseVignettes('rps10package')

prepare_reads(
  maxN = 0, 
  data_directory = "inst/extdata", 
  output_directory = "~/rps10_output", 
  tempdir_id = "run1",
  overwrite_existing = FALSE)

cut_trim(
  analysis_setup,
  cutadapt_path="/opt/homebrew/bin/cutadapt",
  overwrite_existing = FALSE)


make_asv_abund_matrix(
  analysis_setup,
  verbose = TRUE,
  overwrite_existing = FALSE)

assignTax(
  analysis_setup,
  asv_abund_matrix,
  barcode = "rps10_its", 
  retrieve_files=FALSE,
  its_db="sh_general_release_dynamic_29.11.2022.fasta",
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


