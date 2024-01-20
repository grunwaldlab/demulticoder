devtools::load_all("/Users/masudermann/demulticoder")
devtools::document()

#install.packages("devtools")
#devtools::install_github("grunwaldlab/demulticoder")

#sessionInfo()
#devtools::install
#devtools::build_rmd("vignettes/Introduction.Rmd")
#browseVignettes('rps10package')

prepare_reads(
  maxN = 0, 
  data_directory = "~/demulticoder/inst/extdata", 
  output_directory = "~/testing_package", 
  tempdir_id = "run1",
  overwrite_existing = TRUE)

cut_trim(
  analysis_setup,
  cutadapt_path="/opt/homebrew/bin/cutadapt",
  overwrite_existing = TRUE)

make_asv_abund_matrix(
  analysis_setup,
  overwrite_existing = TRUE)

assignTax(
  analysis_setup,
  asv_abund_matrix,
  retrieve_files=TRUE,
  overwrite_existing=TRUE)

asv_matrix_to_taxmap_phyloseq(save_outputs=TRUE)

library(testthat)
check()

#R CMD build

#usethis::use_readme_rmd()
#devtools::build_readme()
#usethis::use_test("rps10package2")
#usethis::use_test("PrepareReads_CountPrimers")
#usethis::use_vignette("Package_Workflow")


#Need to be clear for ASV inference step, for DADA2 needs greater than 1 sample!!!
#make sure input requirements are clear and there are templates and examples
#IDtaxa-requires specific db formats but this could solve some issues as well, and maybe make others? 