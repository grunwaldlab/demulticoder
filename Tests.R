devtools::load_all("/Users/masudermann/rps10package3")
devtools::document()
#devtools::install
#devtools::build_rmd("vignettes/Introduction.Rmd")
#browseVignettes('rps10package')

#TODO-make separate data directory vs. working directory-otherwise test folder gets messy
#Revise instruction and documentation so more consistent and clearer
#If everything comes together finally add functionality for other barcodes. How hard to integrate ITSx? 

outputs <- prepare_reads(maxN = 0, data_directory = "inst/extdata", output_directory = "~/output_package", tempdir_id = "run2")

#directory_path<-"inst/extdata" ##choose a directory for all downstream steps
#directory_path_temp <- file.path(tempdir(), paste0("run10_", Sys.Date()))
#dir.create(directory_path_temp)
#primer_path <-file.path(directory_path, "primer_info.csv") ##modify .csv name or keep this name
#metadata_path <-file.path(directory_path,"metadata.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
#cutadapt_path<-"/opt/homebrew/bin/cutadapt"

#dir_paths <- 
#  setup_directories(
#    data_directory = "inst/extdata",
#    output_directory = "~/outputs",
#    tempdir_id = "run1"
#  )


#data_tables <- prepare_reads(
#  directory_path = dirs$data_directory,
#  directory_path_temp = dirs$temp_directory,
#  primer_path = dirs$primer_path,
#  metadata_path = dirs$metadata_path,
#  maxN = 0,
#  multithread = TRUE,
#  overwrite_existing = TRUE
#)



cut_trim(
  directory_path,
  directory_path_temp,
  cutadapt_path,
  verbose = TRUE,
  maxEE = 2,
  truncQ = 5,
  minLen = 200,
  maxLen = 297,
  minCutadaptlength = 50
)

asv_abund_matrix <-
  make_asv_abund_matrix(
    directory_path,
    directory_path_temp,
    minOverlap = 15,
    maxMismatch = 2,
    verbose = TRUE,
    multithread = TRUE
  )

summary <- assignTax(
  directory_path,
  directory_path_temp,
  data_tables,
  asv_abund_matrix,
  multithread = TRUE,
  barcode = "rps10", 
  retrieve_files=TRUE
)

#retrieve reads if user would like

obj_dada <-
  asvmatrix_to_taxmap(asv_abund_matrix,
                      min_read_depth = 10,
                      minimum_bootstrap = 75)

phylo_obj <- taxmap_to_phyloseq(obj_dada)


library(testthat)

check()

#R CMD build

#usethis::use_readme_rmd()
#devtools::build_readme()
#usethis::use_test("rps10package2")
#usethis::use_test("PrepareReads_CountPrimers")
#usethis::use_vignette("Package_Workflow")


