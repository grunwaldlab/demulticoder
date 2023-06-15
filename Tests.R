devtools::load_all("~/rps10package")
document()

directory_path<-"~/rps10package/raw_data/rps10_ITS" ##choose a directory for all downstream steps
primer_path <-file.path(directory_path, "primer_info.csv") ##modify .csv name or keep this name
metadata_path <-file.path(directory_path,"metadata.csv") ##modify .csv name or keep this name. The sample_name in the metadata sheet needs to match the first part (before first underscore), of the zipped raw FASTQ files
cutadapt_path<-"/opt/homebrew/bin/cutadapt"

data_tables <-
  prepare_reads(
    directory_path,
    primer_path,
    metadata_path,
    maxN = 0,
    multithread = TRUE
  )

cut_trim(
  directory_path,
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
    minOverlap = 15,
    maxMismatch = 2,
    verbose = TRUE,
    multithread = TRUE
  )

summary <- assignTax(
  directory_path,
  data_tables,
  asv_abund_matrix,
  multithread = TRUE,
  barcode = "rps10_its",
)

obj_dada <-
  asvmatrix_to_taxmap(asv_abund_matrix,
                      min_read_depth = 10,
                      minimum_bootstrap = 75)

phylo_obj <- taxmap_to_phyloseq(obj_dada)


library(testthat)

check()

#R CMD build

usethis::use_readme_rmd()
devtools::build_readme()

usethis::use_test("rps10package2")
usethis::use_test("PrepareReads_CountPrimers")


devtools::build_vignettes()
