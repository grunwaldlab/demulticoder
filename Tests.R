devtools::load_all("/Users/masudermann/demulticoder")
devtools::document()

prepare_reads(
  maxN = 0, 
  data_directory = "~/demulticoder/inst/extdata", 
  output_directory = "~/testing_package5", 
  tempdir_id = "test11",
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
set.seed(2)
assign_tax(
  analysis_setup,
  asv_abund_matrix,
  tryRC = FALSE,
  db_its = "fungidb_larger3.fasta",
  db_rps10 = "oomycetedb.fasta",
  db_16s = "bacteriadb.fasta",
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

load("~/testing_package4/test11/asvabund_matrixDADA2_rps10.RData")
set.seed(2)
taxa<-assignTaxonomy(
  asv_abund_matrix,
  "~/demulticoder/inst/extdata/oomycetedb.fasta",
  minBoot = 0,
  tryRC = FALSE,
  outputBootstraps = TRUE,
  taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  multithread = TRUE,
  verbose = TRUE
)

set.seed(2)
load("~/testing_package5/test11/")
taxa<-assignTaxonomy(
  asv_abund_matrix,
  "~/test11/sixteenS_reference_db.fa",
  minBoot = 0,
  tryRC = FALSE,
  outputBootstraps = TRUE,
  #taxLevels = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Reference"),
  multithread = TRUE,
  verbose = TRUE
)
write.table(taxa,file="~/test_rps10original2.out",sep="\t",quote=F)
