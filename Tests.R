devtools::load_all("~/demulticoder")
devtools::document()

prepare_reads(
  maxN = 0, 
  data_directory = "~/demulticoder/inst/extdata", 
  output_directory = "~/testing_package10", 
  tempdir_id = "tempdir",
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
assign_tax(
  analysis_setup,
  asv_abund_matrix,
  db_its = "fungidb_larger3.fasta",
  db_rps10 = "oomycetedb.fasta",
  db_16s = "bacteriadb.fasta",
  retrieve_files=TRUE,
  overwrite_existing=TRUE)

convert_asv_matrix_to_objs(save_outputs=TRUE, overwrite=TRUE)

load("~/demulticoder/inst/extdata/obj_dada_its.RData") 
obj_dada_its<-obj_dada


heat_tree(obj_dada_its,
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs,
          node_color_axis_label = "ASV count",
          node_size_axis_label = "Total Abundance of Taxa",
          layout = "da", initial_layout = "re")


#TODO-cleanup functions again and revise documentation
#fix align file
#IDtaxa implementation for silva or unite?
#BLAST addition? 

library(testthat)
check()

usethis::use_pkgdown()

load("~/testing_package8/test14/asvabund_matrixDADA2_rps10.RData")
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

set.seed(1)
load("~/testing_package10/tempdir/asvabund_matrixDADA2_sixteenS.RData")
taxa<-assignTaxonomy(
  asv_abund_matrix,
  "~/bacteriadb.fasta",
  minBoot = 0,
  tryRC = TRUE,
  outputBootstraps = TRUE,
  multithread = TRUE,
  verbose = TRUE
)


load("~/testing_package10/tempdir/TaxMatrix_rps10.RData")
load("~/testing_package10/test14/.RData")

write.table(taxa,file="~/soxteemS_norename.out",sep="\t",quote=F)
write.table(tax_results,file="~/rps10_NAs_refdb.out",sep="\t",quote=F)
