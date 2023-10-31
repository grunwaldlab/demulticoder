#' Filter ASV abundance matrix and convert to taxmap object
#'
#' @param asv_abund_matrix The final ASV abundance matrix
#' @param min_read_depth Threshold for ASVs to remove if number of
#' reads is less than this value across all samples-todo check on this
#' @param minimum_bootstrap Threshold for bootstrap support value
#' for taxonomic assignments. Below designated minimum bootstrap
#' threshold, taxnomoic assignments will be set to N/A.
#' @param pid_species If percent identity is below this value, at the species
#' level, taxonomic assignment will be set to N/A
#' @param pid_genus If percent identity is below this value, at the genus
#' level, taxonomic assignment will be set to N/A
#' @param pid_family If percent identity is below this value, at the family
#' level, taxonomic assignment will be set to N/A
#' @return ASV matrix converted to taxmap object
#' @export asv_matrix_to_taxmap
#' @examples
#' directory_path<-"~/Desktop/rps10test"
#' primer_path <-file.path(directory_path, "primer_info.csv")
#' metadata_path <-file.path(directory_path,"metadata.csv")
#' cutadapt_path<-"/opt/homebrew/bin/cutadapt"
#' data_tables <-
#' prepare_reads(
#' directory_path,
#' primer_path,
#' metadata_path,
#' maxN = 0,
#' )
#' cut_trim(
#' directory_path,
#' cutadapt_path,
#' verbose = TRUE,
#' maxEE = 2,
#' truncQ = 5,
#' minLen = 200,
#' maxLen = 297,
#' minCutadaptlength = 50
#')
#' asv_abund_matrix <-
#' make_asv_abund_matrix(
#' directory_path,
#' minOverlap = 15,
#' maxMismatch = 2,
#' verbose = TRUE,
#' )
#' summary <- assignTax(
#' directory_path,
#' data_tables,
#' asv_abund_matrix,
#' multithread = TRUE,
#' barcode = "rps10",
#' database_rps10 = "oomycetedb.fasta",
#' database_its = "sh_general_release_dynamic_22.08.2016.fasta"
#' )
#' obj_dada <-
#' asv_matrix_to_taxmap(
#' asv_abund_matrix,
#' min_read_depth = 10,
#' minimum_bootstrap = 75
#' )
asv_matrix_to_taxmap <- function(min_read_depth=0, minimum_bootstrap=50, pid_species=0, pid_genus=0, pid_family=0){
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  
  abund_matrix <- read_csv(file.path(directory_path,'final_asv_abundance_matrix.csv'))
  is_low_abund <- rowSums(abund_matrix[, data_tables$metadata$samplename_barcode]) < min_read_depth
  abundance <- filter(abund_matrix, ! is_low_abund)
  pid_cutoffs <- list(species = pid_species, genus = pid_genus, family = pid_family)
  
  abundance$dada2_tax <- map_chr(strsplit(abundance$dada2_tax, ';'), function(x) {
    map_chr(strsplit(x, '--'), function(parts) {
      if (parts[[3]] != "ASV" && as.numeric(parts[[2]]) <= minimum_bootstrap) {
        parts[[1]] <- "Unknown"
      }
      return(paste(parts, collapse = '--'))
    }) %>% paste(collapse = ';')
  })
  obj_dada <- metacoder::parse_tax_data(abundance, class_cols = 'dada2_tax', class_sep = ';', include_tax_data = TRUE,
                                        class_regex = '^(.+)--(.+)--(.+)$',
                                        class_key = c(taxon = 'taxon_name', boot = 'info', rank = 'taxon_rank'))
  names(obj_dada$data) <- c('abund', 'score')
  
  taxmap_path <-
    file.path(directory_path, "taxmap_obj.RData")
  save(obj_dada, file = taxmap_path)
  
  assign("obj_dada", obj_dada, envir = .GlobalEnv)
}

#' Convert taxmap object to Phyloseq object (metacoder wrapper function)
#'
#' @param taxmap_obj ASV matrix converted to taxmap object using asv_matrix_to_taxmap function
#' @return Taxmap object converted to a phyloseq object
#' @export taxmap_to_phyloseq
#' @examples
#' directory_path<-"~/Desktop/rps10test"
#' primer_path <-file.path(directory_path, "primer_info.csv")
#' metadata_path <-file.path(directory_path,"metadata.csv")
#' cutadapt_path<-"/opt/homebrew/bin/cutadapt"
#' data_tables <-
#' prepare_reads(
#' directory_path,
#' primer_path,
#' metadata_path,
#' maxN = 0,
#' )
#' cut_trim(
#' directory_path,
#' cutadapt_path,
#' verbose = TRUE,
#' maxEE = 2,
#' truncQ = 5,
#' minLen = 200,
#' maxLen = 297,
#' minCutadaptlength = 50
#')
#' asv_abund_matrix <-
#' make_asv_abund_matrix(
#' directory_path,
#' minOverlap = 15,
#' maxMismatch = 2,
#' verbose = TRUE,
#' )
#' summary <- assignTax(
#' directory_path,
#' data_tables,
#' asv_abund_matrix,
#' multithread = TRUE,
#' barcode = "rps10",
#' database_rps10 = "oomycetedb.fasta",
#' database_its = "sh_general_release_dynamic_22.08.2016.fasta"
#' )
#' obj_dada<-asv_matrix_to_taxmap(
#' asv_abund_matrix,
#' min_read_depth=10,
#' minimum_bootstrap=75
#' )
#' phylo_obj<-taxmap_to_phyloseq(obj_dada)
taxmap_to_phyloseq <- function(taxmap_obj=obj_dada) {
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  
  taxmap_obj$data$otu_table=taxmap_obj$data$abund[,-2:-4]
  taxmap_obj$data$otu_table$otu_id=paste0('ASV', 1:nrow(taxmap_obj$data$otu_table))
  taxmap_obj$data$sample_data=data_tables$metadata
  
  phylo_obj<-metacoder::as_phyloseq(
    taxmap_obj,
    sample_data = obj_dada$data$sample_data,
    sample_id_col = "samplename_barcode")
  
  assign("phylo_obj", phylo_obj, envir = .GlobalEnv)
  
  phylo_path <-
    file.path(directory_path, "phylo_obj.RData")
  save(phylo_obj, file = phylo_path)
}
