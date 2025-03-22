#' Filter ASV abundance matrix and convert to taxmap and phyloseq objects
#' @importFrom utils modifyList read.table stack
#' @param analysis_setup An object containing directory paths and
#'   data tables, produced by the `prepare_reads` function
#' @param min_read_depth ASV filter parameter. If mean read depth of across all
#'   samples is less than this threshold, ASV will be filtered.
#' @param minimum_bootstrap Set threshold for bootstrap support value for taxonomic
#'   assignments. Below designated minimum bootstrap threshold, taxnomoic
#'   assignments will be set to N/A
#'@param save_outputs Logical, indicating whether to save the resulting phyloseq
#'  and taxmap objects. If TRUE, the objects will be saved; if FALSE, they will
#'  only be available in the global environment. Default is FALSE.
#'
#' @return ASV matrix converted to taxmap object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert final matrix to taxmap and phyloseq objects for downstream analysis steps
#' analysis_setup <- prepare_reads(
#'   data_directory = system.file("extdata", package = "demulticoder"),
#'   output_directory = tempdir(),
#'   tempdir_path = tempdir(),
#'   tempdir_id = "demulticoder_run_temp",
#'   overwrite_existing = TRUE
#' )
#' cut_trim(
#' analysis_setup,
#' cutadapt_path="/usr/bin/cutadapt",
#' overwrite_existing = TRUE
#' )
#' make_asv_abund_matrix(
#' analysis_setup,
#' overwrite_existing = TRUE
#' )
#' assign_tax(
#' analysis_setup,
#' asv_abund_matrix,
#' retrieve_files=FALSE,
#' overwrite_existing=TRUE
#' )
#' objs<-convert_asv_matrix_to_objs(
#' analysis_setup
#' )
#' }

convert_asv_matrix_to_objs <-
  function(analysis_setup,
           min_read_depth = 0,
           minimum_bootstrap = 0,
           save_outputs = FALSE) {
    data_tables <- analysis_setup$data_tables
    output_directory_path <-
      analysis_setup$directory_paths$output_directory
    
    files <-
      list.files(path = output_directory_path,
                 pattern = "^final_asv_abundance_matrix_.*\\.csv$",
                 full.names = TRUE)
    suffixes <-
      gsub("^final_asv_abundance_matrix_(.*)\\.csv$",
           "\\1",
           basename(files))
    unique_suffixes <- unique(suffixes)
    
    result_list <- list()
    
    for (suffix in unique_suffixes) {
      taxmap_name <- paste0("taxmap_obj_", suffix)
      phyloseq_name <- paste0("phylo_obj_", suffix)
      
      taxmap_path <-
        file.path(output_directory_path, paste0(taxmap_name, ".RData"))
      phyloseq_path <-
        file.path(output_directory_path, paste0(phyloseq_name, ".RData"))
      
      if (suffix == "rps10") {
        abundance <-
          readr::read_csv(file.path(
            output_directory_path,
            paste0('final_asv_abundance_matrix_', suffix, '.csv')
          ))
        abundance$dada2_tax <-
          gsub(";oomycetedb_[0-9]+--[0-9]+--Reference",
               "",
               abundance$dada2_tax)
      } else {
        abundance <-
          readr::read_csv(file.path(
            output_directory_path,
            paste0('final_asv_abundance_matrix_', suffix, '.csv')
          ))
      }
      
      is_low_abund <-
        rowSums(abundance[, grepl(paste0("_", suffix, "$"), colnames(abundance))]) < min_read_depth
      abundance <- dplyr::filter(abundance,!is_low_abund)
      abundance$dada2_tax <-
        purrr::map_chr(strsplit(abundance$dada2_tax, ';'), function(x) {
          paste(sapply(strsplit(x, '--'),
                       function(parts) {
                         if (parts[[3]] != "ASV" &&
                             as.numeric(parts[[2]]) <= minimum_bootstrap) {
                           parts[[1]] <- "Unsupported"
                         }
                         paste(parts, collapse = '--')
                       }),
                collapse = ';')
        })
      taxmap_obj <-
        metacoder::parse_tax_data(
          abundance,
          class_cols = 'dada2_tax',
          class_sep = ';',
          include_tax_data = TRUE,
          class_regex = '^(.+)--(.+)--(.+)$',
          class_key = c(
            taxon = 'taxon_name',
            boot = 'info',
            rank = 'taxon_rank'
          )
        )
      names(taxmap_obj$data) <- c('abund', 'score')
      
      taxmap_obj$data$otu_table = taxmap_obj$data$abund[,-3:-4]
      filtered_sample_data <-
        data_tables$metadata[data_tables$metadata$primer_name == suffix,]
      taxmap_obj$data$sample_data = filtered_sample_data
      
      
      phylo_obj <- metacoder::as_phyloseq(
        taxmap_obj,
        sample_data = taxmap_obj$data$sample_data,
        sample_id_col = "samplename_barcode",
        otu_id_col = "asv_id"
      )
      
      result_list[[paste0("taxmap_", suffix)]] <- taxmap_obj
      result_list[[paste0("phyloseq_", suffix)]] <- phylo_obj
      
      
      save(taxmap_obj, file = taxmap_path)
      save(phylo_obj, file = phyloseq_path)
      
      cat("For", suffix, "dataset", "\n")
      cat("Taxmap object saved in:", taxmap_path, "\n")
      cat("Phyloseq object saved in:", phyloseq_path, "\n")
      cat("ASVs filtered by minimum read depth:", min_read_depth, "\n")
      cat(
        "For taxonomic assignments, if minimum bootstrap was set to:",
        minimum_bootstrap,
        "assignments were set to 'Unsupported'",
        "\n"
      )
      cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    }
    
    return(result_list)
  }