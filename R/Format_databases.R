#' General function to format the database based on barcode type
#'
#' @param analysis_setup A list containing directory paths and data tables, produced by the `prepare_reads` function.
#' @param barcode The barcode for which the database should be formatted
#' @return A formatted database based on the specified barcode type
#' @keywords internal
#'
format_database <- function(analysis_setup, barcode, db_its, db_rps10, db_16S, db_other) {
  if (barcode == "rps10") {
    return(format_database_rps10(analysis_setup, db_rps10))
  } else if (barcode == "its") {
    return(format_database_its(analysis_setup, db_its))
  } else if (barcode == "16s") {
    return(format_database_16s(analysis_setup, db_its))
  } else if (barcode == "other") {
    return(format_database_other(analysis_setup, db_its))
  } else {
    stop("Barcode not recognized")
  }
}

#' Create modified reference rps10 database for downstream analysis
#'
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' metadata, and primer_info files
#' @param database_rps10 The name of the database
#' @return A rps10 database that has modified headers and is output in the reference_databases folder.
#' @keywords internal
#'
format_database_rps10 <-function(analysis_setup, database_rps10){
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "rps10_reference_db.fa")
  db_rps10 <- read_fasta(file.path(data_path, database_rps10))
  rps10_data <- str_match(names(db_rps10), pattern = "name=(.+)\\|strain=(.+)\\|ncbi_acc=(.+)\\|ncbi_taxid=(.+)\\|oodb_id=(.+)\\|taxonomy=(.+)$")
  colnames(rps10_data) <- c("header", "name", "strain", "ncbi_acc", "ncbi_taxid", "oodb_id", "taxonomy")
  rps10_data <- as_tibble(rps10_data)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = 'cellular_organisms;', replacement = '', fixed = TRUE)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = 'Eukaryota', replacement = 'Eukaryota;Heterokontophyta', fixed = TRUE)
  binomial <- map_chr(str_split(rps10_data$taxonomy, pattern = ';'), `[`, 7)
  species <- map_chr(str_split(binomial, pattern = '_'), `[`, 1)
  unique(species)
  rps10_data$taxonomy <- map_chr(seq_along(rps10_data$taxonomy), function(index) {
    sub(rps10_data$taxonomy[index], pattern = binomial[index], replacement = paste0(species[index], ';', binomial[index]))
  })
  rps10_data$taxonomy <- paste0(rps10_data$taxonomy, ';', 'oodb_', seq_along(rps10_data$taxonomy))
  rps10_data$taxonomy <- paste0(rps10_data$taxonomy, ';')
  rps10_data$taxonomy <- trimws(rps10_data$taxonomy)
  db_rps10 <- trimws(db_rps10)
  #optional-need to decide on finalized database format
  stopifnot(all(str_count(rps10_data$taxonomy, pattern = ";") == 9))
  species_count <- table(map_chr(strsplit(rps10_data$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('species', 'Number of sequences')
  write_csv(count_table, file = file.path(directory_path, "rps10_species_count_table.csv"))
  write_lines(paste0(">", rps10_data$taxonomy, "\n", db_rps10), file = database_path)
  return(rps10_data)
}


#' An ITS database that has modified headers and is output in the reference_databases folder.
#'
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' metadata, and primer_info files
#' @param database_its The name of the database
#' @return An ITS database that has modified headers and is output in the reference_databases folder.
#' @keywords internal
#'
format_database_its <-function(analysis_setup, database_its){
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "its_reference_db.fa")
  db_its <- read_fasta(file.path(data_path, database_its))
  its_data <- str_match(names(db_its), pattern = "(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)$")
  colnames(its_data) <- c("header", "name", "ncbi_acc", "unite_db", "db", "taxonomy")
  its_data <- as_tibble(its_data)
  its_data$taxonomy <- gsub(its_data$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  its_data$taxonomy <- paste0('Eukaryota;', its_data$taxonomy)
  its_data$taxonomy <- gsub(its_data$taxonomy, pattern = 'Stramenopila;Oomycota', replacement = 'Heterokontophyta;Stramenopiles', fixed = TRUE)
  its_data$taxonomy <- paste0(its_data$taxonomy, ';', 'unite_', seq_along(its_data$taxonomy))
  its_data$taxonomy <- gsub(its_data$taxonomy, pattern = "[a-z]__", replacement = '')
  its_data$taxonomy <- paste0(its_data$taxonomy, ';')
  its_data$taxonomy <- trimws(its_data$taxonomy)
  #Fix after checking out later analysis
  stopifnot(all(str_count(its_data$taxonomy, pattern = ";") == 9))
  species_count <- table(map_chr(strsplit(its_data$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Species', 'Number of sequences')
  write_csv(count_table, file = file.path(directory_path, "its_species_count_table.csv"))
  write_lines(paste0(">", its_data$taxonomy, "\n", db_its), file = database_path)
  return(its_data)
}
