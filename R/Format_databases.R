#' General function to format the database based on barcode type
#'
#' @param analysis_setup A list containing directory paths and data tables, produced by the `prepare_reads` function.
#' @param barcode The barcode for which the database should be formatted
#' @return A formatted database based on the specified barcode type
#' @keywords internal
#'
format_database <- function(analysis_setup, barcode, db_its, db_rps10, db_16s, db_other1, db_other2) {
  if (barcode == "rps10") {
    return(format_db_rps10(analysis_setup, db_rps10))
  } else if (barcode == "its") {
    return(format_db_its(analysis_setup, db_its))
  } else if (barcode == "sixteenS") {
    return(format_db_16s(analysis_setup, db_16s))
  } else if (barcode == "other1") {
    return(format_db_other(analysis_setup, db_other1))
  } else if (barcode == "other2") {
    return(format_db_other(analysis_setup, db_other2))
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
#' @param db_rps10 The name of the database
#' @return A rps10 database that has modified headers and is output in the reference_databases folder.
#' @keywords internal
#'
format_db_rps10 <-function(analysis_setup, db_rps10){
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "rps10_reference_db.fa")
  db_rps10 <- read_fasta(file.path(data_path, db_rps10))
  data_rps10 <- str_match(names(db_rps10), pattern = "name=(.+)\\|strain=(.+)\\|ncbi_acc=(.+)\\|ncbi_taxid=(.+)\\|oodb_id=(.+)\\|taxonomy=(.+)$")
  colnames(data_rps10) <- c("header", "name", "strain", "ncbi_acc", "ncbi_taxid", "oodb_id", "taxonomy")
  data_rps10 <- as_tibble(data_rps10)
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = 'cellular_organisms;', replacement = '', fixed = TRUE)
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = 'Eukaryota', replacement = 'Eukaryota;Heterokontophyta', fixed = TRUE)
  binomial <- map_chr(str_split(data_rps10$taxonomy, pattern = ';'), `[`, 7)
  species <- map_chr(str_split(binomial, pattern = '_'), `[`, 1)
  unique(species)
  data_rps10$taxonomy <- map_chr(seq_along(data_rps10$taxonomy), function(index) {
    sub(data_rps10$taxonomy[index], pattern = binomial[index], replacement = paste0(species[index], ';', binomial[index]))
  })
  data_rps10$taxonomy <- paste0(data_rps10$taxonomy, ';', 'oodb_', seq_along(data_rps10$taxonomy))
  data_rps10$taxonomy <- paste0(data_rps10$taxonomy, ';')
  data_rps10$taxonomy <- trimws(data_rps10$taxonomy)
  db_rps10 <- trimws(db_rps10)
  #optional-need to decide on finalized database format
  stopifnot(all(str_count(data_rps10$taxonomy, pattern = ";") == 9))
  species_count <- table(map_chr(strsplit(data_rps10$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('species', 'Number of sequences')
  write_csv(count_table, file = file.path(directory_path, "species_count_table_rps10.csv"))
  write_lines(paste0(">", data_rps10$taxonomy, "\n", db_rps10), file = database_path)
  return(data_rps10)
}

#' An ITS database that has modified headers and is output in the reference_databases folder.
#'
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' metadata, and primer_info files
#' @param db_its The name of the database
#' @return An ITS database that has modified headers and is output in the reference_databases folder.
#' @keywords internal
#'
format_db_its <- function(analysis_setup, db_its) {
  # Extract directory paths from analysis_setup
  dir_paths <- analysis_setup$dir_paths
  
  # Set directory paths
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  
  # Define database path
  database_path <- file.path(directory_path_temp, "its_reference_db.fa")
  
  # Read fasta file
  db_its <- read_fasta(file.path(data_path, db_its))
  
  # Extract header information
  data_its <- str_match(names(db_its), pattern = "(.*)\\|(.*)\\|(.*)\\|(.*)\\|(.*)$")
  colnames(data_its) <- c("header", "name", "ncbi_acc", "unite_db", "db", "taxonomy")
  data_its <- as_tibble(data_its)
  
  # Identify indices of headers containing ";unidentified;"
  unidentified_indices <- grepl("__unidentified;", data_its$taxonomy)
  
  # Filter out headers containing ";unidentified;"
  filtered_data_its <- data_its[!unidentified_indices, ]
  
  # Apply modifications to the taxonomy column
  filtered_data_its$taxonomy <- gsub(filtered_data_its$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  filtered_data_its$taxonomy <- paste0('Eukaryota;', filtered_data_its$taxonomy)
  filtered_data_its$taxonomy <- gsub(filtered_data_its$taxonomy, pattern = 'Stramenopila;Oomycota', replacement = 'Heterokontophyta;Stramenopiles', fixed = TRUE)
  filtered_data_its$taxonomy <- paste0(filtered_data_its$taxonomy, ';', 'unite_', seq_along(filtered_data_its$taxonomy))
  filtered_data_its$taxonomy <- gsub(filtered_data_its$taxonomy, pattern = "[a-z]__", replacement = '')
  filtered_data_its$taxonomy <- paste0(filtered_data_its$taxonomy, ';')
  filtered_data_its$taxonomy <- trimws(filtered_data_its$taxonomy)
  
  # Extract the corresponding sequences
  filtered_db_its <- db_its[!unidentified_indices]
  
  # Extract sequences as character vector
  sequences <- as.character(filtered_db_its)
  
  # Fix after checking out later analysis
  stopifnot(all(str_count(filtered_data_its$taxonomy, pattern = ";") == 9))
  species_count <- table(map_chr(strsplit(filtered_data_its$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Species', 'Number of sequences')
  
  # Write species count table to CSV
  write_csv(count_table, file = file.path(directory_path, "species_count_table_its.csv"))
  
  # Write modified headers and corresponding sequences to database file
  write_lines(paste0(">", filtered_data_its$taxonomy, "\n", sequences), file = database_path)
  
  return(filtered_data_its)
}

#' An 16s database that has modified headers and is output in the reference_databases folder.
#'
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' metadata, and primer_info files
#' @param db_16s The name of the database
#' @return An 16s database that has modified headers and is output in the reference_databases folder.
#' @keywords internal
#'
format_db_16s <-function(analysis_setup, db_16s){
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "sixteenS_reference_db.fa")
  db_16s <- read_fasta(file.path(data_path, db_16s))
  
  data_16s <- tibble(
    taxonomy = str_replace(names(db_16s), ">", "")
  )
  colnames(data_16s) <- "taxonomy"
  data_16s <- as_tibble(data_16s)
  data_16s$taxonomy <- gsub(data_16s$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  filter_condition <- grepl("Bacteria", data_16s$taxonomy)
  # Exclude entries with Chloroplast or Mitochondria in the header
  filter_condition <- filter_condition & !grepl("Chloroplast|Mitochondria", data_16s$taxonomy)
  data_16s <- data_16s[filter_condition, ]
  
  data_16s$taxonomy <- paste0('Prokaryota;',  data_16s$taxonomy)
  
  add_NA_to_taxonomy <- function(taxonomy) {
    assignments <- strsplit(taxonomy, ";")[[1]]
    missing_assignments <- 8 - length(assignments)
    
    if (missing_assignments > 0) {
      assignments <- c(assignments, rep("NA", missing_assignments))
    }
    
    return(paste(assignments, collapse = ";"))
  }
  
  data_16s$taxonomy <- sapply(data_16s$taxonomy, add_NA_to_taxonomy)
  data_16s$taxonomy <- paste0(data_16s$taxonomy, ';', 'silva_', seq_along(data_16s$taxonomy))
  data_16s$taxonomy <- gsub(data_16s$taxonomy, pattern = "[a-z]__", replacement = '')
  data_16s$taxonomy <- paste0(data_16s$taxonomy, ';')
  data_16s$taxonomy <- trimws(data_16s$taxonomy)
  
  stopifnot(all(str_count(data_16s$taxonomy, pattern = ";") == 9))
  genus_count <- table(sapply(strsplit(data_16s$taxonomy, ";"), function(x) {
    if (length(x) >= 3) {
      return(x[length(x) - 2])
    } else {
      return(NA)
    }
  }))
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  
  write_csv(count_table, file = file.path(directory_path, "genus_count_table_16S.csv"))
  write_lines(paste0(">",  data_16s$taxonomy, "\n", db_16s), file = database_path)
  return(data_16s)
}

#' An other, user-specified database that is initially in UNITE fungal db format
#'
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' metadata, and primer_info files
#' @param db_other1 The name of the database
#' @return An other1 database that has modified headers and is output in the reference_databases folder.
#' @keywords internal
#'
format_db_other1 <-function(analysis_setup, db_other1){
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "other1_reference_db.fa")
  db_other1 <- read_fasta(file.path(data_path, db_other1))
  data_other1 <- str_match(names(db_other1), pattern = "(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)$")
  colnames(data_other1) <- c("header", "name", "ncbi_acc", "unite_db", "db", "taxonomy")
  data_other1 <- as_tibble(data_other1)
  data_other1$taxonomy <- gsub(data_other1$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_other1$taxonomy <- paste0('Eukaryota;', data_other1$taxonomy)
  data_other1$taxonomy <- gsub(data_other1$taxonomy, pattern = 'Stramenopila;Oomycota', replacement = 'Heterokontophyta;Stramenopiles', fixed = TRUE)
  data_other1$taxonomy <- paste0(data_other1$taxonomy, ';', 'unite_', seq_along(data_other1$taxonomy))
  data_other1$taxonomy <- gsub(data_other1$taxonomy, pattern = "[a-z]__", replacement = '')
  data_other1$taxonomy <- paste0(data_other1$taxonomy, ';')
  data_other1$taxonomy <- trimws(data_other1$taxonomy)
  #Fix after checking out later analysis
  stopifnot(all(str_count(data_other1$taxonomy, pattern = ";") == 9))
  species_count <- table(map_chr(strsplit(data_other1$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Species', 'Number of sequences')
  write_csv(count_table, file = file.path(directory_path, "species_count_table_other1.csv"))
  write_lines(paste0(">", data_other1$taxonomy, "\n", db_other1), file = database_path)
  return(data_other1)
}

#' An other, user-specified database that is initially in UNITE fungal db format
#'
#' @param directory_path The path to the directory containing the fastq,
#' metadata, and primer_info files
#' @param directory_path_temp User-defined temporary directory to place reads throughout the workflow
#' metadata, and primer_info files
#' @param db_other2 The name of the database
#' @return An other database that has modified headers and is output in the reference_databases folder.
#' @keywords internal
#'
format_db_other2 <-function(analysis_setup, db_other2){
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "other2_reference_db.fa")
  db_other2 <- read_fasta(file.path(data_path, db_other2))
  data_other2 <- str_match(names(db_other2), pattern = "(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)$")
  colnames(data_other2) <- c("header", "name", "ncbi_acc", "unite_db", "db", "taxonomy")
  data_other2 <- as_tibble(data_other2)
  data_other2$taxonomy <- gsub(data_other2$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_other2$taxonomy <- paste0('Eukaryota;', data_other2$taxonomy)
  data_other2$taxonomy <- gsub(data_other2$taxonomy, pattern = 'Stramenopila;Oomycota', replacement = 'Heterokontophyta;Stramenopiles', fixed = TRUE)
  data_other2$taxonomy <- paste0(data_other2$taxonomy, ';', 'unite_', seq_along(data_other2$taxonomy))
  data_other2$taxonomy <- gsub(data_other2$taxonomy, pattern = "[a-z]__", replacement = '')
  data_other2$taxonomy <- paste0(data_other2$taxonomy, ';')
  data_other2$taxonomy <- trimws(data_other2$taxonomy)
  #Fix after checking out later analysis
  stopifnot(all(str_count(data_other2$taxonomy, pattern = ";") == 9))
  species_count <- table(map_chr(strsplit(data_other2$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Species', 'Number of sequences')
  write_csv(count_table, file = file.path(directory_path, "species_count_table_other2.csv"))
  write_lines(paste0(">", data_other2$taxonomy, "\n", db_other2), file = database_path)
  return(data_other2)
}