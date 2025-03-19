utils::globalVariables(c("data_other1", "db_other1", "data_other2", "db_other2"))

#' Create modified reference rps10 database for downstream analysis
#' @param data_tables The data tables containing the paths to read files, metadata, primer sequences
#' @param data_path Path to the data directory
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param temp_directory_path User-defined temporary directory to place reads
#'   throughout the workflow metadata, and primer_info files
#' @param db_rps10 The name of the database
#' @return A rps10 database that has modified headers and is output in the
#'   reference_databases folder.
#' @keywords internal

format_db_rps10 <- function(data_tables, data_path, output_directory_path, temp_directory_path, db_rps10) {
  
  database_path <- file.path(temp_directory_path, "rps10_reference_db.fa")
  
  db_rps10 <- metacoder::read_fasta(file.path(data_path, db_rps10))
  
  data_rps10 <- stringr::str_match(names(db_rps10), pattern = "name=(.+)\\|strain=(.+)\\|ncbi_acc=(.+)\\|ncbi_taxid=(.+)\\|oodb_id=(.+)\\|taxonomy=(.+)$")
  colnames(data_rps10) <- c("header", "name", "strain", "ncbi_acc", "ncbi_taxid", "oodb_id", "taxonomy")
  data_rps10 <- tibble::as_tibble(data_rps10)
  
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = 'cellular_organisms;', replacement = '', fixed = TRUE)
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = 'Eukaryota;Stramenopiles', replacement = 'Stramenopila;Oomycota', fixed = TRUE)
  
  binomial <- purrr::map_chr(stringr::str_split(data_rps10$taxonomy, pattern = ';'), `[`, 6)
  species <- purrr::map_chr(stringr::str_split(binomial, pattern = '_'), `[`, 1)
  
  data_rps10$taxonomy <- purrr::map_chr(seq_along(data_rps10$taxonomy), function(index) {
    sub(data_rps10$taxonomy[index], pattern = binomial[index], replacement = paste0(species[index], ';', binomial[index]))
  })
  
  data_rps10$taxonomy <- trimws(data_rps10$taxonomy)
  data_rps10$taxonomy <- paste0(data_rps10$taxonomy, ';')
  
  readr::write_lines(paste0(">", data_rps10$taxonomy, "\n", db_rps10), file = database_path)
  
  species_count <- table(binomial)
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- tibble::tibble(count_table)
  names(count_table) <- c('Species', 'Number of sequences')
  
  readr::write_csv(count_table, file = file.path(output_directory_path, "species_count_table_rps10.csv"))
  readr::write_lines(paste0(">", data_rps10$taxonomy, "\n", db_rps10), file = database_path)
  
  return(data_rps10)
}
#' An ITS database that has modified headers and is output in the
#' reference_databases folder
#' @param data_tables The data tables containing the paths to read files, metadata, primer sequences
#' @param data_path Path to the data directory
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param temp_directory_path User-defined temporary directory to place reads
#'   throughout the workflow metadata, and primer_info files
#' @param db_its The name of the database
#' @return An ITS database that has modified headers and is output in the
#'   reference_databases folder.
#' @keywords internal
format_db_its <- function(data_tables, data_path, output_directory_path, temp_directory_path, db_its) {
  database_path <- file.path(temp_directory_path, "its_reference_db.fa")
  
  db_its <- metacoder::read_fasta(file.path(data_path, db_its))
  data_its <- tibble::tibble(header = names(db_its), sequence = db_its)
  
  data_its$taxonomy <- stringr::str_extract(data_its$header, "(?<=\\|)[^|]+$")
  data_its$taxonomy <- gsub("[A-Za-z]__", "", data_its$taxonomy)
  data_its$taxonomy <- gsub(" ", "_", data_its$taxonomy)
  
  data_its$taxonomy <- sapply(data_its$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  data_its$taxonomy <- trimws(data_its$taxonomy)
  data_its$taxonomy <- paste0(data_its$taxonomy, ';')
  

  data_its$genus <- ifelse(
    sapply(strsplit(data_its$taxonomy, ";"), length) >= 7,
    sapply(strsplit(data_its$taxonomy, ";"), function(x) x[6]),
    NA
  )
  
  genus_count <- table(data_its$genus)
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- tibble::tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  
  
  readr::write_csv(count_table, file = file.path(output_directory_path, "genus_count_table_its.csv"))
  readr::write_lines(paste0(">", data_its$taxonomy, "\n", data_its$sequence), file = database_path)
  return(data_its)
}

#' An 16S database that has modified headers and is output in the
#' reference_databases folder
#' @param data_tables The data tables containing the paths to read files, metadata, primer sequences
#' @param data_path Path to the data directory
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param temp_directory_path User-defined temporary directory to place reads
#'   throughout the workflow metadata, and primer_info files
#' @param db_16S The name of the database
#'
#' @return An 16S database that has modified headers and is output in the
#'   reference_databases folder
#'
#' @keywords internal
format_db_16S <- function(data_tables, data_path, output_directory_path, temp_directory_path, db_16S) {

  database_path <- file.path(temp_directory_path, "r16S_reference_db.fa")
  
  db_16S <- metacoder::read_fasta(file.path(data_path, db_16S))
  data_16S <- tibble::tibble(taxonomy = names(db_16S), sequence = db_16S)
  data_16S$taxonomy <- gsub(data_16S$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_16S$taxonomy <- trimws(data_16S$taxonomy)
  data_16S$taxonomy <- stringr::str_replace(data_16S$taxonomy, "^((?:[^;]*;){5})([^;]+);([^;]+);$", "\\1\\2;\\2_\\3;")
  
  data_16S$taxonomy <- sapply(data_16S$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  data_16S$taxonomy <- trimws(data_16S$taxonomy)
  data_16S$taxonomy <- paste0(data_16S$taxonomy, ';')
  
  data_16S$genus <- ifelse(
    sapply(strsplit(data_16S$taxonomy, ";"), length) >= 7,
    sapply(strsplit(data_16S$taxonomy, ";"), function(x) x[6]),
    NA
  )
  
  genus_count <- table(data_16S$genus)
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- tibble::tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  
  readr::write_csv(count_table, file = file.path(output_directory_path, "genus_count_table_16S.csv"))
  readr::write_lines(paste0(">", data_16S$taxonomy, "\n", db_16S), file = database_path)
  
  return(data_16S)
}

#' An other, user-specified database that is initially in the format specified by DADA2 with header simply taxonomic levels (kingdom down to species, separated by semi-colons, ;)
#' @param data_tables The data tables containing the paths to read files, metadata, primer sequences
#' @param data_path Path to the data directory
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param temp_directory_path User-defined temporary directory to place reads
#'   throughout the workflow metadata, and primer_info files
#' @param db_other1 The name of the database
#'
#' @return A database (other than SILVA, UNITE, and oomyceteDB that has modified headers and is output in the
#'   reference_databases folder
#' 
#' @return An other database that has modified headers and is output in the reference_databases folder.
#' 
#' @keywords internal
format_db_other1 <-function(data_tables, data_path, output_directory_path, temp_directory_path, db_other1){
  
  database_path <- file.path(temp_directory_path, "other1_reference_db.fa")
  
  db_other1 <- metacoder::read_fasta(file.path(data_path, db_other1))
  data_other1 <- tibble::tibble(taxonomy = names(db_other1), sequence = db_other1)
  data_other1$taxonomy <- gsub(data_other1$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_other1$taxonomy <- trimws(data_other1$taxonomy)
  
  data_other1$taxonomy <- sapply(data_other1$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  
  data_other1$taxonomy <- paste0(data_other1$taxonomy, ";")

  data_other1$genus <- ifelse(
    sapply(strsplit(data_other1$taxonomy, ";"), length) == 7,
    sapply(strsplit(data_other1$taxonomy, ";"), function(x) x[6]),
    ifelse(
      sapply(strsplit(data_other1$taxonomy, ";"), length) == 6,
      sapply(strsplit(data_other1$taxonomy, ";"), function(x) x[6]),
      NA
    )
  )
  genus_count <- table(data_other1$genus)
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- tibble::tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')

  readr::write_csv(count_table, file = file.path(output_directory_path, "genus_count_table_other1.csv"))
  readr::write_lines(paste0(">", data_other1$taxonomy, "\n", db_other1), file = database_path)
  return(data_other1)
}

#' An second user-specified database that is initially in the format specified
#' by DADA2 with header simply taxonomic levels (kingdom down to species,
#' separated by semi-colons, ;)
#' 
#' @param data_tables The data tables containing the paths to read files, metadata, primer sequences
#' @param data_path Path to the data directory
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param temp_directory_path User-defined temporary directory to place reads
#'   throughout the workflow metadata, and primer_info files
#' @param db_other2 The name of the database
#'
#' @return An second database (other than SILVA, UNITE, and oomyceteDB that has modified headers and is output in the
#'   reference_databases folder
#'
#' @keywords internal
format_db_other2 <-function(data_tables, data_path, output_directory_path, temp_directory_path, db_other2){
  
  database_path <- file.path(temp_directory_path, "other2_reference_db.fa")
  
  db_other2 <- metacoder::read_fasta(file.path(data_path, db_other2))
  data_other2 <- tibble::tibble(taxonomy = names(db_other2), sequence = db_other2)
  data_other2$taxonomy <- gsub(data_other2$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_other2$taxonomy <- trimws(data_other2$taxonomy)
  
  data_other2$taxonomy <- sapply(data_other2$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  #data_other2$taxonomy <- paste0(data_other2$taxonomy, ";refdb_", seq_along(data_other2$taxonomy), ";")
  data_other2$taxonomy <- paste0(data_other2$taxonomy, ";")
  
  data_other2$genus <- ifelse(
    sapply(strsplit(data_other2$taxonomy, ";"), length) == 7,
    sapply(strsplit(data_other2$taxonomy, ";"), function(x) x[6]),
    ifelse(
      sapply(strsplit(data_other2$taxonomy, ";"), length) == 6,
      sapply(strsplit(data_other2$taxonomy, ";"), function(x) x[6]),
      NA
    )
  )
  genus_count <- table(data_other2$genus)
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- tibble::tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  
  readr::write_csv(count_table, file = file.path(output_directory_path, "genus_count_table_other2.csv"))
  readr::write_lines(paste0(">", data_other1$taxonomy, "\n", db_other1), file = database_path)
  return(data_other2)
}

#' General functions to format user-specified databases
#' @importFrom utils modifyList read.table stack
#' @param data_tables The data tables containing the paths to read files, metadata, primer sequences
#' @param data_path Path to the data directory
#' @param output_directory_path The path to the directory where resulting files
#'   are output
#' @param temp_directory_path User-defined temporary directory to place reads
#'   throughout the workflow metadata, and primer_info files
#' @param barcode The barcode for which the database should be formatted
#'
#' @return A formatted database based on the specified barcode type
#'
#' @keywords internal
format_database <- function(data_tables, data_path, output_directory_path, temp_directory_path, barcode, db_its, db_rps10, db_16S, db_other1, db_other2) {
  if (barcode == "rps10") {
    return(format_db_rps10(data_tables, data_path, output_directory_path, temp_directory_path, db_rps10))
  } else if (barcode == "its") {
    return(format_db_its(data_tables, data_path, output_directory_path, temp_directory_path, db_its))
  } else if (barcode == "r16S") {
    return(format_db_16S(data_tables, data_path, output_directory_path, temp_directory_path, db_16S))
  } else if (barcode == "other1") {
    return(format_db_other1(data_tables, data_path, output_directory_path, temp_directory_path, db_other1))
  } else if (barcode == "other2") {
    return(format_db_other2(data_tables, data_path, output_directory_path, temp_directory_path, db_other2))
  } else {
    stop("Barcode not recognized: ", barcode)
  }
}