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
    stop("Barcode not recognized: ", barcode)
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
format_db_rps10 <- function(analysis_setup, db_rps10) {
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "rps10_reference_db.fa")

  db_rps10 <- read_fasta(file.path(data_path, db_rps10))
  data_rps10 <- tibble(header = names(db_rps10), sequence = db_rps10)
  data_rps10$taxonomy <- gsub(".*taxonomy=([^|]+).*", "\\1", data_rps10$header)
  
  data_rps10$name <- gsub(".*name=([^|]+)\\|.*", "\\1", data_rps10$header)
  data_rps10$name <- gsub(data_rps10$name, pattern = " ", replacement = "_")
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = 'cellular_organisms;', replacement = '', fixed = TRUE)
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_rps10$taxonomy <- gsub(data_rps10$taxonomy, pattern = 'Eukaryota', replacement = 'Eukaryota;Heterokontophyta', fixed = TRUE)
  data_rps10$taxonomy <- trimws(data_rps10$taxonomy)
  
  data_rps10$taxonomy <- sapply(data_rps10$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  data_rps10$taxonomy <- paste0(data_rps10$taxonomy, ";refdb_", seq_along(data_rps10$taxonomy), ";")
  
  data_rps10$name <-gsub(data_rps10$name, pattern = " ", replacement = "_")
  species_count <- table(data_rps10$name)
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Species', 'Number of sequences')
  
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
  dir_paths <- analysis_setup$dir_paths
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "its_reference_db.fa")
  
  db_its <- read_fasta(file.path(data_path, db_its))
  data_its <- tibble(header = names(db_its), sequence = db_its)
  
  data_its$taxonomy <- str_extract(data_its$header, "(?<=\\|)[^|]+$")
  
  data_its$taxonomy <- gsub(" ", "_", data_its$taxonomy)
  data_its$taxonomy <- gsub(data_its$taxonomy, pattern = 'Stramenopila;Oomycota', replacement = 'Heterokontophyta;Stramenopiles', fixed = TRUE)
  data_its$taxonomy <- gsub(data_its$taxonomy, pattern = "[a-z]__", replacement = '')
  data_its$taxonomy <- paste0(data_its$taxonomy, ";")
  data_its$taxonomy <- trimws(data_its$taxonomy)
  
  data_its$taxonomy <- sapply(data_its$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  data_its$taxonomy <- paste0(data_its$taxonomy, ";refdb_", seq_along(data_its$taxonomy), ";")
  
  data_its$name <- sub("^([^|]+)\\|.*$", "\\1", data_its$header)
  species_count <- table(data_its$name)
  count_table <- as.data.frame(species_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Species', 'Number of sequences')
  
  write_csv(count_table, file = file.path(directory_path, "species_count_table_its.csv"))
  write_lines(paste0(">", data_its$taxonomy, "\n", data_its$sequence), file = database_path)
  return(data_its)
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
format_db_16s <- function(analysis_setup, db_16s) {
  dir_paths <- analysis_setup$dir_paths
  data_tables <- analysis_setup$data_tables
  directory_path <- dir_paths$output_directory
  data_path <- dir_paths$data_directory
  directory_path_temp <- dir_paths$temp_directory
  database_path <- file.path(directory_path_temp, "sixteenS_reference_db.fa")
  
  db_16s <- read_fasta(file.path(data_path, db_16s))
  data_16s <- tibble(taxonomy = names(db_16s), sequence = db_16s)
  data_16s$taxonomy <- gsub(data_16s$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_16s$taxonomy <- trimws(data_16s$taxonomy)
  data_16s$taxonomy <- str_replace(data_16s$taxonomy, "^((?:[^;]*;){5})([^;]+);([^;]+);$", "\\1\\2;\\2_\\3;")
  
  data_16s$taxonomy <- sapply(data_16s$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  data_16s$taxonomy <- paste0(data_16s$taxonomy, ";refdb_", seq_along(data_16s$taxonomy), ";")
  
  data_16s$genus <- ifelse(
    sapply(strsplit(data_16s$taxonomy, ";"), length) >= 7,
    sapply(strsplit(data_16s$taxonomy, ";"), function(x) x[6]),
    NA
  )
  
  genus_count <- table(data_16s$genus)
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  
  write_csv(count_table, file = file.path(directory_path, "genus_count_table_16s.csv"))
  write_lines(paste0(">", data_16s$taxonomy, "\n", db_16s), file = database_path)
  
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
  data_other1 <- tibble(taxonomy = names(db_other1), sequence = db_other1)
  data_other1$taxonomy <- gsub(data_other1$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_other1$taxonomy <- trimws(data_other1$taxonomy)
  
  data_other1$taxonomy <- sapply(data_other1$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  data_other1$taxonomy <- paste0(data_other1$taxonomy, ";refdb_", seq_along(data_other1$taxonomy), ";")
  
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
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')

  write_csv(count_table, file = file.path(directory_path, "genus_count_table_other1.csv"))
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
  data_other2 <- tibble(taxonomy = names(db_other2), sequence = db_other2)
  data_other2$taxonomy <- gsub(data_other2$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  data_other2$taxonomy <- trimws(data_other2$taxonomy)
  
  data_other2$taxonomy <- sapply(data_other2$taxonomy, function(tax) {
    tax_parts <- strsplit(tax, ";")[[1]]
    while (length(tax_parts) < 7) {
      tax_parts <- c(tax_parts, "NA")
    }
    paste(tax_parts, collapse = ";")
  })
  
  data_other2$taxonomy <- paste0(data_other2$taxonomy, ";refdb_", seq_along(data_other2$taxonomy), ";")
  
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
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  
  write_csv(count_table, file = file.path(directory_path, "genus_count_table_other2.csv"))
  write_lines(paste0(">", data_other1$taxonomy, "\n", db_other1), file = database_path)
  return(data_other2)
}
