#' Create modified reference rps10 database for downstream analysis
#'
#' @param directory_path A path to a directory that contains raw data
#' @param database_rps10 The name of the database
#'
#' @return A rps10 database that has modified headers and is output in the reference_databases folder.
#' @export
#'
#' @examples
format_database_rps10 <-function(directory_path, database_rps10){
  database_path <- file.path(directory_path, "rps10_reference_db.fa")
  rps10_db <- read_fasta(file.path(directory_path, database_rps10))
  rps10_data <- str_match(names(rps10_db), pattern = "name=(.+)\\|strain=(.+)\\|ncbi_acc=(.+)\\|ncbi_taxid=(.+)\\|oodb_id=(.+)\\|taxonomy=(.+)$")
  colnames(rps10_data) <- c("header", "name", "strain", "ncbi_acc", "ncbi_taxid", "oodb_id", "taxonomy")
  rps10_data <- as_tibble(rps10_data)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = 'cellular_organisms;', replacement = '', fixed = TRUE)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  rps10_data$taxonomy <- gsub(rps10_data$taxonomy, pattern = 'Eukaryota', replacement = 'Eukaryota;Heterokontophyta', fixed = TRUE)
  binomial <- map_chr(str_split(rps10_data$taxonomy, pattern = ';'), `[`, 7)
  genus <- map_chr(str_split(binomial, pattern = '_'), `[`, 1)
  unique(genus)
  rps10_data$taxonomy <- map_chr(seq_along(rps10_data$taxonomy), function(index) {
    sub(rps10_data$taxonomy[index], pattern = binomial[index], replacement = paste0(genus[index], ';', binomial[index]))
  })
  rps10_data$taxonomy <- paste0(rps10_data$taxonomy, ';', 'oodb_', seq_along(rps10_data$taxonomy))
  rps10_data$taxonomy <- paste0(rps10_data$taxonomy, ';')
  rps10_data$taxonomy <- trimws(rps10_data$taxonomy)
  rps10_db <- trimws(rps10_db)
  #optional-need to decide on finalized database format
  stopifnot(all(str_count(rps10_data$taxonomy, pattern = ";") == 9))
  genus_count <- table(map_chr(strsplit(rps10_data$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  #directory_path <- file.path(directory_path, "reference_databases")
  write_csv(count_table, file = file.path(directory_path, "rps10_genus_count_table.csv"))
  rps10_ref_path <- file.path(directory_path, "rps10_reference_db.fa")
  paste0(">", rps10_data$taxonomy, "\n", rps10_db) %>%
    write_lines(file = rps10_ref_path)
  return(rps10_data)
}


#' An ITS database that has modified headers and is output in the reference_databases folder.
#'
#' @param directory_path A path to a directory that contains raw data
#' @param database_its The name of the database
#'
#' @return An ITS database that has modified headers and is output in the reference_databases folder.
#' @export
#'
#' @examples
format_database_its <-function(directory_path, database_its){
  database_path <- file.path(directory_path, "its_reference_db.fa")
  its_db <- read_fasta(file.path(directory_path, database_its))
  its_data <- str_match(names(its_db), pattern = "(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)$")
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
  genus_count <- table(map_chr(strsplit(its_data$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  write_csv(count_table, file = file.path(directory_path, "its_genus_count_table.csv"))
  its_ref_path <- file.path(directory_path , "its_reference_db.fa")
  paste0(">", its_data$taxonomy, "\n", its_db) %>%
    write_lines(file = its_ref_path)
  return(its_data)
}

#Is there a way to generalize this more?

#' An 16S database that has modified headers and is output in the reference_databases folder.
#'
#' @param directory_path A path to a directory that contains raw data
#' @param database_its The name of the database
#'
#' @return A rps10 database that has modified headers and is output in the reference_databases folder.
#' @export
#'
#' @examples
format_database_16s <-function(directory_path, database_16s){
  database_path <- file.path(directory_path, "its_reference_db.fa")
  its_db <- read_fasta(file.path(directory_path, database_its))
  its_data <- str_match(names(its_db), pattern = "(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)$")
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
  genus_count <- table(map_chr(strsplit(its_data$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  write_csv(count_table, file = file.path(directory_path, "its_genus_count_table.csv"))
  its_ref_path <- file.path(directory_path , "its_reference_db.fa")
  paste0(">", its_data$taxonomy, "\n", its_db) %>%
    write_lines(file = its_ref_path)
  return(its_data)
}
