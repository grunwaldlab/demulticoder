#################################################################################################################
## Database functions-Feb 2022
## Purpose: Make usable databases that people can reference without having to do significant reformatting during the analysis phase 
## Inputs: Raw databases from http://oomycetedb.cgrb.oregonstate.edu/search.html and https://unite.ut.ee/repository.php
## Returns: Reformatted database and a .csv of genus counts for each database. 
## Misc notes: Still in progress. Adapted from Zach Foster's pipeline scripts. Each database has slightly different header formats, and it is a challenge to generalize these functions. 

##example headers
##Rps10 >name=Aphanomyces_invadans|strain=NJM9701|ncbi_acc=KX405005|ncbi_taxid=157072|oodb_id=13|taxonomy=cellular_organisms;Eukaryota;Stramenopiles;Oomycetes;Saprolegniales;Saprolegniaceae;Aphanomyces_invadans
##ITS Unite header >Cryptosphaeria_sp|KU761147|SH1572218.08FU|reps|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Diaporthales;f__Valsaceae;g__Cryptosphaeria;s__Cryptosphaeria_sp

################################################################################################################## 
#load packages 
library(dada2)
library(ShortRead)
library(Biostrings)
library(dplyr)
library(purrr)
library(furrr)
library(tidyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(metacoder)
library(readr)
library(tibble)
library(stringr)
library(taxize)
library(DT)
library(ape)
library(sessioninfo)


create_database <- function(intermediate_path){
  formatted_ref_dir <- file.path(intermediate_path, "reference_databases")
  if (! dir.exists(formatted_ref_dir)) {
    dir.create(formatted_ref_dir)
  }
}

create_database("~/rhod_metabarcoding4/intermediate_data")
#Download and place your desired database into the folder

#include a function to download an database? 

#Format rps10_database first attempt
format_database <-function(database_path, database_name){
  rps10_db <- read_fasta(file.path(database_path, database_name))
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
  write_csv(count_table, file = file.path(database_path, "rps10_genus_count_table.csv"))
  rps10_ref_path <- file.path(database_path, "rps10_reference_db.fa")
  paste0(">", rps10_data$taxonomy, "\n", rps10_db) %>%
    write_lines(file = rps10_ref_path)
  return(rps10_data)
}

#test commands
rps10_data<-format_database("~/rhod_metabarcoding4/databases", "oomycetedb.fa")

## ITS database
format_database2 <-function(database_path, database_name){
  its_db <- read_fasta(file.path(database_path, database_name))
  its_data <- str_match(names(its_db), pattern = "(.+)\\|(.+)\\|(.+)\\|(.+)\\|(.+)$")
  colnames(its_data) <- c("header", "name", "ncbi_acc", "unite_db", "db", "taxonomy")
  its_data <- as_tibble(its_data)
  its_data$taxonomy <- gsub(its_data$taxonomy, pattern = ' ', replacement = '_', fixed = TRUE)
  its_data$taxonomy <- paste0(its_data$taxonomy, ';')
  its_data$taxonomy <- trimws(its_data$taxonomy)
  #Fix after checking out later analysis
  stopifnot(all(str_count(its_data$taxonomy, pattern = ";") == 7))
  genus_count <- table(map_chr(strsplit(its_data$name, split = '_'), `[`, 1))
  count_table <- as.data.frame(genus_count, stringsAsFactors = FALSE)
  count_table <- as_tibble(count_table)
  names(count_table) <- c('Genus', 'Number of sequences')
  write_csv(count_table, file = file.path(database_path, "its_genus_count_table.csv"))
  its_ref_path <- file.path(database_path, "its_reference_db.fa")
  paste0(">", its_data$taxonomy, "\n", its_db) %>%
    write_lines(file = its_ref_path)
  return(its_data)
  
}

#test commands
its_data<-format_database2("~/rhod_metabarcoding4/intermediate_data/reference_databases", "unite.fasta")