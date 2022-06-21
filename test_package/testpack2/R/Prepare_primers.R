#' Take in user primers and creates the complement, reverse, reverse complement of primers in one tibble
#'
#' @param primer_path a path to the csv file that holds the primer info for this project
#'
#' @return A tibble that contains the forward and reverse primers with all complements of primers
#' @export
#'
#' @examples
prepare_primers <- function(primer_path){
  primer_data_path <- file.path(primer_path)
  primer_data <- read_csv(primer_data_path)

  #seperate forward and reverse to make various primers
  forward_primers <- primer_data[, c(1:2)]
  toString(forward_primers)
  reverse_primers <- primer_data[, c(1,3)]

  forward_primers$f_compt <- map_chr(forward_primers$forward, function(x) toString(Biostrings::complement(DNAString(x))))
  forward_primers$f_rev <- map_chr(forward_primers$forward, function(x) toString(Biostrings::reverse(DNAString(x))))
  forward_primers$f_rc <- map_chr(forward_primers$forward, function(x) toString(Biostrings::reverseComplement(DNAString(x))))

  reverse_primers$r_compt <- map_chr(reverse_primers$reverse, function(x) toString(Biostrings::complement(DNAString(x))))
  reverse_primers$r_rev <- map_chr(reverse_primers$reverse, function(x) toString(Biostrings::reverse(DNAString(x))))
  reverse_primers$r_rc <- map_chr(reverse_primers$reverse, function(x) toString(Biostrings::reverseComplement(DNAString(x))))

  #add back together
  primer_data <- forward_primers %>%
    left_join(reverse_primers, by = "primer_name")

  return(primer_data)

}

