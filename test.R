devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
library("demulticoder")
devtools::clean_dll()
#devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
devtools::load_all("~/demulticoder")
library("demulticoder")
devtools::document()
pkgdown::build_site()
pkgdown::preview_site()

usethis::use_github_links(overwrite=TRUE)
pkgdown::build_favicons(overwrite=TRUE)

urlchecker::url_check()

devtools::build_readme()

devtools::check(remote = TRUE, manual = TRUE)

devtools::check_win_devel()

git push


rcmdcheck::rcmdcheck(args = "--as-cran")

pkgdown::build_articles("~/demulticoder/")

use_release_issue()
usethis::use_cran_comments()

usethis::use_article("vignettes/Getting_started.Rmd")
usethis::use_article("vignettes/Documentation.Rmd")
usethis::use_article("vignettes/DADA2_16S_mothur_validation.Rmd")

devtools::build()
devtools::check()
devtools::document()
usethis::use_release_issue()
usethis::use_news_md()
########

output<-prepare_reads(
  data_directory = system.file("extdata", package = "demulticoder"),
  output_directory = "~/demulticoder_test17", 
  overwrite_existing = TRUE
  )

cut_trim(
  output,
  cutadapt_path = "/usr/bin/cutadapt", #CHANGE LOCATION TO YOUR LOCAL INSTALLATION
  overwrite_existing = FALSE)

make_asv_abund_matrix(
  output,
  overwrite_existing = TRUE)

assign_tax(
  output,
  asv_abund_matrix,
  retrieve_files=TRUE,
  overwrite_existing = TRUE)

objs<-convert_asv_matrix_to_objs(output, minimum_bootstrap = 0, save_outputs = TRUE)


