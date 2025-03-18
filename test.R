#devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
devtools::clean_dll()
#devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
devtools::load_all("~/demulticoder")
library("demulticoder")
devtools::document()
devtools::document()
install.packages("dplyr")
library("dplyr")

install.packages("pkgdown")
library("pkgdown")
install.packages(c("metacoder", "sessioninfo", "pkgdown"))
library("metacoder")
library("sessioninfo")
library("pkgdown")

pkgdown::build_site()
pkgdown::preview_site()

pkgdown::build_favicons(overwrite=TRUE)

rcmdcheck::rcmdcheck(args = "--as-cran")

pkgdown::build_articles("~/demulticoder/")

usethis::use_article("vignettes/Getting_started.Rmd")
usethis::use_article("vignettes/Documentation.Rmd")
usethis::use_article("vignettes/DADA2_16S_mothur_validation.Rmd")
system(paste("R CMD Rd2pdf --no-clean", "~/demulticoder"))

devtools::build()
devtools::check()
devtools::document()
