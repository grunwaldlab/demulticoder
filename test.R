#devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
devtools::clean_dll()
devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
devtools::load_all("~/demulticoder")
library("demulticoder")
devtools::document()
devtools::document()
install.packages("dplyr")
library("dplyr")

install.packages("pkgdown")
library("pkgdown")
pkgdown::build_site(check = FALSE)
install.packages(c("metacoder", "sessioninfo", "pkgdown"))
library("metacoder")
library("sessioninfo")
library("pkgdown")

pkgdown::build_site()
pkgdown::preview_site()

pkgdown::build_favicons(overwrite=TRUE)

rcmdcheck::rcmdcheck(args = "--as-cran")