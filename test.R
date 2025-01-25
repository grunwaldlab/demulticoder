#devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
rcmdcheck::rcmdcheck(args = "--as-cran")
#devtools::clean_dll()
#devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
devtools::load_all("~/demulticoder")
library("demulticoder")
devtools::document()

pkgdown::build_site()

