#devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
rcmdcheck::rcmdcheck(args = "--as-cran")
#devtools::clean_dll()
#devtools::install_github("grunwaldlab/demulticoder", force=TRUE)
devtools::load_all("~/demulticoder")
library("demulticoder")
devtools::document()

output<-prepare_reads(
  data_directory = system.file("extdata", package = "demulticoder"),
  output_directory = "~/output_test", 
  tempdir_id = "test_dataset",
  overwrite_existing=TRUE)

cut_trim(
  output,
  cutadapt_path="/opt/homebrew/bin/cutadapt",
  overwrite_existing = TRUE)

make_asv_abund_matrix(
  output,
  overwrite_existing = TRUE)

assign_tax(
  output,
  overwrite_existing = TRUE)

objs<-convert_asv_matrix_to_objs(output, save_outputs=TRUE, overwrite=TRUE, minimum_bootstrap = 0)

obj<-objs$taxmap_its
print(obj)
obj$data$otu_table
obj$data$tax_data <- metacoder::zero_low_counts(obj, data = "otu_table", min_count = 5)
obj$data$tax_data <-metacoder::calc_obs_props(obj, "tax_data")

obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "tax_data",
                                                  cols = obj$data$sample_data$samplename_barcode)

obj$data$tax_occ <- metacoder::calc_n_samples(obj, "tax_abund", groups = obj$data$sample_data$organism, cols = obj$data$sample_data$samplename_barcode)


obj$data$diff_table <- metacoder::compare_groups(obj,
                                                 data = "tax_abund",
                                                 cols = obj$data$sample_data$samplename_barcode,
                                                 groups = obj$data$sample_data$organism)
print(obj$data$diff_table)

set.seed(999)
metacoder::heat_tree(obj, 
                     node_label = taxon_names,
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                     node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
                     node_color_range = c("cyan", "gray", "tan"), # The color palette used
                     node_size_axis_label = "OTU count",
                     node_color_axis_label = "Log 2 ratio of median proportions",
                     layout = "davidson-harel", # The primary layout algorithm
                     initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations


devtools::check()
usethis::use_vignette("my-vignette")
usethis::use_pkgdown()
pkgdown::build_site()
usethis::use_pkgdown_github_pages()

usethis::use_logo("man/figures/demulticoder_logo_iteration1.png")
