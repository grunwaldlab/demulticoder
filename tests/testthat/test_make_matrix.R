if (Sys.which("cutadapt") != "") {
  test_that("make_asv_abund_matrix works correctly", {
    data_directory <- system.file("extdata", package = "demulticoder")
    analysis_setup <- prepare_reads(data_directory, overwrite_existing = TRUE)
    cut_trim(analysis_setup, cutadapt_path = "/usr/bin/cutadapt", overwrite_existing = TRUE)
    make_asv_abund_matrix(analysis_setup, overwrite_existing = TRUE)
    rdata_files <- list.files(analysis_setup$directory_paths$temp_directory, pattern = "^asvabund_matrixDADA2_.*\\.RData$", full.names = TRUE)
    asvplots <- list.files(analysis_setup$directory_paths$output_directory, pattern = "^asv_seqlength_plot_.*\\.RData$", full.names = TRUE)
    expect_true(all(file.exists(rdata_files)))
    expect_true(all(file.exists(asvplots)))
  })
} else {
  message("cutadapt is not installed, skipping tests.")
}