if (Sys.which("cutadapt") != "") {
  test_that("assign_taxonomy works correctly", {
    data_directory <- system.file("extdata", package = "demulticoder")
    analysis_setup <- prepare_reads(data_directory, overwrite_existing = TRUE)
    cut_trim(analysis_setup, cutadapt_path = "/usr/bin/cutadapt", overwrite_existing = TRUE)
    make_asv_abund_matrix(analysis_setup, overwrite_existing = TRUE)
    assign_tax(analysis_setup, overwrite_existing = TRUE)
    final_asv_matrix <- list.files(analysis_setup$directory_paths$output_directory, pattern = "^final_asv_abundance_matrix_*\\.RData$", full.names = TRUE)
    asvplots <- list.files(analysis_setup$directory_paths$output_directory, pattern = "^track_reads_*\\.RData$", full.names = TRUE)
    expect_true(all(file.exists(final_asv_matrix)))
    expect_true(all(file.exists(asvplots)))
  })
} else {
  message("cutadapt is not installed, skipping tests.")
}