if (Sys.which("cutadapt") != "") {
  test_that("cut_trim works correctly", {
    data_directory <- system.file("extdata", package = "demulticoder")
    analysis_setup <- prepare_reads(data_directory, overwrite_existing = TRUE)
    cut_trim(analysis_setup,
             cutadapt_path = "/usr/bin/cutadapt", overwrite_existing = TRUE)
    expect_true(all(file.exists(analysis_setup$data_tables$cutadapt_data$trimmed_path)))
  })
} else {
  message("cutadapt is not installed, skipping tests.")
}


