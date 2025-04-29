test_that("prepare_reads works correctly", {
  analysis_setup <- prepare_reads(data_directory = system.file("extdata", package = "demulticoder"), overwrite_existing = TRUE)
  expect_type(analysis_setup, "list")  
  expect_true(length(analysis_setup) > 0) 
})
