test_that("convert_asv_matrix works correctly", {
  data_directory <- system.file("extdata", package = "demulticoder")
  analysis_setup <- prepare_reads(data_directory, overwrite_existing = TRUE)
  cut_trim(analysis_setup, cutadapt_path = "/usr/bin/cutadapt", overwrite_existing = TRUE)
  make_asv_abund_matrix(analysis_setup, overwrite_existing = TRUE)
  assign_tax(analysis_setup, overwrite_existing = TRUE)
  converted_objs <- convert_asv_matrix_to_objs(analysis_setup)
  expect_type(converted_objs, "list")  
  expect_true(length(analysis_setup) > 0) 
})

