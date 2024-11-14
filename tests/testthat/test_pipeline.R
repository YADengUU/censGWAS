# tests/test_pipeline.R

# runs check_pipeline() to see if the process can be reached correctly

test_that("Pipeline (GWAS, clumping, and refinement) can run in order", {
  dummy_data_path <- file.path("..", "dummy_data")
  step_list <- check_pipeline(pfile = paste0(dummy_data_path, "/dummy"), file_pheno = paste0(dummy_data_path, "/dummy.pheno"),
                          file_covar = paste0(dummy_data_path, "/dummy.covar"))
  # by not providing the GWAS results or selected snps, the steps 1, 1.5, 2 should be run
  expect_equal(step_list, c(1, 1.5, 2))

  step_list <- check_pipeline(pfile = paste0(dummy_data_path, "/dummy"), file_pheno = paste0(dummy_data_path, "/dummy.pheno"),
                          file_covar = paste0(dummy_data_path, "/dummy.covar"), gwas_results = paste0(dummy_data_path, "/dummy.glm.linear"))
  # by providing the gwas results but not the selected snps, the steps 1.5, 2 should be run
  expect_equal(step_list, c(1.5, 2))

  step_list <- check_pipeline(pfile = paste0(dummy_data_path, "/dummy"), file_pheno = paste0(dummy_data_path, "/dummy.pheno"),
                          file_covar = paste0(dummy_data_path, "/dummy.covar"), file_snps_chosen = paste0(dummy_data_path, "/dummy_snps.txt"))
  # by providing the selected snps, only the step 2 should be run
  expect_equal(step_list, c(2))
})

