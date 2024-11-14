# tests/testthat/test_general.R
# runs check_pipeline() to see if the program can check the availability of PLINK files when the inputs are given differently
# generates simulated genotype, covariate, and phenotype data and runs step2_run_tobit to see if the refinement is done without error
test_that("Can take bfiles as input",{
  dummy_data_path <- file.path("..", "dummy_data")
  step_list <- check_pipeline(bfile = paste0(dummy_data_path, "/dummy"), file_pheno = paste0(dummy_data_path, "/dummy.pheno"),
                              file_covar = paste0(dummy_data_path, "/dummy.covar"))
  expect_equal(step_list, c(1, 1.5, 2))
})

test_that("Can take PLINK files (pgen) separately",{
  dummy_data_path <- file.path("..", "dummy_data")
  step_list <- check_pipeline(pgen = paste0(dummy_data_path, "/dummy.pgen"), pvar = paste0(dummy_data_path, "/dummy.pvar"), psam = paste0(dummy_data_path, "/dummy.psam"),
                              file_pheno = paste0(dummy_data_path, "/dummy.pheno"),
                              file_covar = paste0(dummy_data_path, "/dummy.covar"))
  expect_equal(step_list, c(1, 1.5, 2))
})

test_that("Can take PLINK files (bed) seperately",{
  dummy_data_path <- file.path("..", "dummy_data")
  step_list <- check_pipeline(bed = paste0(dummy_data_path, "/dummy.bed"), bim = paste0(dummy_data_path, "/dummy.bim"), fam = paste0(dummy_data_path, "/dummy.fam"),
                              file_pheno = paste0(dummy_data_path, "/dummy.pheno"),
                              file_covar = paste0(dummy_data_path, "/dummy.covar"))
  expect_equal(step_list, c(1, 1.5, 2))
})

test_that("Gives the corresponding warning message if something invalid with PLINK files", {
  dummy_data_path <- file.path("..", "dummy_data")
  warning_message <- check_pipeline(bed = paste0(dummy_data_path, "/dummyWrong.bed"), bim = paste0(dummy_data_path, "/dummyWrong.bim"), fam = paste0(dummy_data_path, "/dummyWrong.fam"),
                              file_pheno = paste0(dummy_data_path, "/dummy.pheno"),
                              file_covar = paste0(dummy_data_path, "/dummy.covar"))
  expect_equal(warning_message, "Problematic genotype files")
})

test_that("The refinement step can run smoothly with truncated measurements and .raw genotype", {

  # can simulate data
  expect_error({
    simulate_data(cens_prop=0.6, betas=c(0.1,0.2), mafs=c(0.4, 0.3), N=25000, to_truncate=TRUE)
  }, regexp = NA)

  simulate_data(cens_prop=0.6, betas=c(0.1,0.2), mafs=c(0.4, 0.3), N=25000, to_truncate=TRUE)

  # can check SNP list format
  snp_df_raw <- read.table("snps.chosen", header = TRUE)
  expect_error({
    check_snp_list(snp_df_raw)
  }, regexp = NA)
  snp_df <- check_snp_list(snp_df_raw)

  # if the snp list is wrong format
  dummy_snp_path <- file.path("..", "dummy_snp")
  snp_df_dummy <- read.table(paste0(dummy_snp_path, "/snps_wrongChrom.chosen"), header = TRUE)
  expect_warning(
    check_snp_list(snp_df_dummy),
    regexp = "must be 1~22, X, XY", fixed = TRUE
  )
  snp_df_dummy <- read.table(paste0(dummy_snp_path, "/snps_wrongFormat.chosen"), header = TRUE)
  expect_warning(
    check_snp_list(snp_df_dummy),
    regexp = "incorrectly formmated", fixed=TRUE)
  snp_df_dummy <- read.table(paste0(dummy_snp_path, "/snps_wrongNuc.chosen"), header = TRUE)
  expect_warning(
    check_snp_list(snp_df_dummy),
    regexp = "reference bases", fixed=TRUE)


  participant_list <- read.table("examples.pheno", header=TRUE)
  participant_list <- data.frame(FID = participant_list$FID, IID=participant_list$IID)
  write.table(participant_list, file="participantsIID.txt", row.names = FALSE, col.names=TRUE, quote = FALSE, sep = "\t")


  pheno_covar_df <- merge_data(file_covar = "examples.covar", file_pheno = "examples.pheno", participants = "participantsIID.txt", output_name = "example_organize", run_linear=TRUE)

  # refinement can run
  expect_error(
    step2_run_tobit(geno_converted = "examples.raw", pheno_covar = pheno_covar_df, snps_to_refine = snp_df),
    regexp = NA)

})

test_that("The refined estimate is more accurate than the primitive one.", {

  # this time let's try with only one SNP; and we also check that the program can work if censored data (not truncated with NA)
  simulate_data(cens_prop=0.6, betas=c(0.2), mafs=c(0.4), N=25000, to_truncate=FALSE)
  snp_df_raw <- read.table("snps.chosen",header=TRUE)# "snps.chosen" is just the file name default in our simulations
  #print(snp_df_raw)
  snp_df<-check_snp_list(snp_df_raw)
  #print(nrow(snp_df))
  participant_list <- read.table("examples.pheno", header=TRUE)
  participant_list <- data.frame(FID = participant_list$FID, IID=participant_list$IID)
  write.table(participant_list, file="participantsIID.txt", row.names = FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
  pheno_covar_df <- merge_data(file_covar = "examples.covar", file_pheno = "examples.pheno", participants = "participantsIID.txt", output_name = "example_organize", run_linear=TRUE)

  # run linear regression on the simulated data
  df_lin <- data.frame(matrix(nrow = nrow(snp_df), ncol = 6)) # ID, intercept, beta, se, stat, p
  colnames(df_lin) <- c("ID", "INTERCEPT", "BETA_LIN","SE_LIN", "STAT", "P_LIN")
  geno_data <- read.table("examples.raw", header =TRUE)
  geno_data <- geno_data[, -c(1, 3:6)] # exclude the FID as well as columns that are not the SNPs
  colnames(geno_data)[2:ncol(geno_data)] <- sub("_.*", "", colnames(geno_data)[2:ncol(geno_data)])
  merge_for_test <- merge(pheno_covar_df, geno_data, by = "IID")
  linear_model <- lm(as.formula(paste0("tf_lin ~ ", colnames(merge_for_test)[ncol(merge_for_test)], " + COVAR1")), data = merge_for_test)
  linear_organized <- data.frame(tidy(linear_model))
  df_lin[1, 2] <- linear_organized[1, "estimate"]
  df_lin[1, 3] <- linear_organized[2, "estimate"]
  df_lin[1, 4] <- linear_organized[2, "std.error"]
  df_lin[1, 5] <- linear_organized[2, "statistic"]
  df_lin[1, 6] <- linear_organized[2, "p.value"]

  step2_run_tobit(geno_converted = "examples.raw", pheno_covar = pheno_covar_df, snps_to_refine = snp_df, output_name = "examples") # here we also test if it can take another output name
  tobit_result <- read.table("examples.refined", header = TRUE)

  expect_true(abs(tobit_result$BETA_TB[1]-0.2) < abs(df_lin$BETA_LIN[1]-0.2)) # we expect that the error between tobit estimates and assumed slope is smaller than the error of linear


})
