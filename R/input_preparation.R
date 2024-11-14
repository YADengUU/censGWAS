
#' Retrieve Genotype Data
#'
#' This function takes PLINK genotype dosage files (.pgen/.bed, .pvar/.bim, .psam/.fam), participant list,
#' and table of desired SNPs,
#' and returns the genotype files (.raw) readable by R.
#'
#' @param snps_chosen REQUIRED. A dataframe of the SNPs to extract
#' with 1st column including the chromosome numbers (integer 1~22, string 'X', 'XY'),
#' 2nd column for the SNP IDs (string of rsid),
#' and 3rd column for the effect alleles (character 'A', 'T', 'G', 'C').
#' @param participants REQUIRED. A string representing the file for the FID and IID, listed in two columns, of the participants to extract.
#' @param file_calls REQUIRED. A string representing the genotype call file, ending with .pgen or .bed.
#' @param file_variant_info REQUIRED. A string representing the variant information files, ending with .pvar or .bim.
#' @param file_sample_info REQUIRED. A string representing the sample information file, ending with .psam or .fam.
#' @param plink_file_version REQUIRED. Integer indicating the version PLINK, 1 or 2. Generated automatically by identifying the genotype file suffices.
#' @param output_name OPTIONAL. A string of the name for the output file, can be the group name or chromosome number.
#' Default is "out". The suffix is added automatically by PLINK.
#' @param delete_unconverted OPTIONAL. Boolean (TRUE or FALSE) to indicate whether remove the extracted but not converted intermediate files or not.
#' Default is TRUE
#' @importFrom utils write.table
#' @returns Genotype files: extracted (.bed or .pgen along with the companions) and converted (.raw) for R for subsequent computations,
#' as well as the a value 0 or NULL to indicate the successfulness of this stage.
#' @export
get_geno <- function(snps_chosen, participants = NULL, file_calls = NULL, file_variant_info = NULL, file_sample_info = NULL,
                     plink_file_version, output_name = "out", delete_unconverted = TRUE){

  message("Trying to extract genotype files for Step 2 ...")

  if (output_name != "out"){
    output_dir <- dirname(output_name)
    if (output_dir=="."){
      output_dir <- getwd()
    }
  }else{
    output_dir <- getwd()
  }


  if (!(is.logical(delete_unconverted) && length(delete_unconverted)==1)){
    warning("Invalid input for delete_unconverted, should be boolean TRUE/FALSE.")
    return(NULL)
  }

  snps <- snps_chosen[, 2]
  current_snps <- paste0(output_name, "_snps.txt")
  # write snp id
  write.table(snps, file = current_snps, row.names = FALSE, col.names=FALSE, quote = FALSE, sep = "\t")

  # select the corresponding PLINK version to extract genotypes
  if (plink_file_version == 2){
    # use pgen
    plink_cmd1 <- paste0("plink2 --pgen ", file_calls, " --psam ", file_sample_info, " --pvar ", file_variant_info, " --keep ", participants, " --extract ", current_snps, " --make-pgen --out ", output_name)
  }else{
    # use bed
    plink_cmd1 <- paste0("plink2 --bed ", file_calls, " --fam ", file_sample_info, " --bim ", file_variant_info, " --keep ", participants, " --extract ", current_snps, " --make-bed --out ", output_name)
  }

  exit_status <- system(plink_cmd1)
  if (exit_status !=0){
    warning("Error: PLINK2 encountered an error and exited with status ", exit_status, ". SNP extraction failed.")
    return(NULL)
  }
  else{
    if (output_dir == "."){
      message("SNPs successfully extracted. Results saved to current working directory ", output_dir, ".")
    }else{
      message("SNPs successfully extracted. Results saved to ", output_dir, ".")
    }
  }

  # specific on the effect allele A1
  snps_with_a1 <- snps_chosen[, c(2, 3)]
  write.table(snps_with_a1, file = current_snps, row.names = FALSE, col.names=FALSE, quote = FALSE, sep = "\t")
  # make the .raw file for R to read
  if (plink_file_version==2){
    plink_cmd2 <- paste0("plink2 --pfile ", output_name, " --export A --export-allele ", current_snps, " --out ", output_name)
  }else{
    plink_cmd2 <- paste0("plink2 --bfile ", output_name, " --export A --export-allele ", current_snps, " --out ", output_name)
  }

  exit_status <- system(plink_cmd2)
  if (exit_status !=0){
    warning("Error: PLINK2 encountered an error and exited with status ", exit_status, ". Genotype format conversion to .raw failed.")
    return(NULL)
  }
  else{
    message("Genotype format successfully converted to .raw. Results saved to ", output_dir, ".")
  }

  # delete the intermediate unconverted genotype files if opt to
  if (delete_unconverted){

    if (plink_file_version==2){
      file1 <- paste0(output_name, ".pgen")
      file2 <- paste0(output_name, ".psam")
      file3 <- paste0(output_name, ".pvar")
    }else{
      file1 <- paste0(output_name, ".bed")
      file2 <- paste0(output_name, ".fam")
      file3 <- paste0(output_name, ".bim")
    }

    files_to_remove <- c(file1, file2, file3)
    files_exist <- file.exists(files_to_remove)
    if (any(files_exist)){
      file.remove(files_to_remove[files_exist])
      message("Intermediate (unconverted) genotype files removed: ", paste(files_to_remove[files_exist], collapse = ", "))
    }else{
      message("No intermediate (unconverted) genotype files found for removal. Potential error at the extraction step.")
    }
  }
  # also delete the intermediate snp list
  file.remove(current_snps)
  return(0)

}


#' Check Input Files
#'
#' This function takes the input paths and checks if any paths are missing or unavailable.
#'
#' @param file_name A string of the path to the file to be checked for availability
#' @param file_type A string of the specific type of file to be checked,
#' including "SNPs for Extraction", "Participant IDs for Extraction", "Genotype Calls (.pgen/.bed)", "Variant Information (.pvar/.bim)", "Sample Information (.psam/.fam)",
#' "Output Directory", "Covariate File", "Phenotype File", "GWAS Result File", "Clumped Result File"
#' @returns A value 0 or NULL indicating if all criteria are met and procedures can continue.
check_input <- function(file_name, file_type){

  if (file_type == "Output Directory"){
    if (!file.exists(file_name)){
      warning("Error: The directory for the output is not valid.")
      return(NULL)
    }
  }else {
    if (missing(file_name) | !file.exists(file_name)){
      warning("Error: ", file_type, " is not valid or missing.")
      return(NULL)
    }
  }
  return(0)

}

#' Check Format of SNP List
#'
#' This function checks format of the tab-/space-separated file including SNPs to extract, if directly provided by the user
#' and returns the converted dataframe if it is correctly formatted with the three columns (Chr., ID, refAllele).
#'
#' @param snp_df The dataframe read from the text (tab-/space-separated) file containing the SNPs to extract
#' @returns Dataframe of the list of selected SNPs, if correctly formatted; otherwise, NULL.
#' @importFrom utils read.delim2
#' @export
check_snp_list <- function(snp_df){
  # check number of columns, the number of rows will be check in main.R
  if (ncol(snp_df) != 3){
    warning("Error: The list of SNPs is incorrectly formmated, should be three columns with chr no., SNP RSID, effect allele")
    return(NULL)
  }

  # tidy some format
  colnames(snp_df) <- c("chr", "ID", "A1")
  valid_chr <- !is.na(snp_df$chr) & snp_df$chr %in% as.character(c(1:23, 'X', 'XY', 'x', 'xy', 'Y', 'y'))
  snp_df$A1 <- toupper(snp_df$A1)

  # check if the data type of each column is correct
  if (!all(valid_chr)){
    warning("Error: The chromosome no. in the SNP list must be 1~23, X, XY (case insensitive), and non-missing.")
    return(NULL)
  }
  if (!all(!is.na(snp_df$A1) & snp_df$A1 %in% c("A", "C", "G", "T"))){
    warning("Error: The reference bases must be A, T, C, or G.")
    return(NULL)
  }

  # output the dataframe if everything looks fine
  return(snp_df)

}

#' Merge Covariate and Phenotype Data Files
#'
#' This function reads the covariate and phenotype files and participant ID
#' and returns the organized dataframe after checking for complete cases
#'
#' @param file_covar REQUIRED. A string of the path to the covariate file (tab- or space-separated) containing the FID, IID, and covariate columns.
#' @param file_pheno REQUIRED. A string of the path to the phenotype file (tab- or space-separated) containing the FID, IID, and censored status (0 for below /1 for above LOD), measurements.
#' @param participants REQUIRED. A string representing the file for the IIDs, listed in a column, of the participants to extract.
#' @param output_name OPTIONAL. A string of the prefix name to include for the output file, can be the group or strata name.
#' Default is "out". Function automatically adds suffix ".covar", ".pheno".
#' @param need_stdz OPTIONAL. Boolean indicating whether standardization of quantitative covariates is needed. Default is FALSE.
#' @param run_linear REQUIRED. Boolean indicating whether to run linear model with PLINK, defined automatically based on if the gwas_results are provided. If TRUE, the phenotype and covariates are outputted as TSV files.
#' @returns Dataframe of the covariates and phenotypes (FID, IID, binary censoring status, measurements before and after censoring, rank-transformed values for linear and Tobit, and the covariates), with complete cases,
#' along with the corresponding tsv files for PLINK use for running linear model, if chosen to run.
#' If the procedures were not able to continue, returns NULL.
#' @importFrom utils read.table write.table
#' @importFrom stats complete.cases
#' @importFrom RNOmni RankNorm
#' @export
merge_data <- function(file_covar = NULL, file_pheno = NULL, participants = NULL, output_name = "out", need_stdz = FALSE, run_linear){

  iid_chosen <- read.table(participants, header = TRUE)
  iid_chosen <- iid_chosen$IID

  message("Checking covariate and phenotype file validity ...")
  if (is.null(check_input(file_covar, "Covariate File")) | is.null(check_input(file_pheno, "Phenotype File"))){
    return(NULL)
  }


  if (output_name != "out"){
    output_dir <- dirname(output_name)
    if (output_dir=="."){
      output_dir <- getwd()
    }
  }else{
    output_dir <- getwd()
  }


  covars <- read.table(file_covar, header = TRUE)
  phenos <- read.table(file_pheno, header = TRUE)
  message("Checking covariate sample availability ...")
  covars_chosen <- check_FID_IID(covars, iid_chosen)
  if (is.null(covars_chosen)){
    return(NULL)
  }
  message("Checking phenotype sample availability ...")
  phenos_chosen <- check_FID_IID(phenos, iid_chosen)
  if (is.null(phenos_chosen)){
    return(NULL)
  }
  rm(covars, phenos) # remove original entries

  # check case completeness
  covars_complete <- covars_chosen[complete.cases(covars_chosen), ]
  # for phenotypes, need to censor the truncated data if they haven't been done so
  if (any(is.na(phenos_chosen[phenos_chosen[,3]==0, 4]))){
    message("The provided measurements are truncated, now censoring with the LOD.")
    measure_temp <- phenos_chosen[,4]
    measure_temp[which(phenos_chosen[,3]==0)] <- min(phenos_chosen, na.rm=TRUE)
    phenos_chosen[, 4] <- measure_temp
  }

  phenos_complete <- phenos_chosen[complete.cases(phenos_chosen), ]
  message(nrow(covars_chosen), " samples in covariate entries; ", nrow(covars_complete), " have complete cases.")
  message(nrow(phenos_chosen), " samples in phenotype entries; ", nrow(phenos_complete), " have complete cases.")
  rm(covars_chosen, phenos_chosen)

  # begin merging dataframes by merging the IIDs
  covar_ids <- covars_complete$IID
  pheno_ids <- phenos_complete$IID
  common_ids <- intersect(covar_ids, pheno_ids)
  message(length(common_ids), " samples are available for both covariate and phenotype data.")
  covars_complete <- covars_complete[covars_complete$IID %in% common_ids, ]
  phenos_complete <- phenos_complete[phenos_complete$IID %in% common_ids, ]

  # process the covariate and phenotype files
  # standardize continuous variables in covariates if indicated
  if (need_stdz){
    message("Standardizing the numerical covariates...")
    numeric_cols <- sapply(covars_complete[, 3:ncol(covars_complete)], is.numeric)
    covars_complete[, 3:ncol(covars_complete)][, numeric_cols] <-
      scale(covars_complete[, 3:ncol(covars_complete)][, numeric_cols])
  }
  colnames(phenos_complete) <- c("FID", "IID", "cens_status", "censored")
  # minumum measurement becomes the transformed LOD
  phenos_complete$tf_tb <- RankNorm(phenos_complete$censored, ties.method = "max") # 5 columns

  if (run_linear){
    # inverse normal rank-based transformation, minimum measurement becomes the imputed average below LOD
    phenos_complete$tf_lin <- RankNorm(phenos_complete$censored, ties.method = "average") # 6 columns
  }

  write.table(covars_complete, file = paste0(output_name, ".covar"), row.names = FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
  write.table(phenos_complete, file = paste0(output_name, ".pheno"), row.names = FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
  return(pheno_covar=merge(phenos_complete, covars_complete[, 2:ncol(covars_complete)], by = "IID"))

}

#' Check Covariate and Phenotype File First Two Columns and Sample Availability
#'
#' This function reads and checks the covariate and phenotype files if they have the 1st and 2nd columns as FID and IID,
#' and if the input contains less cases than the desired sample, if the input contains the cases for the desired sample
#'
#' @param df Dataframe of the phenotype or covariates
#' @param iid_chosen List of IIDs chosen for analyses
#' @returns Dataframe containing the samples with selected IIDs, if available; otherwise, NULL.
check_FID_IID <- function(df, iid_chosen){
  if (!all(colnames(df)[1:2] == c("FID", "IID"))){
    warning("Error: The first two columns of phenotype/covariate file must be FID and IID.")
    return(NULL)
  }

  if (length(iid_chosen) > nrow(df)){
    warning("Warning: The input data has less cases than your chosen sample.")
    return(NULL)
  }

  iid_exist <- iid_chosen %in% df$IID
  num_exist <- sum(iid_exist)
  if (num_exist == 0){
    warning("Error: None of the selected IIDs are present in the input data.")
    return(NULL)
  }

  df_filtered <- df[df$IID %in% iid_chosen, ]
  message(num_exist, " of ", length(iid_chosen), " chosen samples are present in the input data with ", nrow(df), " sample.")

  return(df_filtered)
}
