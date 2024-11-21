#' Whole Pipeline of Linear-Tobit Scheme
#'
#' This function runs Step 1 GWAS with the linear model, Step 1.5 clumps the GWAS results,
#' and conducts Step 2 refinement with the Tobit model.
#' Please note that PLINK2 is used to run the linear model, make sure it is installed in the current computational environment.
#' Also, PLINK2 does not account for genetic relatedness, so another tool (SAIGE, REGENIE, etc.) could be used.
#' The three options for the pipeline (Step 1 + Clumping + Step 2; Clumping + Step 2, given that GWAS is run by the user themselves; Step 2 only, given that the SNPs have been clumped/selected by the user)
#' will be chosen to run automatically based on the users' input files. See examples below.
#' @param pfile OPTIONAL. String. Name of the PLINK2 genotype files, if .pgen, .pvar, .psam have the same name, no suffix needed in this field. For example, with sample.pgen, sample.pvar, sample.psam, simply provide "sample" for this entry.
#' @param bfile OPTIONAL. String. Name of the PLINK1 genotype files, if .bed, .bim, .fam have the same name, no suffix needed in this field. For example, with sample.bed, sample.bim, sample.fam, simply provide "sample" for this entry.
#' @param pgen OPTIONAL. A string representing the genotype call file in PLINK2 format, ending with .pgen, should be provided along with pvar and psam in order for the program to proceed.
#' @param pvar OPTIONAL. A string representing the variant information file, ending with .pvar.
#' @param psam OPTIONAL. A string representing the sample information file, ending with .psam.
#' @param bed OPTIONAL. A string representing the genotype call file in PLINK1 format, ending with .bed, should be provided along with bim and fam in order for the program to proceed.
#' @param bim OPTIONAL. A string representing the variant information file, ending with .bim.
#' @param fam OPTIONAL. A string representing the sample information file, ending with .fam.
#' @param file_pheno REQUIRED.A string of the path to the phenotype file (tab- or space-separated) containing the FID, IID, and censored status (0 for below /1 for above LOD), measurements.
#' @param file_covar REQUIRED. A string of the path to the covariate file (tab- or space-separated) containing the FID, IID, and covariate columns.
#' @param participants REQUIRED. A string of the file including the list of participant FIDs and IIDs, organized in 2 column, with header 'FID', 'IID' (no quote).
#' @param gwas_thres_hwe OPTIONAL. Threshold for Hardy-Weignberg equilibrium, default is 1e-20.
#' @param gwas_thres_maf OPTIONAL. Threshold for MAF of the analyzed SNPs, default is 0.01.
#' @param gwas_thres_geno OPTIONAL. Threshold of missing call rates for filtering out variants. Default is 0.05.
#' @param gwas_ci OPTIONAL. Confidence intervals with the given width to be reported for each beta. Default is 0.95.
#' @param gwas_results OPTIONAL. A string of the file including GWAS outcomes, if already performed by the user.
#' @param snp_id_field OPTIONAL. The field (column) name in the GWAS result files that contains the SNP IDs, default is 'ID'.
#' @param p_field OPTIONAL. The field (column) name in the GWAS result files that contains the p-values, default is 'P'.
#' @param a1_field OPTIONAL. The field (column) name in the GWAS result files that contain the effect allele. Default is 'A1'.
#' @param a1_freq_field OPTIONAL. The field (column) name in the GWAS result files that contain the A1 allele frequency. Defeult is 'A1_FREQ'.
#' @param beta_field OPTIONAL. The field (column) name in the GWAS result files that contain the effect estimate (slope). Default is 'BETA'.
#' @param se_field OPTIONAL. The field (column) name in the GWAS result files that contain the standard error. Default is 'SE'.
#' @param pos_field OPTIONAL. The field (column) name in the GWAS result files that contain the SNP position. Default is 'POS'.
#' @param clump_p1 OPTIONAL. The index variant p-value threshold for clumping. Default is 5e-8.
#' @param clump_p2 OPTIONAL. The SP2 column p-value threshold for clumping. Default is 1e-4.
#' @param clump_r2 OPTIONAL. The r^2 (correlation between genetic variants to reduce redundancy) threshold for clumping. Default is 0.1.
#' @param clump_kb OPTIONAL. The maximum distance from the lead variant for SNPs to be considered to be clumped with. Default is 1000.
#' @param file_snps_chosen OPTIONAL.A string representing the file listing the chosen SNPs with the three columns (Chr., ID, Allele), if user directly wants to test certain SNPs.
#' @param output_name OPTIONAL. A string of the name for the output file, can be the group name or chromosome number.
#' Default is "out". The corresponding suffices will be added automatically by the program or PLINK2.
#' @param delete_unconverted OPTIONAL. Boolean (TRUE or FALSE) to indicate whether remove the extracted but not converted intermediate genotype files or not.
#' Default is TRUE
#' @param standardize_covariates OPTIONAL. Boolean (TRUE or FALSE) to indicate whether the quantitative covariates need to be standardized. Default is FALSE.
#' PLINK2 may fail if the covariates have not been standardized. When PLINK2 reports the error, turn this flag to be TRUE. Or, provide your already standardized covariates instead.
#' @importFrom utils write.table read.table
#' @returns A plaintext file containing the refined estimates for the SNPs.
#' @export
censGWAS <- function(pfile = NULL, bfile = NULL, pgen = NULL, pvar = NULL, psam = NULL, bed = NULL, bim = NULL, fam = NULL,
                     file_pheno = NULL, file_covar = NULL, participants = NULL,
                     gwas_thres_hwe = 1e-20, gwas_thres_maf = 0.01, gwas_thres_geno = 0.05, gwas_ci = 0.95, gwas_results = NULL,
                     snp_id_field = "ID", p_field = "P", a1_field = "A1", a1_freq_field = "A1_FREQ", beta_field = "BETA", se_field = "SE", pos_field = "POS",
                     clump_p1 = 5e-8, clump_p2 = 1e-4, clump_r2 = 0.1, clump_kb = 1000,
                     file_snps_chosen = NULL, output_name = "out", delete_unconverted = TRUE, standardize_covariates = FALSE){

  # check output directory
  message("Checking output directory validity ...")
  if (output_name != "out"){
    output_dir <- dirname(output_name)
    if (output_dir != "."){
      if(is.null(check_input(output_dir, "Output Directory"))){
        return(NULL)
      }
    }else{
      output_dir <- getwd()
    }
  }else{
    output_dir <- getwd()
  }

  # whichever pipeline the user chooses, the plink files need to be checked for validity
  if (is.null(pfile) & is.null(bfile) & is.null(pgen) & is.null(pvar) & is.null(psam) & is.null(bed) & is.null(bim) & is.null(fam)){
    warning("Error: No PLINK genotype files are provided.")
    return(NULL)
  }

  if (!is.null(pfile)){
    file_calls <- paste0(pfile, ".pgen")
    file_variant_info <- paste0(pfile, ".pvar")
    file_sample_info <- paste0(pfile, ".psam")
  }

  if (!is.null(bfile)){
    file_calls <- paste0(bfile, ".bed")
    file_variant_info <- paste0(bfile, ".bim")
    file_sample_info <- paste0(bfile, ".fam")
  }

  if (!is.null(pgen)){
    file_calls <- pgen
  }

  if(!is.null(pvar)){
    file_variant_info <- pvar
  }

  if (!is.null(psam)){
    file_sample_info <- psam
  }

  if(!is.null(bed)){
    file_calls <- bed
  }

  if(!is.null(bim)){
    file_variant_info <- bim
  }

  if(!is.null(fam)){
    file_sample_info <- fam
  }

  message("Checking validity of the input PLINK genotype files ...")

  if (is.null(check_input(file_calls, "Genotype Calls (.pgen/.bed)")) | is.null(check_input(file_variant_info, "Variant Information (.pvar/.bim)")) | is.null(check_input(file_sample_info, "Sample Information (.psam/.fam)"))){
    return(NULL)
  }

  if (!grepl("\\.(pgen|bed)$", file_calls)){
    warning("Error: genotype calls should be in .pgen or .bed format.")
    return(NULL)
  }
  if (!grepl("\\.(pvar|bim)$", file_variant_info)) {
    warning("Error: variant info should be in .pvar or .bim format.")
    return(NULL)
  }
  if (!grepl("\\.(psam|fam)$", file_sample_info)) {
    warning("Error: sample info should be in '.psam' or '.fam' format.")
    return(NULL)
  }


  # identify PLINK version
  plink_file_version <- if(grepl("\\.pgen$", file_calls)) 2 else 1

  # read phenotype and covariate files
  # if gwas_results is provided, it does not need to run GWAS here; same for if the snps have already been chosen
  if (!is.null(file_snps_chosen)){
    need_gwas <- FALSE
  }else{
    if(!is.null(gwas_results)){
      need_gwas <- FALSE
    }else{
      need_gwas <- TRUE
    }
  }

  pheno_covar <- merge_data(file_covar = file_covar, file_pheno = file_pheno, participants = participants, output_name = output_name, need_stdz = standardize_covariates, run_linear = need_gwas)
  if (is.null(pheno_covar)){
    return(NULL)
  }
  # now the participant list is updated,
  participants_list <- data.frame(FID =pheno_covar[,1], IID=pheno_covar[,2])
  # output the new participant list
  write.table(participants_list, file = paste0(output_name, ".participants"), row.names = FALSE, col.names=TRUE, quote = FALSE, sep = "\t")
  participants_update <- paste0(output_name, ".participants") # use later for the genotype extraction

  # the program should look for the .covar and .pheno files in the output directory
  base_name <- basename(output_name)
  file_covar_prepared <- list.files(path = output_dir, pattern = paste0(base_name, ".*\\.covar$"), full.names = TRUE)
  file_pheno_prepared <- list.files(path = output_dir, pattern = paste0(base_name, ".*\\.pheno$"), full.names = TRUE)

  # step 1 linear model
  if (need_gwas){
    if(is.null(step1_plink(pheno_prepared = file_pheno_prepared, covar_prepared = file_covar_prepared, file_calls = file_calls, file_variant_info = file_variant_info, file_sample_info = file_sample_info,
                           thres_hwe = gwas_thres_hwe, thres_maf = gwas_thres_maf, thres_geno = gwas_thres_geno, ci = gwas_ci, plink_file_version = plink_file_version, output_name = output_name))){
      return(NULL)
    }
    # at this point we have the gwas results, the names are automatically set with the suffix .glm.linear
    #gwas_results <- paste0(output_name, ".tf_lin.glm.linear")
    gwas_results <- list.files(path = output_dir, pattern = paste0("^", base_name, ".*\\.glm\\.linear$"), full.names = TRUE)
  }

  # step 1.5 clumping
  # if GWAS file was provided, but the file_snps_chosen was not, we also need clumping
  need_clumping <- if(is.null(file_snps_chosen)) TRUE else FALSE
  if (need_clumping){

    if (!need_gwas){

      if(is.null(check_input(gwas_results, "GWAS Result File"))){
        return(NULL)
      }

      message("GWAS results are provided, but SNPs have not been selected. Now perform clumping ...")

    }
    clumped_snps <- clump_gwas(gwas_results = gwas_results, p1 = clump_p1, p2 = clump_p2, r2 = clump_r2, kb = clump_kb,
                               file_calls = file_calls, file_variant_info = file_variant_info, file_sample_info = file_sample_info,
                               plink_file_version = plink_file_version, output_name = output_name,
                               snp_id_field = snp_id_field, p_field = p_field, a1_field = a1_field, a1_freq_field = a1_freq_field, beta_field = beta_field, se_field = se_field, pos_field = pos_field)

    if(is.null(clumped_snps)){
      warning("Unable to proceed for refinement because no SNPs were left after clumping.")
      return(NULL)
    }
    snps_chosen <- snps_to_extract(clumped_snps)
    snps_to_refine <- clumped_snps# this is for the output df that will be concatenated with the refined results
  } else{
    # if the file of selected snps are provided (do not need clumping), we need to check if it is empty and if there are other problems
    snp_df <- read.table(file_snps_chosen, header = TRUE, ,colClasses = "character")

    if (!nrow(snp_df) > 0){
      warning("Error: The list of SNPs is empty.")
      return(NULL)
    }

    snps_chosen <- check_snp_list(snp_df)
    snps_to_refine <- snps_chosen# this is for the output df that will be concatenated with the refined results

  }

  # step 2 refinement
  # extract and convert the genotypes
  if(is.null(get_geno(snps_chosen = snps_chosen, participants = participants_update, file_calls = file_calls, file_variant_info = file_variant_info, file_sample_info = file_sample_info,
                      plink_file_version = plink_file_version, output_name = output_name, delete_unconverted = delete_unconverted))){
    return(NULL)
  }

  file_geno_converted <- list.files(path = output_dir, pattern = paste0(base_name, ".*\\.raw$"), full.names = TRUE)

  step2_run_tobit(geno_converted = file_geno_converted, pheno_covar = pheno_covar, snps_to_refine = snps_to_refine, output_name = output_name)
  return(0)
}

