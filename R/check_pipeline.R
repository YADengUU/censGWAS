#' Test the pipeline of censGWAS() without running PLINK analyses
#'
#' This function tests whether the pipeline of censGWAS() is designed correctly, by taking different input files,
#' except for the genotype files in PLINK format. To ensure efficiency, GWAS or clumping will not be conducted.
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
#' @param participants REQUIRED. A string of the file including the list of participant IIDs, organized in one column, with header 'IID' (no quote).
#' @param gwas_thres_hwe OPTIONAL. Threshold for Hardy-Weignberg equilibrium, default is 1e-20.
#' @param gwas_thres_maf OPTIONAL. Threshold for MAF of the analyzed SNPs, default is 0.01.
#' @param gwas_thres_geno OPTIONAL. Threshold of missing call rates for filtering out variants. Default is 0.05.
#' @param gwas_ci OPTIONAL. Confidence intervals with the given width to be reported for each beta or odds-ratio. Default is 0.95.
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
#' @param file_snps_chosen OPTIONAL.A string representing the file listing the chosen SNPs in one column, if user directly wants to test certain SNPs.
#' @param output_name OPTIONAL. A string of the name for the output file, can be the group name or chromosome number.
#' Default is "out". The corresponding suffices will be added automatically by the program or PLINK2.
#' @param delete_unconverted OPTIONAL. Boolean (TRUE or FALSE) to indicate whether remove the extracted but not converted intermediate genotype files or not.
#' Default is TRUE
#' @param standardize_covariates OPTIONAL. Boolean (TRUE or FALSE) to indicate whether the quantitative covariates need to be standardized. Default is FALSE.
#' PLINK2 may fail if the covariates have not been standardized. When PLINK2 reports the error, turn this flag to be TRUE. Or, provide your already standardized covariates instead.
#' @returns A list of steps that has been reached to.
check_pipeline <- function(pfile = NULL, bfile = NULL, pgen = NULL, pvar = NULL, psam = NULL, bed = NULL, bim = NULL, fam = NULL,
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
        return("Problematic Output Directory")
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
    return("No PLINK files")
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
    return("Problematic genotype files")
  }

  if (!grepl("\\.(pgen|bed)$", file_calls)){
    warning("Error: genotype calls should be in .pgen or .bed format.")
    return("Error: genotype calls should be in .pgen or .bed format.")
  }
  if (!grepl("\\.(pvar|bim)$", file_variant_info)) {
    warning("Error: variant info should be in .pvar or .bim format.")
    return("Error: variant info should be in .pvar or .bim format.")
  }
  if (!grepl("\\.(psam|fam)$", file_sample_info)) {
    warning("Error: sample info should be in '.psam' or '.fam' format.")
    return("Error: sample info should be in '.psam' or '.fam' format.")
  }


  # identify PLINK version
  plink_file_version <- if(grepl("\\.pgen$", file_calls)) 2 else 1

  step_list <- c() # take record of which steps have been run: 1, 1.5, 2

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

  # step 1 linear model
  if (need_gwas){
    message("Reaching the step of running GWAS. Skipping the actual analyses during the test.")
    gwas_results <- paste0(output_name, ".glm.linear")
    step_list <- c(step_list, 1)
  }

  # step 1.5 clumping
  # if GWAS file was provided, but the file_snps_chosen was not, we also need clumping
  need_clumping <- if(is.null(file_snps_chosen)) TRUE else FALSE
  if (need_clumping){

    if (!need_gwas){

      if(is.null(check_input(gwas_results, "GWAS Result File"))){
        return("Don't need GWAS but have problems in GWAS result files")
      }

      message("GWAS results are provided, but SNPs have not been selected. Now perform clumping ...")

    }
    message("Reaching the step of clumping. Skipping the actual clumping during the test.")
    step_list <- c(step_list, 1.5)

  } else{
    # if the file of selected snps are provided (do not need clumping), we need to check if it is empty and if there are other problems
    message("Has been directly provided the SNPs to test, proceeeding.")
  }

  # step 2 refinement
  # extract and convert the genotypes
  message("Proceeeding to refinement. Actually analyses is not run in testing pipeline.")
  step_list <- c(step_list, 2)
  return(step_list)
}


