% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/step1_linear_clump.R
\name{step1_plink}
\alias{step1_plink}
\title{Run Linear Model with PLINK2}
\usage{
step1_plink(
  pheno_prepared = NULL,
  covar_prepared = NULL,
  file_calls = NULL,
  file_variant_info = NULL,
  file_sample_info = NULL,
  thres_hwe = 1e-20,
  thres_maf = 0.01,
  thres_geno = 0.05,
  ci = 0.95,
  plink_file_version,
  output_name = "out"
)
}
\arguments{
\item{pheno_prepared}{REQUIRED. The path to the phenotype file, prepared by merge_data() which has been filtered for the selected IIDs.}

\item{covar_prepared}{REQUIRED. The path to the covariate file, prepared by merge_data() which has been filtered for the selected IIDs.}

\item{file_calls}{REQUIRED. A string representing the genotype call file, ending with .pgen or .bed.}

\item{file_variant_info}{REQUIRED. A string representing the variant information files, ending with .pvar or .bim.}

\item{file_sample_info}{REQUIRED. A string representing the sample information file, ending with .psam or .fam.}

\item{thres_hwe}{OPTIONAL. Threshold for Hardy-Weignberg equilibrium, default is 1e-20.}

\item{thres_maf}{OPTIONAL. Threshold for MAF of the analyzed SNPs, default is 0.01.}

\item{thres_geno}{OPTIONAL. Threshold of missing call rates for filtering out variants. Default is 0.05.}

\item{ci}{OPTIONAL. Confidence intervals with the given width to be reported for each beta or odds-ratio. Default is 0.95.}

\item{plink_file_version}{REQUIRED. Version of plink files, default is 2}

\item{output_name}{OPTIONAL. A string of the name for the output file, can be the group name or chromosome number.
Default is "out". The suffix "glm.linear" will be added automatically by PLINK.}
}
\value{
GWAS outcomes as text files, generated by PLINK2.
}
\description{
This function runs the linear model GWAS, as Step 1, using PLINK2 with the unextracted genotype files and the prepared covariate and phenotype files
and outputs the GWAS outcomes
}
