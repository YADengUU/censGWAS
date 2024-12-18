% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_simulation.R
\name{simulate_data}
\alias{simulate_data}
\title{Simulate Data (Phenotype, Covariates, Genotypes) for Testing Purpose}
\usage{
simulate_data(cens_prop, betas, mafs, N, to_truncate)
}
\arguments{
\item{cens_prop}{A non-negative floating-point numeric value indicating the proportion of measurements below LOD.}

\item{betas}{A vector list of floating-point numeric values indicating the effect size of the SNP, with linear association assumed.}

\item{mafs}{A vector list of floating-point numeric values less than 0.5 indicating the minor allele frequency of the SNP.}

\item{N}{An integer indicating the sample size.}

\item{to_truncate}{Boolean indicating whether to have truncated values (NA) or censored. Just for testing the function merge_data().}
}
\value{
Three textfiles - phenotypes.txt, covariates.txt, genotypes.raw - containing the simulated data,
typically for testing step2. In genotypes.raw, a dummy RS-ID of the SNP will be created as the SNP ID.
}
\description{
This function takes the desired proportion below LOD (limit of detection), effect size and MAF of the SNP, and sample sizes
as the inputs and generates corresponding phenotype, covariate, and genotype files for testing the pipelines.
}
