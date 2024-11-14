# censGWAS: Perform GWAS on Phenotype Measurements Subject to Limit of Detection for Biobank-Scale Dataset
`censGWAS` is an R package for conducting genome-wide association study (GWAS) analyses on phenotypes subject to the limit of detection (LOD), which makes the measurements left-truncated or censored. Taking detection power, accuracy and efficiency into account, we have proposed the Linear-Tobit scheme. This package runs the first step, GWAS with linear model, by PLINK2 ([Chang et al 2015](https://www.cog-genomics.org/plink/2.0)), using which it also clumps the results to select the genetic variants, i.e., single nucleotide polymorphisms (SNPs), for further refinement to mitigate the impacts of biased distribution by the Tobit model implemented through `censReg` ([Henningsen 2024](https://cran.r-project.org/web/packages/censReg/censReg.pdf)).

While the pipeline can be run all at once, you may choose to start from clumping by providing your own GWAS results, or, just the refinement by providing the list of interested SNPs, since we are aware of the various choices of GWAS tools for the association testing and clumping, see instruction details below.

## Installation Requirements
To use this package, the following dependencies are required, including the `devtools` for package installation, additional R packages for data processing and analyses, and the software PLINK2.

### 1. R package dependencies
Ensure you have the following R packages installed:
```r
# Install devtools (if you haven't) for installing this package from GitHub
install.packages("devtools")

# install dependencies for data processing and analyses
install.packages("RNOmni", "censReg", "broom")
```
### 2. Installing PLINK2
PLINK2 is necessary for processing genetic data. Even if you choose to run GWAS and clumping with external tools, PLINK2 is needed to extract and convert genetic data into readable formats by R in the refinement step. You can download from the [official website of PLINK2](https://www.cog-genomics.org/plink/2.0/) and follow the instructions for your operating system. If you are on UPPMAX or similar HPC servers, you may install by:
```bash
module load bioinfo-tools plink2
```
which by default retrieves the latest version.

### 3. Installing censGWAS
After the prerequisites are fulfilled, you can install this package directly from GitHub in R using `devtools`:
```r
devtools::install_github("YADengUU/censGWAS")
```

## Getting started
Once the package is installed, you can use the core functions to perform your analyses. Here is a quick example that performs the whole pipeline - GWAS with linear model, clump the hits, and generate the refined Tobit estimates:
```r
library(censGWAS)

censGWAS(pfile = "chr1",
         file_pheno = "phenotypes.txt",
         file_covar = "covariates.txt",
         participants = "participants.txt")
```
in which: 

- `pfile` stands for the filename of three genotype files (genotype calls - .pgen, variant information - .pvar, sample information - .psam); in case they do not share the same name, you can specify separately with the parameters `pgen`, `pvar`, and `psam`, e.g., `pgen = "chr1.pgen"`, `pvar = "chr1.pvar"`, `psam = "ukb_project.psam"`. Similarly, if your data are in PLINK1 format instead, use `bfile`, or, `bed = "chr1.bed"`, `bim = "chr1.bim"`, `fam = "ukb_project.fam"`.

- `file_pheno` is the file including the phenotype measurements.

- `file_covariates` is the file with covariate data.

- `participants` is the file with the list of participants for the analyses.

Formatting details for the input files are described in later sections.

By default, the primitive GWAS is conducted on SNPs that have Hardy-Weignberg equilibrium exact p-value above 1e-20, minor allele frequencies (MAF) above 0.01, and missing call rates below 0.05; for each effect estimate it reports the 95% confidence intervals. These can be customized with the parameters `gwas_thres_hwe`, `gwas_thres_maf`, `gwas_thres_geno`, and `gwas_ci`. For example, you can add `gwas_thres_maf = 0.001`, `gwas_ci = 0.9`.

Sometime PLINK2 terminates with the warning of covariate scales varying to widely for numerical stability and requests to standardize the numerical covariates. You may choose to include this step by including the flag `standardize_covariates = TRUE` when using censGWAS().

### Already performed GWAS, need clumping and refinement
If GWAS was already performed, you can directly provide the GWAS result file by adding the parameter `gwas_results`:
```r
library(censGWAS)

censGWAS(pfile = "chr1",
         file_pheno = "phenotypes.txt",
         file_covar = "covariates.txt",
         participants = "participants.txt",
         gwas_results = "chr1.tf_lin.glm.linear")
```
When the package detects that `gwas_results` is provided, the pipeline will start from clumping the results and refine the estimates for the selected SNPs.

Depending on the software, the key information in GWAS results, including the SNP ID and position, the effect allele and its frequency, the effect estimate, standard error, and p-value, may have column/field names differ from the default (as in PLINK2), "ID", "POS", "A1", "A1_FREQ", "BETA", "SE", and "P", respectively. You can customize the names in the parameters `snp_id_field`, `pos_field`, `a1_field`, `a1_feq_field`, `beta_field`, `se_field`, and `p_field`, e.g., `snp_id_field = "rsid"`, `a1_field = "Allele2"`, and `p_field = "p.value"`, etc.

For clumping, the following default values are used for the index variant p-value threshold $p1 = 5e-8$, the SP2 column p-value threshold $p2 = 1e-4$, the threshold of correlation between genetic variants $r^2=0.1$, the maximum distance from the lead variant for SNPs $kb=1000$. These could be customized with `clump_p1`, `clump_p2`, `clump_r2`, and `clump_kb`: for example, adding the flag `clump_kb=250` to `censGWAS()` in the usage above.

### Have specific SNPs chosen, only need refinement
If you have also conducted your own clumping, or simply have specific SNPs of interest, you can organize them into one plaintext file and provide it through the parameter `file_snps_chosen`. In this case, you do not need to provide the GWAS results (but genotype files and phenotype/covariates are still required), and the package will start from the refinement procedure.
```r
library(censGWAS)

censGWAS(pfile = "chr1",
         file_pheno = "phenotypes.txt",
         file_covar = "covariates.txt",
         participants = "participants.txt",
         file_snps_chosen = "chr1_snps.chosen")
```
The final result of the package is the refined estimates (along with the standard error and p-values) organized in one file, by default the output will be named "out.refine". You can customize the output name for better clarity, e.g., `output_name = "groupA_chr1"`, the suffix will be added automatically.

In (perhaps rare) case that the memory space is scarce and you realize that the extracted but unconverted genotype files occupy unwanted space, you may choose to delete them as the subsequent refinement calculations won't need the unconverted files by having the flag `delete_unconverted = TRUE`.

## Input file formatting
To ensure the package runs smoothly, the input files should be in space- or tab-delimited formats as described below.

### Phenotype file
The phenotype data should be organized in 4 columns: FID and IID of participants, censored status (1 if above LOD, 0 otherwise), and the measurements. The package detects whether the measurements are truncated (`NA` for those below LOD) or censored (already assigned those below LOD with the LOD value); if they are truncated, the package modifies them to be censored.

- Truncated
```plaintext
FID IID cens_status measurements
123 101 1 3.8
121 100 1 3.9
119 103 0 NA
```
- Censored
```plaintext
FID IID cens_status measurements
123 101 1 3.8
121 100 1 3.9
119 103 0 3.0
```

In the preprocessing step, the package performs inverse normal rank-based transformation for the censored values, implemented by `RNOmni` ([McCaw 2022](https://cran.r-project.org/web/packages/RNOmni/RNOmni.pdf)). The modified phenotype file will be stored in your output directory.

### Covariate file
Like the phenotypes, the covariates are also organized as columns after the first two columns with the FID and IID. For example:
```plaintext
FID IID age smoking
123 101 39 never
121 100 32 occasional
119 103 46 previous
```

### Participant list
This is just the file with two columns - FID and IID - of the participants selected for your study. The package checks whether they are included in the phenotype/covariate files and generated the updated list after preprocessing step in which only the ones with complete cases will be forwarded.
```plaintext
FID IID
123 101
121 100
119 103
```

### SNPs chosen
In the case that you already have chosen the SNPs for refining the estimates, you can organize them in three columns including the chromosome number, SNP ID, and the effect allele:
```plaintext
chr ID  A1
1 rs123 G
1 rs101 A
1 rs111 T
```
The package checks if your list is empty, if it is not correctly formulized in three columns, if the chromosome numbers are correct (1~23, X/XY/Y, case insensitive), and if your nucleotide bases are correct (A, T, G, C) before extracting the desired genotype information.
## References
- Henningsen A (2024). censReg: Censored Regression (Tobit) Models. R package version 0.5-39, https://r-forge.r-project.org/projects/sampleselection.
- Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4.
- McCaw Z (2024). RNOmni: Rank Normal Transformation Omnibus Test. R package version 1.0.1.3, https://github.com/zrmacc/rnomni.

