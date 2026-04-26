# KinLmSBL

A two-stage GWAS algorithm that integrates single marker scanning and sparse Bayesian learning to improve the performance of detection.

## Description

It is a two-stage GWAS algorithm. The first stage of KinLmSBL involves constructing a single-locus mixed linear model and strategically implementing a matrix transformation for polygenic background correction. Next, general linear regression is performed to scan high-throughput markers, identifying a small number of potentially associated SNPs for the second stage. Finally, a multi-locus model is constructed, and SBL is performed to identify the significant loci with FDR correction.

### Two-Stage Analysis:

1. **First Stage**: It constructs a mixed linear model and corrects for the polygenic background using FASTmrEMMA's matrix transformation. It then employs rMVP's MVP.GLM function to scan high-throughput SNPs, screening for a small number of potentially associated SNPs for the next stage. 

2. **Second Stage**: It constructs a hierarchical linear model, applies the SBL algorithm to estimate marker effects and construct the Wald test statistic, and implements FDR correction. 

Eventually, significant loci related to the trait are identified.

## Installation

### System Requirements

- R >= 3.5.0
- C++11 compiler
- Required R packages: KinLmSBL_1.0.0.tar.gz, rMVP, data.table, Matrix, doParallel, Rcpp, RcppArmadillo, RcppEigen

### From source:

```
# Install dependencies first
install.packages(c("KinLmSBL_1.0.0.tar.gz", "rMVP", "data.table", "Matrix", "doParallel", 
                   "Rcpp", "RcppArmadillo", "RcppEigen"))

```

# load the the KinLmSBL core function
Rcpp::sourceCpp("sbl.cpp")
Rcpp::sourceCpp("kinCorrect.cpp")
Rcpp::sourceCpp("lmm_diago_fit.cpp")
source("emma_REMLE.r", encoding = "UTF-8") 

## Usage

```r
library(KinLmSBL)

# Run GWAS analysis
result <- KinLmSBL(
  Y = "phenotype.csv",        # Phenotype file
  GD = "genotype.012.csv",    # Genotype file (012 format)
  KI = "kinship.csv",         # Kinship matrix
  CV = "covariates.csv",      # Covariate file (e.g., PCs)
  file.output = TRUE,         # Save results to file
  dir = getwd()               # Output directory
)

# View results
print(result)
```

## Input Data Format

### Phenotype Data
Phenotype data should be provided in CSV format with a header. The second column contains phenotypic values. Missing values should be coded as NA or NaN. Any missing values in the dataset should be preprocessed by the user beforehand.

### Genotype Data
Genotype data should be provided in CSV format using 0/1/2 coding. The first four columns contain marker information.

### Kinship matrix
The kinship matrix should be a square matrix. Row and column names must correspond to individual IDs in the phenotype data.

### Population Structure
Population structure (covariates) should be provided as a matrix with individuals in rows and covariates in columns. The order of individuals must match the phenotype file.

## Output

Returns a data frame with significant SNPs containing:

- **SNP**: SNP ID
- **Chr**: Chromosome 
- **Pos**: Position on chromosome
- **P.value**: FDR-adjusted P-value
- **MAF**: Minor allele frequency
- **r2(%)**: Variance explained by the SNP
- **Effect**: Effect size

When `file.output = TRUE`, results are saved to `KinLmSBL.result.csv`.

## Methods

This package uses efficient C++ implementations from:

- **sblgwas**: from Wang et al.(2019) and Wu Zhenghan
- **kinCorrect**: from FASTmrEMMA (Wen et al.,2018)
- **lmm_diago_fit**: from Package gaston
- **emma_REMLE**: from FASTmrEMMA (Wen et al.,2018) and Package gaston

## References

- Zhang, J., Wu, Z., Cai, M., Liu, K., Han, X., Liu, C., Han, G., & Wen, Y. (2026). Integrated single marker scanning and sparse Bayesian learning improves performance of detection for GWAS. Plant methods, 22(1), 40. https://doi.org/10.1186/s13007-026-01516-7

- Wang, M., & Xu, S. (2019). A coordinate descent approach for sparse Bayesian learning in high dimensional QTL mapping and genome-wide association studies. Bioinformatics (Oxford, England), 35(21), 4327–4335. https://doi.org/10.1093/bioinformatics/btz244

- Wen, Y. J., Zhang, H., Ni, Y. L., Huang, B., Zhang, J., Feng, J. Y., Wang, S. B., Dunwell, J. M., Zhang, Y. M., & Wu, R. (2018). Methodological implementation of mixed linear models in multi-locus genome-wide association studies. Briefings in bioinformatics, 19(4), 700–712. https://doi.org/10.1093/bib/bbw145

## License

## Author

Designed by Yangjun Wen
Written by Zhenghan Wu, Gang Han, Mingzhi Cai, and Jin Zhang
