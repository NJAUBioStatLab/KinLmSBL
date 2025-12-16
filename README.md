# KinLmSbl
Kinship-based Linear Model with Sparse Bayesian Learning for GWAS

KinLmSbl is a two-stage genome-wide association study (GWAS) pipeline that integrates a kinship-corrected linear mixed model with sparse Bayesian learning (SBL) for SNP detection. The method effectively controls population structure and polygenic background and is suitable for large-scale GWAS analyses.

## Features
- Kinship correction using REML (FASTmrEMMA framework)
- Population structure adjustment using principal components (Q matrix)
- Fast genome-wide GLM screening via rMVP
- Sparse Bayesian Learning (SBL) for SNP effect estimation
- False Discovery Rate (FDR) correction
- Scalable to large genotype datasets

## Requirements
R >= 4.0.0

Required R packages:
rMVP, data.table, doParallel, Rcpp, RcppArmadillo, RcppEigen

## Required Source Files
The following source files must be placed in the working directory before running the pipeline.

C++ source files:
- sbl3.cpp: Sparse Bayesian Learning implementation
- lmm_diago_fit.cpp: from the gaston package
- kinCorrect.cpp: from FASTmrEMMA (Wen et al., 2018)

R source file:
- emma_REMLE.r: from FASTmrEMMA (Wen et al., 2018)

## Input Files
fileGen: Genotype file (CSV). The first four columns contain SNP metadata (e.g., chromosome, SNP ID, position), and the remaining columns contain the genotype matrix. SNPs are stored by column.
filePhe: Phenotype file (RAW format). Column 6 contains phenotype values and columns 7 onward contain SNP names.
fileKin: Kinship matrix file (CSV), an N × N square matrix. Sample order must be consistent with genotype and phenotype files.
filePS: Population structure file (CSV). The first column contains sample IDs and the remaining columns contain principal components.
fileMAF: SNP minor allele frequency file (.frq).

## Usage
Set the working directory:
setwd("/path/to/KinLmSbl")

Load required packages and source files:
library(rMVP)
library(data.table)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
Rcpp::sourceCpp("sbl3.cpp")
Rcpp::sourceCpp("lmm_diago_fit.cpp")
Rcpp::sourceCpp("kinCorrect.cpp")
source("emma_REMLE.r", encoding = "UTF-8")

Run KinLmSbl:
KinLmSbl(
  dir = getwd(),
  fileGen = "snp.all_plant_height.csv",
  filePhe = "rice01.raw",
  fileKin = "kk.csv",
  filePS = "fixedpc.csv",
  fileMAF = "MAF_check.frq"
)

## Method Description
Population structure is first adjusted by regressing phenotypes on principal components. Variance components are then estimated using REML, and phenotypes and genotypes are corrected using the kinship matrix. A genome-wide GLM scan is performed using rMVP, and SNPs with p < 0.01 are retained. Sparse Bayesian Learning is subsequently applied to the selected SNPs, and final significant SNPs are identified using FDR correction (FDR < 0.05).

## Output
The pipeline generates a single output file named KinLmSbl.result.csv. The output contains the following columns: Chromosome (Chr), SNP ID (SNP), physical position (Position), minor allele frequency (MAF), percentage of phenotypic variance explained (r2(%)), SNP effect (beta), and adjusted p-value (pvalue). If no significant SNPs are detected, an empty file with column headers is generated.

## Notes
Sample order must be consistent across all input files. Large datasets require sufficient memory. Parallel computation is currently disabled (cpu = 1). Genotype matrices are internally transposed for computational efficiency.

## References
Wen Y.J. et al. (2018). FASTmrEMMA: A fast multi-locus GWAS method.
Zhou X. and Stephens M. (2012). Genome-wide efficient mixed-model analysis.
Tipping M.E. (2001). Sparse Bayesian Learning and the Relevance Vector Machine.

## Authors
Wuzhenghan – Sparse Bayesian Learning implementation, 
Caimingzhi
Wu Zhenghan – Pipeline integration and implementation

## License
This project is intended for academic research use. Please cite the relevant references if used in publications.
