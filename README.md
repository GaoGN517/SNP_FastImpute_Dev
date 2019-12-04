
# SNPFastImpute

<!-- badges: start -->
<!-- badges: end -->

The goal of SNPFastImpute is to impute missing values in SNP data files. 

## Installation

You can install the development version of SNPFastImpute from github using devtools:

``` r
devtools::install_github("GaoGN517/689_SNP_FastImpute")
```
(currently private.)

## Basic Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SNPFastImpute)

## Read a vcf file as a matrix
## filename <- "data/Test.vcf" ## your file path
## output_df <- vcf2df(filename) ## convert to a matrix

## Introduce some missing values into the original matrix to test the performance
data("SNP_orig_sub")
## original NA ratio is 0.018

## Make the final missing ratio to be 20%
SNP_NA_df02 <- NA_Generator(SNP_orig_sub, 0.2)

## Make another matrix with final missing ratio to be 5%
SNP_NA_df005 <- NA_Generator(SNP_orig_sub, 0.05)







```

