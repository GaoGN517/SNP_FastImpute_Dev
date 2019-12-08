SNPFastImpute Package
================
Guannan Gao

  - [SNPFastImpute](#snpfastimpute)
      - [Installation](#installation)
      - [Basic Example](#basic-example)

# SNPFastImpute

<!-- badges: start -->

<!-- badges: end -->

The goal of SNPFastImpute is to impute missing values in SNP data files.

## Installation

You can install the development version of SNPFastImpute from github
using devtools:

``` r
devtools::install_github("GaoGN517/689_SNP_FastImpute")
```

    ## Downloading GitHub repo GaoGN517/689_SNP_FastImpute@master

    ## 
    ##   
       checking for file ‘/private/var/folders/vq/ffxn7ctj7j19fhtjfb26rwj00000gn/T/RtmpCvc8NT/remotesffa12447751e/GaoGN517-689_SNP_FastImpute-1ade028/DESCRIPTION’ ...
      
    ✔  checking for file ‘/private/var/folders/vq/ffxn7ctj7j19fhtjfb26rwj00000gn/T/RtmpCvc8NT/remotesffa12447751e/GaoGN517-689_SNP_FastImpute-1ade028/DESCRIPTION’
    ## 
      
    ─  preparing ‘SNPFastImpute’:
    ## 
      
       checking DESCRIPTION meta-information ...
      
    ✔  checking DESCRIPTION meta-information
    ## 
      
    ─  installing the package to process help pages
    ## 
      
    ─  saving partial Rd database (3.4s)
    ## 
      
    ─  checking for LF line-endings in source and make files and shell scripts
    ## 
      
    ─  checking for empty or unneeded directories
    ## ─  looking to see if a ‘data/datalist’ file should be added
    ## 
      
         NB: this package now depends on R (>= 3.5.0)
    ## 
      
         WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects:  'SNPFastImpute/vignettes/SNP_missing_value_imputation_process_cache/html/predict  WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects: for different missing ratio before  and after considering  WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects:  correlation_c73172829437dfac00fa71a5e27d4945.RData'  WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects:  'SNPFastImpute/vignettes/SNP_missing_value_imputation_process_cache/html/predict  WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects: for different missing ratio before  and after considering  WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects:  correlation_c73172829437dfac00fa71a5e27d4945.rdx'
    ## 
      
    ─  building 'SNPFastImpute_1.0.0.9000.tar.gz'
    ## 
      
       Warning in utils::tar(filepath, pkgname, compression = compression, compression_level = 9L,  :
    ##      storing paths of more than 100 bytes is not portable:
    ##      'SNPFastImpute/vignettes/SNP_missing_value_imputation_process_cache/html/predict for different missing ratio before and after considering correlation_c73172829437dfac00fa71a5e27d4945.RData'
    ## 
      
       Warning in utils::tar(filepath, pkgname, compression = compression, compression_level = 9L,  :
    ##      storing paths of more than 100 bytes is not portable:
    ##      'SNPFastImpute/vignettes/SNP_missing_value_imputation_process_cache/html/predict for different missing ratio before and after considering correlation_c73172829437dfac00fa71a5e27d4945.rdb'
    ## 
      
       Warning in utils::tar(filepath, pkgname, compression = compression, compression_level = 9L,  :
    ##      storing paths of more than 100 bytes is not portable:
    ##      'SNPFastImpute/vignettes/SNP_missing_value_imputation_process_cache/html/predict for different missing ratio before and after considering correlation_c73172829437dfac00fa71a5e27d4945.rdx'
    ## 
      
       Warning in utils::tar(filepath, pkgname, compression = compression, compression_level = 9L,  :
    ##      using GNU extension for long pathname
    ## 
      
       
    ## 

I am using a windows to develop this package and I used another windows
system to test the installation. I realized very late that this package
may met some specific problem when installing on mac.  
As I currently not using C++ features, I decided to remove C related
parts to make the package more easy to run.

## Basic Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SNPFastImpute)

## Read a vcf file as a matrix
## filename <- "data/Test.vcf" ## your file path
## vcf_df <- read.table(filename) ## read in your vcf data 
## Here we just load the dataset in vcf format.
data(vcf_df)
output_df <- vcf2df(vcf_df)
## There would be warning message for NAs, which is caused by adding missing positions
## of SNPs, which is what we should do. So here NAs does not mean problem. 

## Introduce some missing values into the original matrix to test the performance
data("SNP_orig_sub")
## this is the SNP matrix with the original SNP types.
## original NA ratio is 0.018

## Make the final missing ratio to be 20%
SNP_NA_df02 <- NA_Generator(SNP_orig_sub, 0.2)
## This is an object list of four elements, 
## SNP_NA_df, 
## NA_percent_orig,
## NA_percent_generate,
## NP_generate_positions

## Make another object with final missing ratio to be 5%
SNP_NA_df005 <- NA_Generator(SNP_orig_sub, 0.05)

## Now we have the original matrix, the two objects with additional missing values.
ls()
# [1] "SNP_NA_df005" "SNP_NA_df02"  "SNP_orig_sub"

## We can use the following function to create an object that can be used in the imputation 
## function for each SNP in the matrix.
Create_Single_SNP_Object(SNP_NA_df005$SNP_NA_df, 2, size = 20)

## The main function of this package is to use the imputation function to fill in the missing
## values.

## We can just do the filling for the original matrix
system.time(
  predict_df <- Impute_GenoType_XGBoost(SNP_orig_sub, size = 10)
)
##   user  system elapsed 
##   1.890   0.037   1.957 
## If imputation time is long for large files, this is also an paralleled version of this function:
## Impute_GenoType_XGBoost. 

## We can also perform the filling on the matrix which we introduced additional 
## missing values. And then see how our method performed on predictions. 
system.time(
  df_fill02 <- Impute_GenoType_XGBoost(SNP_NA_df02$SNP_NA_df)
)
##   user  system elapsed 
## 11.455   0.265  12.257 
 
system.time(
  df_fill005 <- Impute_GenoType_XGBoost(SNP_NA_df005$SNP_NA_df)
)
##    user  system elapsed 
##   10.245   0.142  10.804 

NA_positions02 <- SNP_NA_df02$NP_generate_positions
NA_positions005 <- SNP_NA_df005$NP_generate_positions

classification_error(SNP_orig_sub, df_fill02, NA_positions02) ## 0.071
classification_error(SNP_orig_sub, df_fill005, NA_positions005) ## 0.061
```
