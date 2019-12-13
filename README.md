
# SNPFastImpute

<!-- badges: start -->
<!-- badges: end -->

The goal of SNPFastImpute is to impute missing values in SNP data files. 

## 1. Installation

You can install the development version of SNPFastImpute from github using devtools:


``` {r installation, include = F}
devtools::install_github("GaoGN517/689_SNP_FastImpute")
```

## 2. Basic Usage

Load the package:
```{r attatch, results='hide'}
library(SNPFastImpute)
```

### 2.1. Read a vcf file as a matrix
For example, if your data file name is "data/Test.vcf". Then the following codes can read in the vcf file.
```
filename <- "data/Test.vcf" # your file path
vcf_df <- read.table(filename) # read in your vcf data 
```
Here we just load the dataset in vcf format.

```{r, echo=TRUE, results='hide', warning=FALSE}
data(vcf_df)
output_df <- vcf2df(vcf_df)
```
There would be warning message for NAs, which is caused by adding missing positions of SNPs, which is what we should do. So here NAs do not mean problem.

### 2.2. Introduce some missing values into the original matrix to test the performance
```{r}
data("SNP_orig_sub")
```
This is the SNP matrix with the original SNP types with original NA ratio 0.018.

Make the final missing ratio to be 20%
```{r}
SNP_NA_df02 <- NA_Generator(SNP_orig_sub, 0.2)
```
This is an object list of four elements:
* SNP_NA_df
* NA_percent_orig
* NA_percent_generate
* NP_generate_positions

Make another object with final missing ratio to be 5%
```{r}
SNP_NA_df005 <- NA_Generator(SNP_orig_sub, 0.05)
```

Now we have the original matrix, the two objects with additional missing values.
We can use the following function to create an object that can be used in the imputation function for each SNP in the matrix.
### 2.3. Create SNP_Object
This object contains all the information we need for the imputation of the target single SNP column with missing values. 
```{r}
Create_Single_SNP_Object(SNP_NA_df005$SNP_NA_df, 2, size = 20)
```

### 2.4. SNP missing value imputation:
We can just do the filling for the original matrix
```{r}
predict_df <- Impute_GenoType_XGBoost(SNP_orig_sub, size = 10)
```


If imputation time is long for large files, this is also an paralleled version of this function:
```{r}
predict_df_para <- Impute_GenoType_XGBoost_para(SNP_orig_sub, size = 10)
```

### 2.5. Calculation of classification error
We can perform the filling on the matrix which we introduced additional missing values. And then see how our method performed on predictions. 
```{r}
## Imputation
df_fill02 <- Impute_GenoType_XGBoost(SNP_NA_df02$SNP_NA_df)
df_fill005 <- Impute_GenoType_XGBoost(SNP_NA_df005$SNP_NA_df)
## Find introduced NA positions
NA_positions02 <- SNP_NA_df02$NP_generate_positions
NA_positions005 <- SNP_NA_df005$NP_generate_positions
## Calcluate classification error
classification_error(SNP_orig_sub, df_fill02, NA_positions02) ## 0.071
classification_error(SNP_orig_sub, df_fill005, NA_positions005) ## 0.061
```



