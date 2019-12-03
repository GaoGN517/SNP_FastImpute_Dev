library(usethis)
## Read the full dataset containing with SNPs with NAs and not with NAs
fulldata <- SNPFastImpute::vcf2df("testvcf.vcf")

SNP_orig_sub <- fulldata[, 1:500]
use_data(SNP_orig_sub)

Test_df <- readRDS("Test_df.rds")
use_data(Test_df)


SNP_NA_df02 <- NA_Generator(SNP_orig_sub, 0.2)
use_data(SNP_NA_df02)

SNP_NA_df005 <- NA_Generator(SNP_orig_sub, 0.05)
use_data(SNP_NA_df005)
