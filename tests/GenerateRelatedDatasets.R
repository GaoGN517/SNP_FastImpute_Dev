## Read the full dataset containing with SNPs with NAs and not with NAs
fulldata <- SNPFastImpute::vcf2df("testvcf.vcf")

SNP_orig_sub <- fulldata[, 1:500]
use_data(SNP_orig_sub)

