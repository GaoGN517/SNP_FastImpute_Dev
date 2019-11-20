####################################################################################
## This is a Test File for the dataset from Yue Xing, with 20% of missing data 
## randomly generated for each SNP. 
####################################################################################
source("Codes/FunctionsSNP.R")

############## main ######################################################################
start_time <- Sys.time()
chr1A <- vcf2gt(filename = "../Data/chr1A.vcf")
check1_time <- Sys.time() 

cat("Time to convert vcf file is ", check1_time - start_time, "\n")

#chr1A <- read.table("chr1A_gt.txt")
XGBoost_datasets <- Generate_XGBoost_Dataset(chr1A, 100)
#filename <- "chr1A.vcf"
cat("Time to generate xgboost dataset is ", Sys.time() - check1_time, "\n")

check2_time = Sys.time()
#cl <- makeCluster(detectCores())
out <- Impute_GenoType_XGBoost(XGBoost_datasets)
#stopCluster(cl)
end_time <- Sys.time()
cat("Time to build model and finish prediction is ", end_time - check2_time, "\n")

cat(paste0("The total runing time is ", end_time - start_time))

########## debug ###########################################################################
corr <- cor(chr1A, use = "pairwise")
corr[1:5, 1:5]
#corr1 <- cor(chr1A, use = "complete.obs")
#corr1[1:5, 1:5]
#filename = "chr1A.vcf"
#SNP <- 
#grep("chr1A_part1:2660229", SNP) # 90,129,181
#SNP[182] # "chr1A_part1:2660258","chr1A_part1:3383087","chr1A_part1:3582726"
#x <- XGBoost_dataset[[1]][[91]]

#x <- XGBoost_dataset[[1]]$'chr1A_part1:3383087'
# xgb.DMatrix  dim: 68 x 100  info: label  colnames: yes
#XGBoost_datasets[[1]]$'chr1A_part1:2660258'$test_df_label