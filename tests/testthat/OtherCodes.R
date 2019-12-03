
saveRDS(SNP_orig_sub, file="data/SNP_orig_sub.rds")

SNP_df <- subdata
percent <- 0.05
SNP_df_NA_list <- NA_Generator(SNP_df = SNP_df, percent = percent)
SNP_df_NA <- SNP_df_NA_list$SNP_NA_df
introduced_NA_positons <- SNP_df_NA_list$NP_generate_positions

#params <- list(nrounds = 50, booster = "gbtree", objective = "multi:softprob", num_class = 3, eval_metric = "mlogloss")
microbenchmark::microbenchmark(
  predict_SNP_df <- Impute_GenoType_XGBoost(SNP_df_NA, size = 20, nrounds = 50),
  times = 5
)

#df <- SNP_df_NA
NA_orig <- SNP_df[introduced_NA_positons]
NA_pred <- predict_SNP_df[introduced_NA_positons]

classification_error(table(NA_orig, NA_pred))
## 0.03851852, nround = 100, size = 10

# 43 rows, 19999 columns
## show which SNPs contains NAs. 
which(apply(fulldata, 2, function(x) sum(is.na(x)) > 0))
## The first SNP with NA is the 1315th SNP (column).
sum(is.na(fulldata[1:20, 1:500])) # 32
subdata <- fulldata[1:43, 1:1000]
which(is.na(apply(subdata, 2, sum)))
dta_noNA <- t(read.table("testvcf.noNA.txt"))
colnames(dta_noNA) <- dta_noNA[1, ]
dta_noNA <- dta_noNA[-1, ]
dta_noNA[1:5, 1:5]
class(dta_noNA) <- "numeric"
sum(is.na(dta_noNA))


sub_dta_noNA <- dta_noNA[, 1:500]
#class(sub_dta_noNA) <- "numeric"
which(apply(dta_noNA, 2, function(x){
  sum(is.na(x))
})>10, useNames = F)

sum(is.na(sub_dta_noNA))
SNP_df <- sub_dta_noNA
percentage <- 0.2
sub_dta_NA <- NA_Generator(sub_dta_noNA, 0.2)
sum(apply(sub_dta_NA, 2, function(x){
  is.na(x)
}))
