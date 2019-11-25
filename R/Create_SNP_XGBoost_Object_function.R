#' Input a dataframe, generate an object for the SNP missing value imputation with xgboost.
#'
#' @param df A dataframe of the SNPs which we need to impute missing values
#' @param size A windows size which we control for our imputation. Default value is 50. 
#'
#' @return a list of list of 3 contents:  
#' 1. XGBoost_dataset  
#' 1.1. xgb_train: X_train  
#' 1.2. xgb_test: X_test  
#' 1.3. test_df_label  
#' 1.4. snp_NA  
#' 1.5. pred: Empty now, to store the predicted values for each SNP  
#' 1.6. error: Empty now, to store the prediction errors for each SNP  
#' 1.7. xgb_pred  
#' 1.8. SNP  
#'   
#' 2. df_NA  
#' 3. df_noNA  
#' 
#' @export
#'
#' @examples
#' df <- readRDS("Test_df.rds")
#' xgboost_snp_obj <- Create_SNP_XGBoost_Object(df)
#' 
#' xgboost_snp_obj <- Create_SNP_XGBoost_Object(df, size = 100)
#' 
Create_SNP_XGBoost_Object <- function(df, size = 50) {
  ## To save time, we only impute SNP columns with NAs in the dataset.
  NA_cols <- which(apply(df, 2, function(x) sum(is.na(x))> 0))
  df_sub <- df[,NA_cols]
  df_rest <- df[, -NA_cols]
  N <- ncol(df_sub)
  
  #df_noNA_sub <- df_noNA[, NA_cols]
  set.seed(101)
  XGBoost_dataset <- list()
  for (i in 1:N) {
    XGBoost_dataset[[i]] <- list()
    snp_NA <- which(is.na(df_sub[,i]))
    snp_notNA <- which(!is.na(df_sub[,i]))
    
    samp <- sample(length(snp_notNA), size = length(snp_notNA)*0.2)
    train <- snp_notNA[-samp]
    test <- snp_notNA[samp]
    if (sum(snp_NA) != 0) { ## double check that this column has NAs to impute
      ## Create train and test data
      if (i == 1)  range = c((i+1): size)
      else if (i == N)  range = c((i-size): N)
      else if (i <= size/2) range = c(1:(i-1), (i+1): size)
      else if ((i + size/2) >= N)  range = c((i-size):(i-1), (i+1): N)
      else range = c((i-size/2):(i-1), (i+1): (i + size/2))
      
      train_df <- as.matrix(df_sub[train, range])
      test_df <- as.matrix(df_sub[test, range])
      pred_df <- as.matrix(df_sub[snp_NA, range])
      #pred_df <- as.matrix(df_sub[, range])
      
      if(ncol(test_df)==1) test_df <- t(test_df)
      if(ncol(pred_df)==1) pred_df <- t(pred_df)
      
      ## Create labels
      train_df_label <- as.numeric(df_sub[train, i])
      test_df_label <- as.numeric(df_sub[test, i])
      pred_df_label <- as.numeric(df_sub[snp_NA, i]) # change from df_noNA_sub to df_sub
      
      ## Prepare matrices
      XGBoost_dataset[[i]][[1]] <- xgboost::xgb.DMatrix(data = train_df, label = train_df_label)
      XGBoost_dataset[[i]][[2]] <- xgboost::xgb.DMatrix(data = test_df, label = test_df_label)
      XGBoost_dataset[[i]][[3]] <- test_df_label
      XGBoost_dataset[[i]][[4]] <- snp_NA
      XGBoost_dataset[[i]][[5]] <- list()
      XGBoost_dataset[[i]][[6]] <- list()
      XGBoost_dataset[[i]][[7]] <- xgboost::xgb.DMatrix(data = pred_df, label = pred_df_label)
      XGBoost_dataset[[i]][[8]] <- colnames(df_sub)[i]
      
      names(XGBoost_dataset[[i]]) <- c("xgb_train", "xgb_test", "test_df_label", "snp_NA","pred","error", "xgb_pred","SNP")
    }
  }
  names(XGBoost_dataset) <- colnames(df_sub)
  out <- list()
  out[[1]] <- XGBoost_dataset
  out[[2]] <- df_sub
  out[[3]] <- df_rest
  names(out) <- c("XGBoost_dataset","df_NA","df_noNA")
  return(out)
}