#' Impute missing SNPs for the input dataset. 
#'
#' @param df The original dataset including the missing SNPs to be imputed.    
#' @param size The windows size to use as the training dataset for each SNP, default: 10.   
#' @param nrounds Number of fitting rounds, default: 100. 
#' @param num_class Number of classes of response variable (types of SNPs), default: 3.  
#' 
#' params A list of parameters for the xgboost model building. 
#' Default: nrounds = 100, booster = "gbtree", objective = "multi:softprob", 
#' num_class = 3, eval_metric = "mlogloss".   
#'    
#' 
#' @details 
#' In our model, we try to use the types of SNPs around each missing SNP to predict the missing value. 
#' For each missing value, we need to use the size n of SNPs around it as predictors,
#' and use the non-missing samples for this SNP position as the training dataset. 
#' 
#' @return The predicted missing genotypes.
#' @export
#'
#' @examples
#' data("Test_df")
#' predict_df <- Impute_GenoType_XGBoost(Test_df, size = 10)
#' ## May take several seconds to finish.
#' ## Should return a dataset where the missing values are filled by predicted values.  
#' 
#' @references 
#' 
#' Kabisch, Maria, Ute Hamann, and Justo Lorenzo Bermejo. "Imputation of 
#' missing genotypes within LD-blocks relying on the basic coalescent and 
#' beyond: consideration of population growth and structure." BMC genomics 18.1 (2017): 798.
#' 
#' Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System",
#' 22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, \url{https://arxiv.org/abs/1603.02754}
#' 
Impute_GenoType_XGBoost <- function(df, size = 10, num_class = 3, nrounds = 100){
  ## To save time, we only impute SNP columns with NAs in the dataset.
  ## Get the dimension of df
  n <- nrow(df)
  p <- ncol(df)
  ## check whether there are enough 
  if (size > p) stop("The size of the windows should be smaller than the number of SNPs.")
  ## Read the parameters for model building.

  NA_cols <- which(apply(df, 2, function(x) sum(is.na(x))> 0))
  n_NA_cols <- length(NA_cols)
  ## Check whether there are missing values in the dataset.
  if (n_NA_cols == 0) {
    print("There is no missing value in the data set.")
  }
  else {
    ## Create an object to store the dataframe where missing SNPs are replaced by predicted SNPs.
    df_fill <- df
    error_vec <- rep(NA, n_NA_cols)
    for(col in NA_cols) {
      single_SNP_Obj <- Create_Single_SNP_Object(df = df, a = col, size = size)
      train_xgboost <- xgboost::xgb.DMatrix(data = single_SNP_Obj$train_data, label = single_SNP_Obj$train_label)
      pred_xgboost <- xgboost::xgb.DMatrix(data = single_SNP_Obj$pred_data, label = single_SNP_Obj$pred_label)
      xgb_model <- xgboost::xgb.train(data = train_xgboost, num_class = num_class, nrounds = nrounds)
      single_SNP_Obj$pred_label <- predict(xgb_model, newdata = single_SNP_Obj$pred_data)
      df_fill[single_SNP_Obj$NA_positions, single_SNP_Obj$SNP_position] <- single_SNP_Obj$pred_label
    }
  }
  ## Return the column index with missing values, and the corresponding SNP_Objects for these SNPs.
  #return(list(df_fill = df_fill, error = error_vec))
  return(df_fill)
}
