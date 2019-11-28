#' Impute missing SNPs for the input dataset. 
#'
#' @param XGBoost_dataset The XGBoost data type dataset including the missing SNPs
#' @param nrounds number of rounds in xgboost model training. 
#' @param booster
#' @param objective
#' @param num_class How many classes to be predicted (in this case, equal to 3).
#' @param eval_metric
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
#' df <- readRDS("Test_df.rds")
#' predict_df <- Impute_GenoType_XGBoost(df)
#' 
#' size = 10
#' nrounds = 100
#' num_class = 3
#' booster = "gbtree"
#' objective = "multi:softprob"
#' eval_metric = "mlogloss"
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
Impute_GenoType_XGBoost <- function(df, df_test = NULL, size = 10, nrounds = 100, booster = "gbtree", objective = "multi:softprob", num_class = 3, eval_metric = "mlogloss"){
  ## To save time, we only impute SNP columns with NAs in the dataset.
  ## Get the dimension of df
  n <- nrow(df)
  p <- ncol(df)
  ## check whether there are enough 
  if (size > p) stop("The size of the windows should be smaller than the number of SNPs.")
  ## Read the parameters for model building.
  params <- list(booster = booster, objective = objective, num_class = num_class, eval_metric = eval_metric)
  
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
    
    for(i in 1:n_NA_cols) {
      single_SNP_Obj <- Create_Single_SNP_Object(df = df, a = i, size = size)
      train_xgboost <- xgboost::xgb.DMatrix(data = single_SNP_Obj$train_data, label = single_SNP_Obj$train_label)
      pred_xgboost <- xgboost::xgb.DMatrix(data = single_SNP_Obj$pred_data, label = single_SNP_Obj$pred_label)
      xgb_model <- xgboost::xgb.train(data = train_xgboost, num_class = num_class, nrounds = nrounds)
      single_SNP_Obj$pred_label <- predict(xgb_model, newdata = single_SNP_Obj$pred_data)
      #error_vec[i] <- single_SNP_Obj$error
      ## Filling the the NA positions
      df_fill[single_SNP_Obj$NA_positions, single_SNP_Obj$SNP_position] <- single_SNP_Obj$pred_label
    }
  }
  ## Return the column index with missing values, and the corresponding SNP_Objects for these SNPs.
  #return(list(df_fill = df_fill, error = error_vec))
  return(df_fill)
}


Impute_GenoType_XGBoost1 <- function(df, size = 10, nrounds = 100, booster = "gbtree", objective = "multi:softprob", num_class = 3, eval_metric = "mlogloss") {
  ## Get how many SNPs includes missing values needing to be imputed. 
  XGBoost_dataset <- Create_SNP_XGBoost_Object(df, size = size)
  test_length = length(XGBoost_dataset$NA_cols)
  
  ## insert the parameters for 
  params <- list(booster = booster, objective = objective, num_class = num_class, eval_metric = eval_metric)
  
  multi_SNP_Obj <- XGBoost_dataset$Multiple_SNP_Object
  ## Model and prediction
  for (i in 1: test_length) {
    single_SNP_Obj <- multi_SNP_Obj[[i]]
    if(single_SNP_Obj$model_fit == T) {
      train_data <- single_SNP_Obj$train_xgboost
      xgb_model <- xgboost::xgb.train(data = train_data, params = params, num_class = 3)
      
    }
    
  }
  XGBoost_dataset_sub <- lapply(XGBoost_dataset[[1]][c(1:test_length)], function(x) {
    xgb_model <- xgboost::xgb.train(params = params, data =x[[1]], nrounds = 100)
    # Predict for validation set
    xgb_val_tests <- predict(xgb_model, newdata = x[[2]])
    xgb_val_out <- matrix(xgb_val_tests, nrow = 3, ncol = length(xgb_val_tests) / 3) %>% 
      t() %>%
      data.frame() %>%
      mutate(max = max.col(., ties.method = "last"), label = x[[3]]) 
    
    x[[5]] <- xgb_val_out$max
    
    # Confustion Matrix
    xgb_val_conf <- table(true = x[[3]], pred = x[[5]])
    x[[6]] <- classification_error(xgb_val_conf)
    #df_sub[x[[4]], i] <- x[[5]]-1
    
    if (x[[6]] > 0) {
      cat("XGB Validation Classification Error Rate of", x[[8]], " is: ", x[[6]], "\n")
      cat("Test should be: ",x[[3]],"\n")
      cat("But turns out to be: ", x[[5]]-1,"\n")
    }
    else {
      cat("Passed the test. \n")
      cat("XGB Validation Classification Error Rate of", x[[8]], " is: ", x[[6]], "\n")
      cat("Test should be: ",x[[3]],"\n")
      cat("But turns out to be: ", x[[5]]-1,"\n")
      xgb_preds <- predict(xgb_model, newdata = x[[7]])
      xgb_preds_out <- matrix(xgb_preds, nrow = 3, ncol = length(xgb_preds) / 3) %>% 
        t() %>%
        data.frame() %>%
        mutate(max = max.col(., ties.method = "last")) 
      df_sub[x[[4]], x[[8]]] <- xgb_preds_out$max-1
      cat("Prediction is ",df_sub[x[[4]], x[[8]]],"\n")
    }
    cat("\n")
    return(x)
  })
  
  ## Show the results: 
  bind <- cbind(df_sub, XGBoost_dataset[[3]])
  mean_error <- mean(as.numeric(unlist(lapply(XGBoost_dataset_sub, function(x) x[[6]]))))
  cat(paste0("The mean classification error is ", mean_error, ".\n"))
  
  return(df_sub)
}