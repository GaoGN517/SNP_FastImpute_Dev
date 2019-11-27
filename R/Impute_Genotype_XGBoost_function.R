#' Impute missing SNPs for the input dataset. 
#'
#' @param XGBoost_dataset The XGBoost data type dataset including the missing SNPs
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
#' SNP_XGBoost_Obj <- readRDS("Test_xgboost_obj.rds")
#' predict_df <- Impute_GenoType_XGBoost(SNP_XGBoost_Obj)
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
Impute_GenoType_XGBoost <- function(XGBoost_dataset) {
  ####################################################################################
  ## In the next step I should change the size selection to this function. 
  ####################################################################################
  ## Greb the number of 
  test_length = length(XGBoost_dataset[[1]])
  params <- list(booster = "gbtree", objective = "multi:softprob", num_class = 3, eval_metric = "mlogloss")
  df_sub <- XGBoost_dataset[[2]]
  ## Model and prediction
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