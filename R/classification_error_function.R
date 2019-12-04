####################################################################################
## Function to compute classification error 
####################################################################################

#' Classification error calculator.
#'
#' @param df1 The original dataset that contains only true missing values, not introduced 
#'         missing values. 
#' @param df2 The second dataset that has all missing values filled in.
#' @param NA_positions The positions where we introduced NA in the original dataset.
#'
#' @return number of incorrectly classified individuals / total number of individuals
#' @export
#'
#' @examples
#' data("SNP_orig_sub")
#' data("SNP_NA_df02")
#' df_fill <- Impute_GenoType_XGBoost(SNP_NA_df02$SNP_NA_df)
#' NA_positions <- SNP_NA_df02$NP_generate_positions
#' classification_error(SNP_orig_sub, df_fill, NA_positions)
#' 
classification_error <- function(df1, df2, NA_positions) {
  ## Calculate the confusion matrix
  confusion_matrix <- as.matrix(table(df1[NA_positions], df2[NA_positions]))
  
  ## Calculate the classification error. 
  error <- 1 - sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  return (error)
}


