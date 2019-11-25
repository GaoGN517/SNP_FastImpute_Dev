####################################################################################
## Function to compute classification error 
####################################################################################

#' Title Given a confusion matrix, return the classification error.
#'
#' @param confusion_matrix The confusion matrix for prediction summary using the model.
#'
#' @return number of incorrectly classified individuals / total number of individuals
#' @export
#'
#' @examples
#' A <- matrix(c(12, 24, 25, 30))
#' classification_error(A)
#' # 0.8681319
#' 
#' classification_error(diag(c(12, 24)))
#' # 0
#' 
classification_error <- function(confusion_matrix) {
  ##################################################################################
  ## input: the confusion matrix calculated from the test dataset
  ## output: the classification error rate of the prediction using our trained model
  ## on the test dataset. 
  ##################################################################################
  confusion_matrix <- as.matrix(confusion_matrix)
  
  error <- 1 - sum(diag(confusion_matrix)) / sum(confusion_matrix)
  
  return (error)
}