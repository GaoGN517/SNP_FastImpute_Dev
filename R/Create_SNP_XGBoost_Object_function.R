

#' For a single SNP, create an SNPFastImpute Object
#'
#' @param df the dataframe containing NAs, p columns of SNPs, n rows of samples.
#' @param a the column indicator of the SNP in the dataset.
#' @param size the windows size around the SNP to use as predictor variables. 
#' 
#' 
#' @details Basically the function do two different jobs. 
#' Using the known values for each SNP to predict the missing values for that 
#'    SNP. 
#'
#'
#' @return an object for the corresponding SNP.
#' 
#' If the SNP has missing value, then this object is a list with 6 elements:
#' 
#' 1. model_fit: indicator whether we need to fit a model for the SNP.
#' 2. SNP_position: position of the SNP.
#' 3. NA_positions: position of the missing values.
#' 4. train_data: samples that are not missing for this SNP. 
#' 5. train_lable: the values of the non-missing samples for this SNP
#' 6. pred_data: samples that are missing for this SNP.
#' 7. pred_label: an empty vector to store the future predicted labels. 
#' 
#' 
#' @export
#'
#' @examples
#' data("Test_df")
#' single_SNP_Obj <- Create_Single_SNP_Object(Test_df, 1, 10)
#' ## return a list of SNP related information for model building.
#' 
#' Create_Single_SNP_Object(Test_df, 50, 200)
#' ## Error message. 
#' ## Stop as the size of the windows are larger than the number of SNPs in the dataframe. 
#' 
#' Create_Single_SNP_Object(matrix(NA, 10, 5), 1, 1)
#' ## Should print a warning message that all samples for this SNP are NA's.
#' ## return a list containing a model_fit value equal to false. 
#' 
#' 
Create_Single_SNP_Object <- function(df, a, size) {
  ## Store the dimension of df
  n <- nrow(df)
  p <- ncol(df)
  ## check whether there are enough 
  if (size > p) stop("The size of the windows should be smaller than the number of SNPs.")
  
  ## Based on the size, decide the range of the predictor variables. 
  if (a == 1)  range = c((a+1): size)
  else if (a == p)  range = c((a-size): p)
  else if (a <= size/2) range = c(1:(a-1), (a+1): size)
  else if ((a + size/2) >= p)  range = c((a-(size - (p-a))):(a-1), (a+1): p)
  else range = c((a-size/2):(a-1), (a+1): (a + size/2))
  
  ## For this SNP, get which samples has missing value, which does not. 
  NA_sample <- which(is.na(df[, a]))
  NA_length <- length(NA_sample)

  if (NA_length == n) {
    warning(paste0("All samples for the No.", a, 
                 " SNP are NA's, we cannot build a model to predict the genotype for this SNP."))
    return(list(model_fit = F))
  }
  else {
    
    ## Dataset to train the model
    train_data <- df[-NA_sample, range, drop = F]
    ## Labels to train the model 
    train_label <- as.numeric(df[-NA_sample, a])
    
    ## Dataset to do prediction
    pred_data <- df[NA_sample, range, drop = F]
    pred_label <- as.numeric(df[NA_sample, a])
    
    ## return a list of information:
    return(list(model_fit = T, 
                SNP_position = a, 
                NA_positions = NA_sample,
                train_data = train_data, 
                train_label = train_label,
                pred_data = pred_data, 
                pred_label = pred_label))
  }
}




#' Input a dataframe, generate an object for the SNP missing value imputation with xgboost.
#'
#' @param df A dataframe of the SNPs which we need to impute missing values
#' @param size A windows size which we control for our imputation. Default value is 50.
#' 
#' @details 
#' This function takes in a matrix (p columns of SNPs, n rows of samples) with NAs and 
#' returns an object for our imputation function.  
#'
#' @return a list of list of 2 contents:  
#' 
#' 1. NA_cols: column index with missing values  
#' 
#' 2. Multiple_SNP_Object: corresponding SNP_Objects for these SNPs
#' 
#' @export
#'
#' @examples
#' df <- readRDS("Test_df.rds")
#' ## This is a short matrix containing 100 SNPs with 20 samples. 
#' 
#' xgboost_snp_obj <- Create_SNP_XGBoost_Object(df)
#' 
#' xgboost_snp_obj <- Create_SNP_XGBoost_Object(df, size = 200)
#' 
#' 
#' 
Create_SNP_XGBoost_Object <- function(df, size = 10) {
  ## To save time, we only impute SNP columns with NAs in the dataset.
  ## Get the dimension of df
  n <- nrow(df)
  p <- ncol(df)
  ## check whether there are enough 
  if (size > p) stop("The size of the windows should be smaller than the number of SNPs.")
  
  NA_cols <- which(apply(df, 2, function(x) sum(is.na(x))> 0))
  n_NA_cols <- length(NA_cols)
  ## Check whether there are missing values in the dataset.
  if (n_NA_cols == 0) {
    print("There is no missing value in the data set.")
  }
  else {
    ## Create an object to store the SNP objects for all SNPs (columns) with missing values.
    Multiple_SNP_Object <- list()
    for(i in 1:n_NA_cols) {
      Multiple_SNP_Object[[i]] <- Create_Single_SNP_Object(df = df, a = i, size = size)
    }
  }
  ## Return the column index with missing values, and the corresponding SNP_Objects for these SNPs.
  return(list(NA_cols = NA_cols, Multiple_SNP_Object = Multiple_SNP_Object))
}

