

#' For a single SNP, create an SNPFastImpute Object
#'
#' @param df the dataframe containing NAs, p columns of SNPs, n rows of samples.
#' @param a the column indicator of the SNP in the dataset.
#' @param size the windows size around the SNP to use as predictor variables. 
#' @param cor.matrix A matrix storing the correlation of all SNPs in the dataframe. 
#' Defualt is NULL, which is to just use the SNPs around the target SNP to build model. 
#' When given this matrix, we pick the top n = size highest correlated columns to build model.
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
#' 8. windows_range: the range of surrounding SNPs used for model building.
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
#' corr <- cor(Test_df, method = "spearman", use = "pairwise.complete.obs")
#' Create_Single_SNP_Object(Test_df, 3, 20, cor.matrix = corr)
#' Create_Single_SNP_Object(Test_df, 3, 20)
#' 
Create_Single_SNP_Object <- function(df, a, size, cor.matrix = NULL) {
  ## Store the dimension of df
  n <- nrow(df)
  p <- ncol(df)
  ## check whether there are enough 
  if ((size + 1) > p) stop("The size of the windows should be smaller than the number of SNPs besides the one to predict.")
  
  if (!is.null(cor.matrix)) {
    corr.a <- cor.matrix[a, -1]
    corr.len <- sum(!is.na(corr.a))
    if(corr.len >= size) {
      range <- order(cor.matrix[a, ])[2:(size+1)]
      range.part2.size <- 0
    }
    else {
      range.part1 <- which(!is.na(corr.a))
      range.part2.size <- size - corr.len
    }
  }
  else {
    range.part1 <- vector()
    range.part2.size <- size
  }
  if (range.part2.size != 0) {
    position_vec <- seq(1, p, by = 1)
    if (length(range.part1)==0) rest.positions <- position_vec
    else rest.positions <- position_vec[-range.part1]
    a.position.in.rest <- rest.positions[rest.positions == a]
    rest.len <- length(rest.positions)
    ## Based on the size, decide the range of the predictor variables.
    if (a.position.in.rest == 1)  range.part2.pos <- c((a.position.in.rest+1): range.part2.size)
    else if (a.position.in.rest == rest.len)  range.part2.pos <- c((a.position.in.rest-range.part2.size): rest.len)
    else if (a.position.in.rest <= range.part2.size/2) range.part2.pos <- c(1:(a.position.in.rest-1), (a.position.in.rest+1): range.part2.size)
    else if ((a.position.in.rest + range.part2.size/2) >= rest.len)  range.part2.pos <- c((a.position.in.rest-(range.part2.size - (rest.len-a.position.in.rest))):(a.position.in.rest-1), (a.position.in.rest+1): rest.len)
    else range.part2.pos <- c((a.position.in.rest-range.part2.size/2):(a.position.in.rest-1), (a.position.in.rest+1): (a.position.in.rest + range.part2.size/2))
    range.part2 <- rest.positions[range.part2.pos]
    range <- unique(append(range.part1, range.part2))
  }
  
  
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
                pred_label = pred_label, 
                windows_range = range))
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
#' 
#'
#' @examples
#' 
#' ## This is a short matrix containing 100 SNPs with 20 samples. 
#' data("Test_df")
#' xgboost_snp_obj <- Create_SNP_XGBoost_Object(df)
#' 
#' 
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

