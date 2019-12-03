
#' Function to introduce additional NAs for prediction error calculation. 
#'
#' @param SNP_df Original data of SNP. It may contain some SNP
#' @param percent The final percentage of total NA's we want in the dataset. 
#'
#' @return a list of variables:  
#' 1. SNP_NA_df: the dataset with additional SNPs.   
#' 2. NA_percent_orig: original missing percentable in the dataset.   
#' 3. NA_percent_generate: introduced NA percentage.   
#' 4. NP_generate_positions: the positions of the introduced NA's.   
#' @export
#'
#' @examples
#' SNP_df <- readRDS("data/SNP_orig_sub.rds")
#' ## original NA ratio is 0.018
#' SNP_NA_df <- NA_Generator(SNP_df, 0.2)
#' ## Introduced another 18% of NAs. 
#' # saveRDS(SNP_NA_df, file = "data/SNP_NA_df_list.rds")
#' SNP_NA_df <- NA_Generator(SNP_df, 0.01)
#' ## Should report an error: alread has 1.8% NA's, higher than the 
#' ## target percentage. 
#' 
NA_Generator <- function(SNP_df, percent) {
  m <- nrow(SNP_df)
  n <- ncol(SNP_df)
  ## Get the NA positions
  NA_positions <- which(is.na(SNP_df))
  
  ## Get the not NA positions and counts.
  notNA_positions <- which(!is.na(SNP_df))
  notNA_len <- length(notNA_positions)
  
  ## Calculate the original percentage of NA's
  NA_percent_orig <- length(NA_positions) / (m * n)
  
  ## Calculate the percentage and counts of NA we need to generate. 
  NA_percent_generate <- percent - NA_percent_orig
  size = notNA_len * NA_percent_generate
  
  ## If the size of the number of NA to be generated is smaller than 1, then
  ## we cannot introduce any more NA's under current required total percentage.
  if(size < 1) stop(paste("There are already", NA_percent_orig * 100,
                  "% of NAs in the dataset. Cannot generate more."))
  
  ## Introduce random NAs to the original dataset.
  NP_generate_positions <- sample(notNA_positions, size = size)
  SNP_NA_df <- SNP_df
  SNP_NA_df[NP_generate_positions] <- NA
  
  ## Return the list of related information. 
  return(list(SNP_NA_df = SNP_NA_df, NA_percent_orig = NA_percent_orig, 
              NA_percent_generate = NA_percent_generate, NP_generate_positions = NP_generate_positions))
}
