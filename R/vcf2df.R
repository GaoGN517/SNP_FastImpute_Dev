#' Read a vcf file, output the corresponding dataframe. 
#'
#' @param filename The vcf file we want to read and generate the dataframe we need to 
#' generate the xgboost dataframe.
#'
#' @return A dataframe which we need to generate the xgboost data structure. 
#' 
#' @details 
#' The input vcf file is directly from raw sequencing data.
#' It contains (n + 9) columns (Information for each SNP and corresponding SNP values for the samples)
#' and p rows (SNP positions). 
#' 
#' Starting from the 10th column are the information for the first sample. 
#' So we first remove the first 9 columns.
#' 
#' For our imputation, we need the p SNPs as features (columns), n samples as rows, so we need
#' to transpose the dataframe. 
#' 
#' In the input dataset, there are 2 values indicating the SNP types for each SNP position as there
#' are two alleles:
#' 0 (Wild type) and 1 (Mutate type). So the values can be "0/0", "0/1", "1/0", "1/1". Some of the
#' values might be missing. 
#' 
#' We sum up the two values at each position to one value to represent the corresponding SNP type. 
#' 
#' In the output data, each unit is the corresponding SNP type:
#' (1) 0: both alleles are mutations;
#' (2) 1: one of the alleles is a mutation, the other is wild type;
#' (3) 2: both alleles are wild type; 
#' (4) NA: at least one of the SNP type of the two alleles is missing. We need to predict the value
#'     for this position. 
#' 
#' 
#' @export
#'
#' @examples
#' filename <- "Test.vcf"
#' output_df <- vcf2df(filename)
#' ## This dataset has 112 samples and 338 SNP positions.
#' ## The original file has 121 columns and 338 rows.
#' 
#' ## Output should be a dataset with 112 rows and 338 columns. 
#' 

vcf2df <- function(filename) {
  ## Read in the original vcf file as a dataframe.
  df <- read.table(filename)
  
  ## Subset SNP type related information starting from the 10th column of the original dataframe.
  ## Transpose the dataframe so that SNPs are columns and samples are rows. 
  df_gt <- t(df[,10:ncol(df)]) 
  
  ## Greb the 0/1 values for SNP information
  df_gt <- gsub(":.*","",df_gt)
  df_gt <- gsub("\\|","/",df_gt)
  df_gt <- gsub("/"," ",df_gt)
  
  
  write.table(df_gt,"in.txt",sep=" ", quote=FALSE,row.names=FALSE,col.names=FALSE)
  df_gt=read.table("in.txt")
  df_gt=as.matrix(df_gt)
  
  class(df_gt) <- "numeric"
  
  table_df_gt <- apply(df_gt, 2, table)
  levels_table_df_gt <- unlist(lapply(table_df_gt, function(x) sum(as.numeric(names(x))>1)))
  check <- unique(ceiling(which(levels_table_df_gt>0)/2))
  
  n <- ncol(df_gt)
  df_out <- df_gt[,seq(1,n,by=2)] + df_gt[,seq(2,n,by=2)]
  colnames(df_out) <- df$V3
  df_out <- df_out[, -check]
  return(df_out)
}