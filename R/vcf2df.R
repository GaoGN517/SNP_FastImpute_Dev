#' Read a vcf file, output the corresponding dataframe. 
#'
#' @param filename The vcf file we want to read and generate the dataframe we need to 
#' generate the xgboost dataframe.
#'
#' @return A dataframe which we need to generate the xgboost data structure. 
#' @export
#'
#' @examples
#' filename <- "Test.vcf"
#' output_df <- vcf2df(filename)
#' 
#' 
vcf2df <- function(filename) {
  ## Read in the original vcf file as a dataframe.
  df <- read.table(filename)
  ## Greb information from the original dataframe.
  df_gt <- t(df[,10:ncol(df)]) 
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