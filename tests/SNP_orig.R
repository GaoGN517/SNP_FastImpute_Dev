#' SNP subset data with low percent of missing values. 
#'
#' Data generated from the original larger dataset, with 
#'
#' @docType data
#'
#' @usage data(SNP_orig)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#' @keywords datasets
#'
#' @references 
#'
#' 
#'
#' @examples
#' data(SNP_orig)
#' dims(SNP_orig)
#' sum(is.na(SNP_orig)) / length(SNP_orig)
#' 
"SNP_orig"