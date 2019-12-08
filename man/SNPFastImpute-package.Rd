\name{SNPFastImpute-package}
\alias{SNPFastImpute}
\alias{SNPFastImpute}
\docType{package}
\title{
\packageTitle{SNPFastImpute}
}
\description{
\packageDescription{SNPFastImpute}
}
\details{
Read in vcf file time or raw sequencing files, and impute missing single SNP values. This package currently containing a few
testing features that introduces additional missing positions to test different methods. 

Main functions: 
\itemize{
  \item Impute_GenoType_XGBoost
  \item Create_Single_SNP_Object
  \item classfication_error
  \item Generate_NA
  \item vcf2df
}
}
\author{
Guannan Gao

\strong{Maintainer}: \packageMaintainer{SNPFastImpute}
}
\references{
Kabisch, Maria, Ute Hamann, and Justo Lorenzo Bermejo. "Imputation of 
missing genotypes within LD-blocks relying on the basic coalescent and 
beyond: consideration of population growth and structure." BMC genomics 18.1 (2017): 798.

Tianqi Chen and Carlos Guestrin, "XGBoost: A Scalable Tree Boosting System",
22nd SIGKDD Conference on Knowledge Discovery and Data Mining, 2016, \url{https://arxiv.org/abs/1603.02754}
}


\keyword{ package }
