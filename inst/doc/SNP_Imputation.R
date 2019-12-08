## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(SNPFastImpute)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(nlme)

## ----vcf, echo=T, tidy=TRUE, tidy.opts=list(width.cutoff=40)------------------
data("full_vcf")
dim(full_vcf)
full_vcf[1:5, c(1:7, 9:14)]

## ----vcf2df, message=FALSE----------------------------------------------------
full_df <- vcf2df(full_vcf)

## ----showdf, echo=F-----------------------------------------------------------
full_df[1:4, 1:4]
sub <- 100

## ----subset-------------------------------------------------------------------
SNP_orig_sub <- full_df[, 1:sub]

## ----introduceNA--------------------------------------------------------------
## Create ratio sequence
ratios <- seq(0.05, 0.25, by = 0.05)
## Create list and vector to store the introduced NA data
SNP_NA_dfs <- list()
ratios_len <- length(ratios)
## Introduce NAs to reach the corresponding ratio
for (i in 1:ratios_len){
  SNP_NA_dfs[[i]] <- NA_Generator(SNP_orig_sub, ratios[i])
  print(SNP_NA_dfs[[i]]$NA_percent_generate)
}
names(SNP_NA_dfs) <- paste("missing", as.character(ratios))

## ----predict for different missing ratio, cache=F-----------------------------
## Create an empty list to store the filled data
df_fills <- list()
## Create an empty vectore to store the processing times.
proctimes <- rep(NA, ratios_len)
## Fill each data matrix using the parallel version of function. 
for(i in 1:ratios_len) {
  timeprocess <- system.time(
      df_fills[[i]] <- Impute_GenoType_XGBoost_para(SNP_NA_dfs[[i]]$SNP_NA_df)
  )
  proctimes[i] <- timeprocess[1]
}

## ----classification error calculation for different missing rate--------------
errors <- rep(NA, ratios_len)
for(i in 1:ratios_len){
  errors[i] <- classification_error(SNP_orig_sub, df_fills[[i]], SNP_NA_dfs[[i]]$NP_generate_positions)
}

## ----calculate missing ratio vs error rate and time, message=FALSE, echo = F----
errors_df <- data.frame(ratios, errors)
times_df <- data.frame(ratios, proctimes)

## ----plot missing ratio vs error rate and time, echo = F, eval = T, fig.width=7, message = F, tidy = T, tidy.opts=list(width.cutoff=35)----
p1 <- ggplot(errors_df, aes(x = ratios, y = errors)) + geom_point() + ylim(c(0, 0.2))
p2 <- ggplot(times_df, aes(x = ratios, y = proctimes)) + geom_point() + ylim(c(0, 0.5)) + ylab("process time")
ggarrange(p1, p2, nrow = 1)

## ----windows size list--------------------------------------------------------
sizes <- seq(10, 50, by = 10)
sizes_len <- length(sizes)
SNP_NA_df02 <- SNP_NA_dfs[[4]]$SNP_NA_df

## ----predict for different windows size, cache=F------------------------------
df_fills <- list()
proctimes <- rep(NA, ratios_len)
for(i in 1:sizes_len) {
  timeprocess <- system.time(
      df_fills[[i]] <- Impute_GenoType_XGBoost_para(SNP_NA_df02, size = sizes[i])
  )
  proctimes[i] <- timeprocess[1]
}

## ----classification error calculation for different windows size--------------
errors <- rep(NA, ratios_len)
for(i in 1:ratios_len){
  errors[i] <- classification_error(SNP_orig_sub, df_fills[[i]], SNP_NA_dfs[[i]]$NP_generate_positions)
}

## ----calculate windows size vs error rate and time, message=FALSE, echo = F----
errors_df <- data.frame(sizes, errors)
times_df <- data.frame(sizes, proctimes)

## ----plot windows size vs error rate and time, echo = F, eval = T, fig.width=7, message = F, tidy = T, tidy.opts=list(width.cutoff=35)----
p1 <- ggplot(errors_df, aes(x = sizes, y = errors)) + geom_point() + ylim(c(0, 0.2))
p2 <- ggplot(times_df, aes(x = sizes, y = proctimes)) + geom_point() + ylim(c(0, 0.5)) + ylab("process time") 
ggarrange(p1, p2, nrow = 1)

## ----predict for different missing ratio before and after considering correlation, cache=F, warning=F----
df_fills <- df_fills_corr <- list()
proctimes <- rep(NA, ratios_len)
for(i in 1:ratios_len) {
  timeprocess <- system.time(
      df_fills[[i]] <- Impute_GenoType_XGBoost_para(SNP_NA_dfs[[i]]$SNP_NA_df),
      df_fills_corr[[i]] <- Impute_GenoType_XGBoost_para(SNP_NA_dfs[[i]]$SNP_NA_df, select.corr = T)
  )
  proctimes[i] <- timeprocess[1]
}

## ----classification error calculation for different missing rate before and after considering correlation----
errors <- error_corr <- rep(NA, ratios_len)
for(i in 1:ratios_len){
  errors[i] <- classification_error(SNP_orig_sub, df_fills[[i]], SNP_NA_dfs[[i]]$NP_generate_positions)
  error_corr[i] <- classification_error(SNP_orig_sub, df_fills_corr[[i]], SNP_NA_dfs[[i]]$NP_generate_positions)
}

## ----calculate missing ratio vs error rate after considering correlation, message=FALSE, echo=F----
errors_df <- data.frame(ratios, errors)
error_corr_df <- data.frame(ratios, error_corr)

## ----plot missing ratio vs error rate after considering correlation, echo = F, eval = T, fig.width=7, message = F, tidy = T, tidy.opts=list(width.cutoff=35)----
p1 <- ggplot(errors_df, aes(x = ratios, y = errors)) + geom_point() + ylim(c(0, 1)) + ylab("error with windows not consider correlation")
p2 <- ggplot(error_corr_df, aes(x = ratios, y = error_corr)) + geom_point() + ylim(c(0, 1)) + ylab("error with windows consider correlation")
ggarrange(p1, p2, nrow = 1)

