setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

library(fda)
library(devtools)
library(fastmatrix)
library(calculus)
library(ReconstPoFD)
library(snowfall)
library(progress)
library(ggplot2)
library(fields)
library(beepr)
library(tidyverse)
library(latex2exp)

rm(list=ls())
graphics.off()
cat("\014")

## Load Data -------------------------------------------------------------------
load('DATA/curves.RData')
load('DATA/t_period.RData')
load('DATA/obs.RData')
load('DATA/T_hp.RData')
load('DATA/xlist-logg.RData')
load('DATA/events.RData')
load('blist_options/log/blistt_latest.RData')

## Function for the evaluation of the MSE for each proposed method
source('methods/Reconstruction/methods_workflow_new.R')

# Logarithm of the period
t.points    <- log10(T.period)
t.points[1] <- -2.5
N <- length(t.points)

log.Thp     <- log10(T_hp)

## CV PER SCELTA DEI PESI -------------------------------------------------------
B     <- 10

method <- 'KLAl'

#vec.par <- c(2,5,10,15,20)
#vec.par <- "inf"
#vec.par <- c(50,100,1000)
#vec.par <- c(5,10,15,20,100)
vec.par <- c(100)
A <- length(vec.par)

MSE_cv <- array(data=0, dim=c(N, B, A))
MSE_glob.list <- list(list(),list(),list(),list(),list())
MSE_glob_bis.list <- list(list(),list(),list(),list(),list())
MSE_part.list <- list(list(),list(),list(),list(),list())
MSE_part_bis.list <- list(list(),list(),list(),list(),list())

beta1.list <- list()
beta2.list <- list()
beta3.list <- list()
beta6.list <- list()
beta7.list <- list()
beta8.list <- list()
beta9.list <- list()

(Start.Time <- Sys.time())
for(a in 1:A) # a=1
{
  fix.par <- vec.par[a]
  for(b in 1:B) # b=1
  {
    print(paste0("Par:", a, " b:",b))
    method_evaluation <- methods_workflow_new(b          = b,
                                              T.period   = T.period,
                                              t.points   = t.points,
                                              breaks     = t.points,
                                              T_hp       = T_hp,
                                              log.Thp    = log.Thp,
                                              curves     = curves,
                                              events     = event.id,
                                              B          = B,
                                              xlist      = xlist,
                                              blist      = blist,
                                              method     = method,
                                              fix.par    = fix.par,
                                              wgts.flag  = TRUE)
    
    MSE_cv[,b,a] <- method_evaluation$MSE_pw
    MSE_glob.list[[a]][[b]] <- method_evaluation$MSE_glob
    MSE_glob_bis.list[[a]][[b]] <- method_evaluation$MSE_glob_bis
    MSE_part.list[[a]][[b]] <- method_evaluation$MSE_part
    MSE_part_bis.list[[a]][[b]] <- method_evaluation$MSE_part_bis
    
    if(b==1)
    {
      beta1.est <- method_evaluation$beta_estimates[[1]]$fd
      beta2.est <- method_evaluation$beta_estimates[[2]]$fd
      beta3.est <- method_evaluation$beta_estimates[[3]]$fd
      beta6.est <- method_evaluation$beta_estimates[[6]]$fd
      beta7.est <- method_evaluation$beta_estimates[[7]]$fd
      beta8.est <- method_evaluation$beta_estimates[[8]]$fd
      beta9.est <- method_evaluation$beta_estimates[[9]]$fd
      
      beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
      beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
      beta3.est$coefs <- beta3.est$coefs %*% rep(1,B)
      beta6.est$coefs <- beta6.est$coefs %*% rep(1,B)
      beta7.est$coefs <- beta7.est$coefs %*% rep(1,B)
      beta8.est$coefs <- beta8.est$coefs %*% rep(1,B)
      beta9.est$coefs <- beta9.est$coefs %*% rep(1,B)
    } else {
      beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
      beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
      beta3.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
      beta6.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[6]]$fd$coefs, ncol=1)
      beta7.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[7]]$fd$coefs, ncol=1)
      beta8.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[8]]$fd$coefs, ncol=1)
      beta9.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[9]]$fd$coefs, ncol=1)
    }
    
  }
  beta1.list[[a]] <- beta1.est
  beta2.list[[a]] <- beta2.est
  beta3.list[[a]] <- beta3.est
  beta6.list[[a]] <- beta6.est
  beta7.list[[a]] <- beta7.est
  beta8.list[[a]] <- beta8.est
  beta9.list[[a]] <- beta9.est
}

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

#name.file <- paste0('Results/a-definition/MSE_20fold_nowgts.RData')
#name.file <- paste0('Results/a-definition/MSE_10fold_oldwgts.RData')
#name.file <- paste0('Results/a-definition/partial-MSE/MSE_10fold_oldwgts.RData')
#name.file <- paste0('Results/a-definition/partial-MSE/MSE_10fold_unwgt.RData')
name.file <- paste0('Results/a-definition/partial-MSE/MSE_10fold_0wgts.RData')

save(vec.par, MSE_cv, MSE_glob.list, MSE_glob_bis.list, MSE_part.list, MSE_part_bis.list, beta1.list, beta2.list, beta3.list,
     beta6.list, beta7.list, beta8.list, beta9.list, file=name.file)

## CV PER RECONSTRUCTION METHOD COMPARISON ---------------------------------------------
B       <- 10
fix.par <- 100

#method <- 'extrapolation'
method <- 'KLAl'

MSE_cv <- array(data=0, dim=c(N, B, 1))
MSE_glob.list <- list(list(),list(),list(),list(),list())
MSE_glob_bis.list <- list(list(),list(),list(),list(),list())
MSE_part.list <- list(list(),list(),list(),list(),list())
MSE_part_bis.list <- list(list(),list(),list(),list(),list())

beta1.list <- list()
beta2.list <- list()
beta3.list <- list()
beta6.list <- list()
beta7.list <- list()
beta8.list <- list()
beta9.list <- list()

(Start.Time <- Sys.time())
for(b in 1:B) # b=1
{
  print(paste0("b: ",b))
  method_evaluation <- methods_workflow_new(b          = b,
                                            T.period   = T.period,
                                            t.points   = t.points,
                                            breaks     = t.points,
                                            T_hp       = T_hp,
                                            log.Thp    = log.Thp,
                                            curves     = curves,
                                            events     = event.id,
                                            B          = B,
                                            xlist      = xlist,
                                            blist      = blist,
                                            method     = method,
                                            fix.par    = fix.par,
                                            wgts.flag  = TRUE)
  
  MSE_cv[,b,1] <- method_evaluation$MSE_pw
  MSE_glob.list[[1]][[b]] <- method_evaluation$MSE_glob
  MSE_glob_bis.list[[1]][[b]] <- method_evaluation$MSE_glob_bis
  MSE_part.list[[1]][[b]] <- method_evaluation$MSE_part
  MSE_part_bis.list[[1]][[b]] <- method_evaluation$MSE_part_bis
  
  if(b==1)
  {
    beta1.est <- method_evaluation$beta_estimates[[1]]$fd
    beta2.est <- method_evaluation$beta_estimates[[2]]$fd
    beta3.est <- method_evaluation$beta_estimates[[3]]$fd
    beta6.est <- method_evaluation$beta_estimates[[6]]$fd
    beta7.est <- method_evaluation$beta_estimates[[7]]$fd
    beta8.est <- method_evaluation$beta_estimates[[8]]$fd
    beta9.est <- method_evaluation$beta_estimates[[9]]$fd
    
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
    beta3.est$coefs <- beta3.est$coefs %*% rep(1,B)
    beta6.est$coefs <- beta6.est$coefs %*% rep(1,B)
    beta7.est$coefs <- beta7.est$coefs %*% rep(1,B)
    beta8.est$coefs <- beta8.est$coefs %*% rep(1,B)
    beta9.est$coefs <- beta9.est$coefs %*% rep(1,B)
  } else {
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta3.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
    beta6.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[6]]$fd$coefs, ncol=1)
    beta7.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[7]]$fd$coefs, ncol=1)
    beta8.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[8]]$fd$coefs, ncol=1)
    beta9.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[9]]$fd$coefs, ncol=1)
  }
  
}
beta1.list[[1]] <- beta1.est
beta2.list[[1]] <- beta2.est
beta3.list[[1]] <- beta3.est
beta6.list[[1]] <- beta6.est
beta7.list[[1]] <- beta7.est
beta8.list[[1]] <- beta8.est
beta9.list[[1]] <- beta9.est

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

name.file <- paste0('Results/rec-method-comparison/MSE_10fold_oldwgts_extrap.RData')

save(MSE_cv, MSE_glob.list, MSE_glob_bis.list, beta1.list, beta2.list, beta3.list,
     beta6.list, beta7.list, beta8.list, beta9.list, file=name.file)
