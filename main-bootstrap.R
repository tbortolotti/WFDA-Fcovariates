setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(fastmatrix)
library(R.matlab)
library(pbmcapply)
library(latex2exp)
library(calculus)
library(ReconstPoFD)
library(tidyverse)
library("xtable")
library(snowfall)
library(psych)
library(progress)
library(beepr)

rm(list=ls())
graphics.off()
cat("\014")

## Load functions -----------------------------------------
source('methods/Regression/weighted_fRegress.R')
source('methods/Regression/my_predict_fRegress.R')
source('methods/Regression/eval_MSE.R')

source('methods/find_obs_inc.R')
source('methods/extrapolation.R')
source('Simulation/methods/create_old_weights.R')
source('methods/wt_bsplinesmoothing.R')

## Load data -------------------------------------------------------------------
load('DATA/curves.RData')
load('DATA/t_period.RData')
load('DATA/obs.RData')
load('DATA/T_hp.RData')
load('DATA/xlist-logg.RData')
load('DATA/data.RData')
load('DATA/events.RData')
load('blist_options/log/blistt_latest.RData')
blist[[6]]$lambda <- 0.1

## Utilities -------------------------------------------------------------------
n <- dim(curves)[2]
q <- length(xlist)
reconst_fcts  <- find_obs_inc(Y = curves)
N <- length(T.period)

fix.par <- 100
t.points <- log10(T.period)
t.points[1] <- -2.5
breaks <- t.points

data <- list(dJB  = dJB,
             MAG  = MAG,
             SoF  = SoF,
             VS30 = VS30)

data.f <- data.frame(dJB = dJB, MAG = MAG, SoF = SoF, VS30 = VS30)

## Preprocessing ---------------------------------------------------------------
extrapolate   <- extrapolation(curves       = curves,
                               t.points     = T.period,
                               T_hp         = T_hp,
                               reconst_fcts = reconst_fcts)
curves.extrap <- extrapolate$curves.rec

## Construction of the weights
wgt       <- create_old_weights(curves.rec    = curves.extrap,
                                t.points      = t.points,
                                breaks        = breaks,
                                fix.par       = fix.par,
                                reconst_fcts  = reconst_fcts,
                                Thp           = log10(T_hp))

## Smoothing
smth             <- wt_bsplinesmoothing(curves   = curves.extrap,
                                        wgts.obs = wgt$wgts.obs,
                                        t.points = t.points,
                                        breaks   = breaks,
                                        lambda   = 1e-5,
                                        set.cb   = FALSE)
curves.extrap.fd <- smth$curves.fd
L <- curves.extrap.fd$basis$nbasis

## Generate a bootstrap sample of functional coefficients ----------------------

# 1. Fit the regression
mod <- weighted_fRegress(y            = curves.extrap.fd,
                         xfdlist      = xlist,
                         betalist     = blist,
                         wgts         = wgt$wgts.fd)

# 2. evaluate the residuals
curves.extrap.hat <- my_predict_fRegress(mod          = mod,
                                         xlist        = mod$xfdlist,
                                         t.points     = t.points)
res <- curves.extrap.hat - curves.extrap.fd
wgts.fd <- wgt$wgts.fd

# 3. Repeatedly sample from the empirical distribution of data
set.seed(140996)
B <- 1000   # Time consuming. One may decide to lower B and do some trials
#B <- 500

obs <- seq(1,n)
B.list <- array(data=0, dim=c(N,q,B))
                
B.mat <- matrix(nrow=N, ncol=q)

for(b in 1:B)
{
  print(paste0("Iteration ", b, " of ", B))
  obs.b <- sample(obs, replace=T)
  
  res.b <- res
  res.b$coefs <- res$coefs[,obs.b]
  curves.extrap.fd.b <- curves.extrap.hat + res.b
  
  mod.b <- weighted_fRegress(y            = curves.extrap.fd.b,
                             xfdlist      = xlist,
                             betalist     = blist,
                             wgts         = wgt$wgts.fd)
  
  for(i in 1:q)
  {
    B.mat[,i] <- eval.fd(t.points, mod.b$betaestlist[[i]]$fd)
  }
  
  B.list[,,b] <- B.mat
  
}

beep()

save(B.list, mod, n, q, L, N, T.period, file='Bootstrap/bootstrap_utils_latestt.RData')

