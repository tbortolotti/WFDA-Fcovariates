setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(fastmatrix)
library(latex2exp)
library(calculus)
library(ReconstPoFD)
library(tidyverse)
library("xtable")
library(snowfall)
library(psych)
library(progress)

rm(list=ls())
graphics.off()
cat("\014")

## Load Data -------------------------------------------------------------------
load('DATA/curves.RData')
load('DATA/t_period.RData')
load('DATA/T_hp.RData')
load('DATA/xlist.RData')
load('DATA/data.RData')
load('DATA/events.RData')

## Load functions --------------------------------------------------------------
source('methods/find_obs_inc.R')
source('methods/interpolation.R')
source('methods/create_weights.R')
source('methods/wt_bsplinesmoothing.R')

## Utilities -------------------------------------------------------------------
n <- dim(curves)[2]
q <- length(xlist)
reconst_fcts  <- find_obs_inc(Y = curves)

## Interpolation
interpolate   <- interpolation(curves       = curves,
                               t.points     = T.period,
                               T_hp         = T_hp,
                               reconst_fcts = reconst_fcts)
curves.interp <- interpolate$curves.rec

## Construction of the weights
loc         <- 0.3
t.points    <- log10(T.period)
t.points[1] <- -3
wgt         <- create_weights(curves        = curves.interp,
                              t.points      = t.points,
                              loc.par       = loc,
                              reconst_fcts  = reconst_fcts,
                              T_hp          = log10(T_hp))

## Smoothing
smth        <- wt_bsplinesmoothing(curves   = curves.interp,
                                   wgts.obs = wgt$wgts.obs,
                                   t.points = t.points,
                                   lambda   = 1e-5,
                                   set.cb   = FALSE)
curves.interp.fd <- smth$curves.fd

## Parallel cross-validation ---------------------------------------------------
source('methods/Regression/lambda_select.R')

B <- 5

MSE_cv_np <- list()
(Start.Time <- Sys.time())
for(b in 3:B)
{
  MSE_cv_np[[b]] <- lambda_select(b         = b,
                                  B         = B,
                                  curves    = curves.interp,
                                  curves.fd = curves.interp.fd,
                                  xlist     = xlist,
                                  t.points  = t.points,
                                  fPar      = smth$fPar,
                                  events    = event.id,
                                  wgts.fd   = wgt$wgts.fd)
  print(paste0('Iteration ', b,' of ', B))
}
sfStop()

save(MSE_cv_np, file='blist_options/MSE_cv_03_redo.RData')

##------------------------------------
End.Time <- Sys.time()
## Run-time:
round(End.Time - Start.Time, 2)
##------------------------------------

load('blist_options/MSE_cv_03_redo.RData')

# per ogni regressore, per ogni lambda voglio valutare l'MSE
fPar <- smth$fPar
exp.vec    <- seq(log10(fPar$lambda)-1, max(log10(fPar$lambda)+6,2), by=1)
lambda.vec <- 10^exp.vec

MSE_list <- list()
for(reg in 1:9)
{
  cum_sum <- numeric(length(lambda.vec))
  for(i in 1:B)
  {
    cum_sum <- cum_sum + MSE_cv_np[[i]][[reg]]
  }
  MSE_list[[reg]] <- cum_sum/B
}

sd.MSE <- numeric(length(MSE_list))
for(reg in 1:9)
{
  sd.MSE[reg] <- sd(MSE_list[[reg]])
}
sd.MSE

for(reg in 1:9)
{
  print(paste0('reg ',reg))
  print(MSE_list[[reg]])
}

MSE_rel <- list()
for(reg in 1:9)
{
  MSE_rel[[reg]] <- MSE_list[[reg]]/min(MSE_list[[reg]])
}

for(reg in 1:9)
{
  print(paste0('reg ',reg))
  print(MSE_rel[[reg]])
}

## From this analysis I see that the most important regressors on which
## to perform lambda selection are 1 (intercept) 2 (low magnitudes) and 3 (high magnitudes)

## For all regressors except for the intercept take the min value and keep it fixed
## from now on


lambda.fixed <- numeric(9)
for(reg in 1:9)
{
  lambda.fixed[reg] <- lambda.vec[which(MSE_list[[reg]] == min(MSE_list[[reg]]))]
}

# I force some coefficients to be a little more rough
lambda.fixed[3] <- lambda.vec[6]
lambda.fixed[5] <- lambda.vec[7]
lambda.fixed[8] <- lambda.vec[4]

# and force others to be smoother
lambda.fixed[6] <- lambda.vec[6]
lambda.fixed[7] <- lambda.vec[6]

save(lambda.fixed, file='blist_options/lambda_fixed_redo.RData')


blist <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)

for(q in 1:9)
{
  blist[[q]]$lambda <- lambda.fixed[q]
}

save(blist, file='blist_options/blist_prova.RData')

