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
source('methods/extrapolation.R')
source('methods/create_weights.R')
source('methods/wt_bsplinesmoothing.R')

## Utilities -------------------------------------------------------------------
n <- dim(curves)[2]
q <- length(xlist)
reconst_fcts  <- find_obs_inc(Y = curves)

loc <- 2
t.points <- T.period

## Set values for smoothing
#breaks <- c(seq(0,0.5,0.05), seq(0.6,1,0.1), seq(2,10,1))
breaks <- c(seq(0,1,0.1), seq(2,10,0.5))

data <- list(dJB  = dJB,
             MAG  = MAG,
             SoF  = SoF,
             VS30 = VS30)

## Extrapolation
extrapolate   <- extrapolation(curves       = curves,
                               t.points     = T.period,
                               T_hp         = T_hp,
                               reconst_fcts = reconst_fcts)
curves.extrap <- extrapolate$curves.rec

## Construction of the weights
wgt         <- create_weights(curves.rec    = curves.extrap,
                              t.points      = t.points,
                              breaks        = seq(0,10,0.1),
                              loc.par       = loc,
                              reconst_fcts  = reconst_fcts,
                              Thp           = T_hp)

## Smoothing
smth        <- wt_bsplinesmoothing(curves   = curves.extrap,
                                   wgts.obs = wgt$wgts.obs,
                                   t.points = t.points,
                                   breaks   = breaks,
                                   lambda   = 1e-5,
                                   set.cb   = FALSE)
curves.extrap.fd <- smth$curves.fd
plot(curves.extrap.fd)

## Cross-validation ------------------------------------------------------------
source('methods/Regression/lambda_select.R')

B <- 5

MSE_cv <- list()
(Start.Time <- Sys.time())
for(b in 1:B)
{
  MSE_cv[[b]] <- lambda_select(b         = b,
                               B         = B,
                               curves    = curves.extrap,
                               curves.fd = curves.extrap.fd,
                               xlist     = xlist,
                               t.points  = t.points,
                               fPar      = smth$fPar,
                               events    = event.id,
                               wgts.fd   = wgt$wgts.fd)
  print(paste0('Iteration ', b,' of ', B))
}

save(MSE_cv, file='blist_options/MSE_cv_2-newestbreaks.RData')

##------------------------------------
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)
##------------------------------------

load('blist_options/MSE_cv_2-newestbreaks.RData')

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
    cum_sum <- cum_sum + MSE_cv[[i]][[reg]]
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
## to perform lambda selection are 1 (intercept) 2 (low magnitudes) and 6 (geometric attenuation)

## For all regressors except for the intercept take the min value and keep it fixed
## from now on


lambda.fixed <- numeric(9)
for(reg in 1:9)
{
  lambda.fixed[reg] <- lambda.vec[which(MSE_list[[reg]] == min(MSE_list[[reg]]))]
}

blist <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)

for(q in 1:9)
{
  blist[[q]]$lambda <- lambda.fixed[q]
}

save(blist, file='blist_options/blist-newestbreaks.RData')

