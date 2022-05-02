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
load('DATA/xlist-logg.RData')
load('DATA/data.RData')
load('DATA/events.RData')

## Load functions --------------------------------------------------------------
source('methods/find_obs_inc.R')

## Utilities -------------------------------------------------------------------
n            <- dim(curves)[2]
q            <- length(xlist)
reconst_fcts <- find_obs_inc(Y = curves)
t.points     <- log10(T.period)
t.points[1]  <- -2.5
breaks       <- t.points

## Maintain only the fully observed curves
curves.full <- curves[,-reconst_fcts]
event.id    <- event.id[-reconst_fcts]
xlist.full  <- list()

for(i in 1:length(xlist))
{
  if (inherits(xlist[[i]], "fd"))
  {
    xlist.full[[i]]  <- xlist[[i]]
    xlist.full[[i]]$coefs  <- xlist[[i]]$coefs[,-reconst_fcts]
    
  } else if (inherits(xlist[[i]], "numeric")) {
    
    xlist.full[[i]]  <- xlist[[i]][-reconst_fcts]
    
  } else if (inherits(xlist[[i]], "matrix" )) {
    
    xlist.full[[i]]  <- xlist[[i]][-reconst_fcts,1]
    
  }
}

## Unweighted Smoothing
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-7,1, by=1)
lambda.vec <- sort(10^esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, curves.full, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
curves.fd <- smooth.basis(t.points, curves.full, fPar)$fd

## Cross-validation ------------------------------------------------------------
source('methods/Regression/lambda_select-unwgt.R')

B <- 10
fPar_reg <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-2)
exp.vec  <- seq(-6, 1, by=1)

MSE_cv <- list()
(Start.Time <- Sys.time())
for(b in 1:B) # b=1
{
  MSE_cv[[b]] <- lambda_select_unwgt(b         = b,
                                     B         = B,
                                     curves    = curves.full,
                                     curves.fd = curves.fd,
                                     xlist     = xlist.full,
                                     t.points  = t.points,
                                     fPar      = fPar_reg,
                                     exp.vec   = exp.vec,
                                     events    = event.id)
  print(paste0('Iteration ', b,' of ', B))
}

save(MSE_cv, file='blist_options/log/MSE_unwgt_10fold.RData')
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

#load('blist_options/log/MSE_unwgt.RData')

# per ogni regressore, per ogni lambda voglio valutare l'MSE
load('blist_options/log/MSE_unwgt_fPar2.RData')

## Per ogni regressore e per ogni lambda voglio tirare fuori il MSE e la sua sd
B <- 10
exp.vec    <- seq(-6, 1, by=1)
lambda.vec <- 10^exp.vec

MSE_list <- list()
MSEvec_list <- list()
MSE.mat <- matrix(data=0, nrow=length(lambda.vec), ncol=B)
for(reg in 1:9)
{
  for(b in 1:B)
  {
    MSE.mat[,b] <- unlist(MSE_cv[[b]][[reg]])
  }
  MSE_list[[reg]] <- MSE.mat
  MSEvec_list[[reg]] <- apply(MSE_list[[reg]], MARGIN=1, FUN=mean)
}

for(reg in 1:9)
{
  print(paste0('Regressor ', reg))
  mse <- apply(MSE_list[[reg]], MARGIN=1, FUN=mean)
  sigma <- apply(MSE_list[[reg]], MARGIN=1, FUN=sd)
  for(l in 1:length(lambda.vec))
  {
    print(paste0('Exp ', exp.vec[l], '   mse:', round(mse[l], digits=5), '  sd:', round(sigma[l], digits=5)))
  }
}


lambda.fixed <- numeric(9)
for(reg in 1:9)
{
  lambda.fixed[reg] <- lambda.vec[which(MSEvec_list[[reg]] == min(MSEvec_list[[reg]]))]
}
lambda.fixed

alex <- lambda.fixed

save(lambda.fixed, file='blist_options/log/lambda_fixedd.RData')

load('blist_options/log/lambda_fixedd.RData')
lambda.fixed
save(lambda.fixed, file='blist_options/log/lambda_fixed.RData')

# Se inizio con fPar con lambda 1e-2 invece che 1e-5 ottengo:
# lambda.fixed <- c(1e-06, 1e-01, 1e-03, 1e-02, 1e-02, 1e-02, 1e-06, 1e+01, 1e-02)

# Se inizio con fPar con lambda 1e-3 invece che 1e-5 ottengo:
# lambda.fixed <- c(1e-05, 1e-02, 1e-03, 1e-03, 1e-03, 1e-03, 1e-06, 1e+00, 1e-03)

blist <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)

for(q in 1:9)
{
  blist[[q]]$lambda <- lambda.fixed[q]
}

save(blist, file='blist_options/log/blistt_fPar3.RData')

## Rifaccio la CV per l'intercetta e c2 ----------------------------------------



source('methods/Regression/lambda_select_single.R')
B <- 5
exp.vec  <- seq(-6, 1, by=1)


## Intercept ------------------------------------------------------
MSE_cv <- list()
(Start.Time <- Sys.time())
for(b in 1:B) # b=1
{
  MSE_cv[[b]] <- lambda_select_single(b         = b,
                                      B         = B,
                                      curves    = curves.full,
                                      curves.fd = curves.fd,
                                      xlist     = xlist.full,
                                      t.points  = t.points,
                                      regs      = c(1),
                                      exp.vec   = exp.vec,
                                      events    = event.id)
  print(paste0('Iteration ', b,' of ', B))
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

View(MSE_cv)

lambda.vec <- 10^exp.vec

MSE_list <- list()
reg <- 1
cum_sum <- numeric(length(lambda.vec))
for(i in 1:B)
{
  cum_sum <- cum_sum + MSE_cv[[i]][[reg]]
}
MSE_list[[reg]] <- cum_sum/B

lambda.fixed <- lambda.vec[which(MSE_list[[reg]] == min(MSE_list[[reg]]))]
lambda.fixed

load('blist_options/log/blistt_latest.RData')
blist[[reg]]$lambda <- lambda.fixed
save(blist, file='blist_options/log/blistt_latest.RData')

## c2 ---------------------------------------------------------------
MSE_cv <- list()
(Start.Time <- Sys.time())
for(b in 1:B) # b=1
{
  MSE_cv[[b]] <- lambda_select_single_regs(b         = b,
                                           B         = B,
                                           curves    = curves.full,
                                           curves.fd = curves.fd,
                                           xlist     = xlist.full,
                                           t.points  = t.points,
                                           regs      = c(7),
                                           exp.vec   = exp.vec,
                                           events    = event.id)
  print(paste0('Iteration ', b,' of ', B))
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

lambda.vec <- 10^exp.vec

MSE_list <- list()
reg <- 7
cum_sum <- numeric(length(lambda.vec))
for(i in 1:B)
{
  cum_sum <- cum_sum + MSE_cv[[i]][[reg]]
}
MSE_list[[reg]] <- cum_sum/B

lambda.fixed <- lambda.vec[which(MSE_list[[reg]] == min(MSE_list[[reg]]))]
lambda.fixed

load('blist_options/log/blistt_latest.RData')
blist[[reg]]$lambda <- lambda.fixed
save(blist, file='blist_options/log/blistt_latest.RData')

