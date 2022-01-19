setwd('/Users/teresabortolotti/Documents/R/WFDA')

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

rm(list=ls())
graphics.off()
cat("\014")

## Load functions -----------------------------------------
source('Regression/f_Regress.R')
source('Regression/my_predict_fRegress.R')
source('Regression/my_predict_fRegress_sing.R')
source('Regression/eval_MSE.R')

source('methods/find_obs_inc.R')
source('methods/interpolation.R')
source('methods/create_weights.R')
source('methods/wt_bsplinesmoothing.R')

source('methods/Regression/weighted_fRegress.R')
source('methods/Regression/my_predict_fRegress.R')
#source('Regression/my_predict_fRegress_sing.R')
source('methods/Regression/eval_MSE.R')

## Load data ----------------------------------------------
load('DATA/data_for_bootstrap.RData')
load('DATA/t_period.RData')
load('DATA/events.RData')


n <- dim(X)[1]
q <- dim(X)[2]
L <- curves.interp.fd$basis$nbasis
TT <- length(T.period)
xlist <- lapply(seq_len(ncol(X)), function(i) X[,i])

## Generate a bootstrap sample of functional coefficients --

# 1. Fit the regression
mod <- f_Regress(y            = curves.interp.fd,
                 xlist        = xlist,
                 betalist     = blist,
                 returnMatrix = FALSE,
                 method       = 'fRegress',
                 sep          = '.',
                 wgts         = wgt$wgts.fd)

# 2. evaluate the residuals
curves.interp.hat   <- my_predict_fRegress(mod          = mod,
                                           xlist        = xlist,
                                           t.points     = T.period,
                                           returnMatrix = FALSE)

res <- curves.interp.hat - curves.interp.fd
wgts.fd <- wgt$wgts.fd

save(mod, file='DATA/mod.RData')
# 3. Repeatedly sample from the empirical distribution of data
set.seed(14091996)
B <- 1000

obs <- seq(1,n)
B.list <- array(data=0, dim=c(TT,q,B))
                
B.mat <- matrix(nrow=TT, ncol=q)

pb <- progress_bar$new(total=B)
#pb$tick(0)
for(b in 1:B)
{
  obs.b <- sample(obs, replace=T)
  X.b   <- X[obs.b,]
  
  res.b <- res
  res.b$coefs <- res$coefs[,obs.b]
  curves.interp.fd.b <- curves.interp.hat + res.b
  wgts.fd.b <- wgts.fd
  wgts.fd.b$coefs <- wgts.fd$coefs[,obs.b]
  
  mod.b <- f_Regress(y            = curves.interp.fd.b,
                     xlist        = xlist,
                     betalist     = blist,
                     returnMatrix = FALSE,
                     method       = 'fRegress',
                     sep          = '.',
                     wgts         = wgt$wgts.fd)
  
  for(i in 1:q)
  {
    B.mat[,i] <- eval.fd(T.period, mod.b$betaestlist[[i]]$fd)
  }
  
  B.list[,,b] <- B.mat
  
  pb$tick()
}


#save(B.list, file='bootstrap/Blist_event.RData')
save(B.list, file='Blist.RData')

## Plot bootstrap ---------------------------------------------
library(wesanderson)
pal <- wes_palette('Cavalcanti1')


load('DATA/data_for_bootstrap.RData')
load('DATA/t_period.RData')
#load('Blist.RData')
#load('bootstrap/Blist_event.RData')
load('Blist.RData')

n <- dim(X)[1]
q <- dim(X)[2]
L <- curves.interp.fd$basis$nbasis
TT <- length(T.period)

coefs.names <- c('Intercept', 'B1', 'B2', 'F1', 'F2','C1','C2','C3','K')

for(j in 1:q)
{
  beta.j.mod <- mod$betaestlist[[j]]$fd
  beta.j.val <- eval.fd(T.period, beta.j.mod)
  
  png(file = paste0("Regressors/reg_",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
  matplot(T.period, B.list[,j,], type='l', col='grey80', main= coefs.names[j], font.main=1,       xlab='Period [s]', ylab='',
       cex.lab=2.3, cex.axis=2, cex.main=2.3)
  lines(T.period, beta.j.val, col=pal[2], lwd=3)
  abline(h=0, col=pal[5], lty=3, lwd=2)
  grid()
  dev.off()
}

## Functional boxplot
library(roahd)


for(j in 1:q)
{
  fdata <- fData(T.period, t(B.list[,j,]))
  png(file = paste0("Boxplots/",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
  fbplot(fdata, Depths='BD')
  dev.off()
}


## Regression comparison
# DATA
load("DATA/ITA18_regressors.RData")
a0 <- ITA18.regressors[,1]
b1 <- ITA18.regressors[,2]
b2 <- ITA18.regressors[,3]
f1 <- ITA18.regressors[,8]
f2 <- ITA18.regressors[,9]
c1 <- ITA18.regressors[,4]
c2 <- ITA18.regressors[,5]
c3 <- ITA18.regressors[,6]
k0 <- ITA18.regressors[,7]
X.ita18 <- cbind(a0,b1,b2,f1,f2,c1,c2,c3,k0)



for(j in 1:q)
{
  beta.j.mod <- mod$betaestlist[[j]]$fd
  beta.j.val <- eval.fd(T.period, beta.j.mod)
  
  png(file = paste0("Regressors_comparison/CB_",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
  matplot(T.period, B.list[,j,], type='l', col='grey80', main= coefs.names[j], font.main=1, xlab='Period [s]', ylab='',
          cex.lab=2.3, cex.axis=2, cex.main=2.3)
  lines(T.period, beta.j.val, col=pal[2], lwd=3)
  lines(T.period, X.ita18[,j], col=pal[1], lwd=3)
  abline(h=0, col=pal[5], lty=3, lwd=2)
  legend(0, -3.8, legend=c("Functional","Scalar"),
         col=c(pal[2],pal[1]), lty=c(1,1), lwd=c(2,2), cex=1.7)
  
  grid()
  dev.off()
}

t.plot <- log10(T.period)
t.plot[1] <- -3
xtick <- seq(-3,1,by=1)

for(j in 1:q)
{
  beta.j.mod <- mod$betaestlist[[j]]$fd
  beta.j.val <- eval.fd(T.period, beta.j.mod)
  
  png(file = paste0("Regressors_comparison/log/CB_",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
  matplot(t.plot, B.list[,j,], type='l', col='grey80', main= coefs.names[j], font.main=1, xlab='Period [s]', ylab='',
          cex.lab=2.3, cex.axis=2, cex.main=2.3, xaxt='n')
  lines(t.plot, beta.j.val, col=pal[2], lwd=4)
  lines(t.plot, X.ita18[,j], col=pal[1], lwd=4)
  abline(h=0, col=pal[5], lty=3, lwd=2)
  axis(side=1, at=xtick, labels = 10^xtick, cex.axis=2)
  grid()
  dev.off()
}