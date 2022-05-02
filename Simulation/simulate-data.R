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

rm(list=ls())
graphics.off()
cat("\014")

## Load Data -------------------------------------------------------------------
load('DATA/curves.RData')
load('DATA/t_period.RData')
load('DATA/obs.RData')
load('DATA/T_hp.RData')
load('DATA/data.RData')
load('DATA/ITA18_regressors.RData')

## Load Functions --------------------------------------------------------------
source('methods/find_obs_inc.R')
source('methods/extrapolation.R')
source('methods/Regression/my_predict_fRegress.R')

## Utilities -------------------------------------------------------------------
n <- dim(curves)[2]
t.points <- log10(T.period)
t.points[1] <- -2.5
t.points <- t.points + 2.5
N <- length(t.points)
reconst_fcts <- find_obs_inc(Y = curves)
breaks <- seq(range(t.points)[1], range(t.points)[2], 0.25)

## Create the covariates -------------------------------------------------------
# Fit the model against the geometric attenuation and the site term only

reg.D <- matrix(data=0, nrow=N, ncol=n)
for(t in 1:N)
{
  h <- ITA18.regressors$h.vec[t]
  R <- sqrt(dJB^2 + h^2)
  reg.D[t,] <- log10(R)
}
reg.S <- ifelse(VS30<=1500, log10(VS30/800), log10(1500/800))


basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-12,1, by=1)
lambda.vec <- sort(10^esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-3)
D.fd <- smooth.basis(t.points, reg.D, fPar)$fd

# d <- D.fd
# d$coefs <- as.matrix(D.fd$coefs[,1], ncol=1)
# plot(d)

onebasis <- create.constant.basis(range(t.points))
reg.S <- as.matrix(reg.S)
S.fd <- fd(matrix(reg.S,1,n), onebasis)

intercept <- as.matrix(rep(1,n))
intercept.fd <- fd(matrix(intercept,1,n), onebasis)

xlist <- list(interc.fd = intercept.fd,
              S.fd      = S.fd,
              D.fd      = D.fd)


## Evaluating the functional coefficients and the residuals of this simplified
## model -----------------------------------------------------------------------

## Extrapolation 
extrapolate   <- extrapolation(curves       = curves,
                               t.points     = T.period,
                               T_hp         = T_hp,
                               reconst_fcts = reconst_fcts)
curves.extrap <- extrapolate$curves.rec

## Smoothing
basis <- create.bspline.basis(rangeval=range(t.points), breaks=t.points, norder=4)
fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-3)
curves.extrap.fd <- smooth.basis(t.points, curves.extrap, fPar)$fd

# plot(curves.extrap.fd)

## B-list
fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-3)
blist <- list(fPar,fPar,fPar)

## Regression and beta estimation
mod <- fRegress(y            = curves.extrap.fd,
                xfdlist      = xlist,
                betalist     = blist,
                returnMatrix = FALSE,
                method       = 'fRegress',
                sep          = '.')

beta_estimates <- mod$betaestlist

## y_hat
curves.extrap.hat <- my_predict_fRegress(mod          = mod,
                                         xlist        = xlist,
                                         t.points     = t.points)

res <- curves.extrap.hat - curves.extrap.fd
res$coefs <- res$coefs/12

reg.info <- list(t.points       = t.points,
                 xlist          = xlist,
                 beta_estimates = beta_estimates,
                 res            = res)
save(reg.info, file='Simulation/DATA/reg_info.RData')

## DATA SIMULATION -------------------------------------------------------------
## Utilities
set.seed(140996)
n.sim <- 100

## Scalar covariate: x1
K.fd <- xlist[[2]]
K.val <- eval.fd(1, K.fd)
mu1 <- mean(K.val)
sd1 <- sd(K.val)
x1 <- rnorm(n.sim, 5*mu1, sd1)

## Functional covariate: x2
x2.fd <- xlist[[3]]
x2.fd$coefs <- xlist[[3]]$coefs[,1:n.sim]



pca <- pca.fd(x2.fd, nharm=6, centerfns=TRUE)
harm1 <- pca$harmonics[1]

x2.tilde <- x2.fd

for(i in 1:n.sim) #i=1
{
  x2i <- x2.fd
  x2i$coefs <- as.matrix(x2.fd$coefs[,i], ncol=1)
  
  score.i <- pca$scores[i,1]
  
  util.fd <- x2i - score.i*harm1
  
  x2.tilde$coefs[,i] <- util.fd$coefs
}

## espandi solo alcuni dei coefficienti e vediamo que pasa
fact.to.exp <- matrix(data=runif(dim(xlist[[3]]$coefs)[1]*n.sim, min=1, max=1.2),
                      ncol=n.sim)

terence <- xlist[[3]]
terence$coefs <- xlist[[3]]$coefs[,1:n.sim]
terence$coefs <- terence$coefs * fact.to.exp


# A <- terence$coefs
# vec <- fact.to.exp
# 
# aiuto <- A * vec
# 
# terence <- xlist[[3]]
# terence$coefs <- sweep(xlist[[3]]$coefs[,1:n.sim], MARGIN=2, fact.to.exp, '*')

# 
# x2.1 <- xlist[[3]]
# x2.1$coefs <- as.matrix(xlist[[3]]$coefs[,1], ncol=1)
# 
# terence.1 <- terence
# terence.1$coefs <- as.matrix(terence$coefs[,1], ncol=1)
# 
# par(mfrow=c(1,2))
# plot(x2.1)
# plot(terence.1)
# 
# par(mfrow=c(1,2))
# plot(xlist[[3]])
# plot(terence)


expanded.coefs <- terence$coefs
expanded.coefs <- x2.tilde$coefs
mu.vec <- as.numeric(apply(X=expanded.coefs, MARGIN=1, FUN=mean))
cov.mat <- var(t(expanded.coefs))
x2.fd$coefs <- t(mvrnorm(n.sim, mu.vec, cov.mat))

dev.off()
plot(x2.fd)
# 
# alex <- x2.fd
# alex$coefs <- as.matrix(x2.fd$coefs[,1], ncol=1)
# plot(alex)

## FPCA of the residuals
pca <- pca.fd(res, nharm=6, centerfns=TRUE)

# # scree plot
# layout(rbind(1,2))
# plot(pca$values, xlab='', ylab='Eigenvalues', pch=19)
# plot(cumsum(pca$values)/sum(pca$values), xlab='', ylab='CPV', ylim=c(0.8,1), pch=19)
# abline(h=0.99, col='red')
# 
# sum(pca$values[1:3])/sum(pca$values)
# 
# # First three FPCs
# layout(cbind(1,2,3))
# plot(pca$harmonics[1], t.points, ylab='FPC1')
# plot(pca$harmonics[2], t.points, ylab='FPC2')
# plot(pca$harmonics[3], t.points, ylab='FPC3')
#
# par(mfrow=c(1,3))
# plot.pca.fd(pca, nx=100, pointplot=TRUE, harm=c(1,2,3), expand=0, cycle=FALSE)

## Costruisco i residui-tilde come somma delle armoniche per gli scores corrispondenti 
## opportunamente ri-estratti

harm1 <- pca$harmonics[1]
harm1$coefs <- harm1$coefs %*% rep(1,n.sim)

harm2 <- pca$harmonics[2]
harm2$coefs <- harm2$coefs %*% rep(1,n.sim)

pca.scores <- pca$scores[,1:2]
scores.mu.vec <- as.numeric(apply(X=pca.scores, MARGIN=2, FUN=mean))
scores.cov.mat <- var(pca.scores)
new.scores <- mvrnorm(n.sim, scores.mu.vec, scores.cov.mat)

# res_tilde <- new.scores[,1]*harm1 + new.scores[,2]*harm2 + new.scores[,3]*harm3 + 
#   new.scores[,4]*harm4 + new.scores[,5]*harm5 + new.scores[,6]*harm6
res_tilde <- new.scores[,1]*harm1 + new.scores[,2]*harm2


## Construction of y -----------------------------------------------------------
beta0 <- beta_estimates[[1]]$fd
#beta0$coefs <- beta_estimates[[1]]$fd$coefs/3
beta1 <- beta_estimates[[2]]$fd
y1 <- beta0 + x1[1]*beta1

par(mfrow=c(1,3))
plot(beta0)
title(main='Intercept')
plot(beta1)
title(main='beta 1')
plot(y1, ylab='values')
title(main='y')

x21 <- x2.fd
x21$coefs <- as.matrix(x2.fd$coefs[,1], ncol=1)
beta2 <- beta_estimates[[3]]$fd
y1.bis <- y1 + x21*beta2

par(mfrow=c(2,2))
plot(x21)
title(main='x2')
plot(beta2)
title(main='beta 2')
plot(y1, ylab='values')
title(main='y')
plot(y1.bis, ylab='values')
title(main='y bis')


res1 <- res_tilde
res1$coefs <- as.matrix(res_tilde$coefs[,1], ncol=1)
y1.final <- y1.bis + res1

par(mfrow=c(2,2))
plot(y1)
title(main='y')
plot(y1.bis)
title(main='y bis')
plot(res1)
title(main='residual')
plot(y1.final)
title(main='y final')

par(mfrow=c(2,3))
plot(beta1, ylab='values')
title(main='beta 1')
plot(beta2, ylab='values')
title(main='beta 2')
plot(x21, ylab='values')
title(main='x2')
plot(y1.bis, ylab='values')
title(main='y')
plot(res1, ylab='values', ylim=c(-1.5,2))
title(main='scaled regression residual')
plot(y1.final, ylab='values')
title(main='y + scaled regression residual')

#y
#residual
#y noisy

beta0 <- beta_estimates[[1]]$fd
#beta0$coefs <- beta0$coefs/3
beta1 <- beta_estimates[[2]]$fd
beta2 <- beta_estimates[[3]]$fd

for(i in 1:n.sim) #i=1
{
  x2i <- x2.fd
  x2i$coefs <- as.matrix(x2.fd$coefs[,i], ncol=1)
  
  resi <- res_tilde
  resi$coefs <- as.matrix(res_tilde$coefs[,i], ncol=1)
  
  yi <- beta0 + x1[i]*beta1 + x2i*beta2 + resi
  yi.nores <- beta0 + x1[i]*beta1 + x2i*beta2
  if(i==1)
  {
    y.fd <- yi
    y.fd$coefs <- yi$coefs %*% rep(1,n.sim)
    
    y.fd.nores <- yi.nores
    y.fd.nores$coefs <- yi.nores$coefs %*% rep(1,n.sim)
  }
  
  y.fd$coefs[,i] <- yi$coefs
  y.fd.nores$coefs[,i] <- yi.nores$coefs
}

dev.off()
plot(y.fd, ylab='values')
title(main='Curves with regression errors')
plot(y.fd.nores, ylab='values')
title(main='Curves')

## Add the smoothing errors ----------------------------------------------------
## Add an iid error to the curves
t.grid <- seq(range(t.points)[1], range(t.points)[2], 0.25)

## CASO 1
## Curva definita parzialmente
left.bound <- 1.5
right.bound <- 3.5
Ui <- runif(n.sim, min=left.bound,max=right.bound)
Pi <- rbinom(n.sim, size=1, prob=0.4)
Pi.logic <- ifelse(Pi==1, TRUE, FALSE)
T_hp <- rep(right.bound, n.sim)
T_hp[Pi.logic] <- Ui[Pi.logic]

y.noisy <- eval.fd(t.grid, y.fd)
for(i in 1:n.sim)
{
  y.noisy[,i] <- y.noisy[,i] + rnorm(length(t.grid), 0, 0.1)
}

y.no.noise <- eval.fd(t.grid, y.fd)
par(mfrow=c(1,2))
matplot(t.grid, y.no.noise, type='l', main='Sampled curves - without noise', ylab='values', xlab='t')
matplot(t.grid, y.noisy, type='l', main='Sampled curves - iid noise', ylab='values', xlab='t')
y.po <- y.noisy
for(i in 1:n.sim)
{
  t.idx <- (t.grid>T_hp[i])
  y.po[t.idx,i] <- NA
}

obs.inc <- (T_hp<right.bound)
matplot(t.grid, y.po[,obs.inc], type='l', main='Partially observed functional data', ylab='values', xlab='t')

perc <- numeric(length(t.grid))
for(j in 1:length(t.grid))
{
  perc[j] <- sum(T_hp>=t.grid[j])/n.sim
}
plot(t.grid, perc, type='l', col='darkblue', lwd=2, ylim=c(0,1),
     main='Fraction of observed curves per instant', xlab='t',ylab='')
abline(h=0.6, lty=2, lwd=2, col='red')
points(t.grid, perc, pch=19)

curves <- y.po
curves.true <- y.noisy
curves.true.fd <- y.fd.nores
t.points <- t.grid
save(t.points, T_hp, curves, curves.true, curves.true.fd, regression.list, file='Simulation/DATA/CASE1.RData')


## CASO 2
## Curva con rumore che aumenta da un certo istante di tempo in poi
## In questo caso aumentiamo la porzione di curve che a un certo punto 
## iniziano a presentare un rumore alto

## Come prima prova, imponiamo che l'aumento di rumorosità si verifichi in un
## istante che è uguale per tutte le curve che soffrono di questo problema

noise.bound <- 1.5
right.bound <- 3.5
Pi <- rbinom(n.sim, size=1, prob=0.5)
Pi.logic <- ifelse(Pi==1, TRUE, FALSE)
T_hp <- rep(right.bound, n.sim)
T_hp[Pi.logic] <- noise.bound
obs.noisy <- (T_hp<right.bound)

y.vn <- eval.fd(t.grid, y.fd)
for(j in 1:length(t.grid))
{
  for(i in 1:n.sim)
  {
    if(T_hp[i]>=t.grid[j])
    {
      y.vn[j,i] <- y.vn[j,i] + rnorm(1,0,0.1)
    } else {
      y.vn[j,i] <- y.vn[j,i] + rnorm(1,0,0.5)
    }
  }
}

par(mfrow=c(1,2))
matplot(t.grid, y.vn[,!obs.noisy], type='l', main='Data with "natural" noise level',
        ylim=c(-4,4))
matplot(t.grid, y.vn, type='l', main='Data with variable noise level - all data',
        ylim=c(-4,4))

curves <- y.vn
t.points <- t.grid
save(t.points, T_hp, curves, regression.list, file='Simulation/DATA/CASE2.RData')





