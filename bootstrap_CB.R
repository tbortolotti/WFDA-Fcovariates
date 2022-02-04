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
library(boot)

rm(list=ls())
graphics.off()
cat("\014")

load('Bootstrap/bootstrap_utils.RData')

overall.conf <- list()
for(j in 1:q)
{
  overall.conf[[j]] <- envelope(mat = t(B.list[,j,]), level=0.95)$overall
}

# overall.conf.50 <- list()
# for(j in 1:q)
# {
#   overall.conf.50[[j]] <- envelope(mat = t(B.list[,j,]), level=0.50)$overall
# }


## Colors ----------------------------------------------------------------------
library(wesanderson)
pal <- wes_palette('Cavalcanti1')

pal2.fade <- adjustcolor(pal[2], alpha = 0.3)
pal3.fade <- adjustcolor(pal[3], alpha = 0.5)
pal4.fade <- adjustcolor(pal[4], alpha = 0.5)

## Regression comparison -------------------------------------------------------
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

coefs.names <- c('Intercept', 'B1', 'B2', 'F1', 'F2','C1','C2','C3','K')
names.coefs <- c('a','b1','b2','f1','f2','c1','c2','c3','k')

t.points <- log10(T.period)
t.points[1] <- -2.4
xtick <- seq(-2,1,by=1)

for(j in 1:q)
{
  beta.j.mod <- mod$betaestlist[[j]]$fd
  beta.j.val <- eval.fd(t.points, beta.j.mod)
  
  y.min <- min(c(overall.conf[[j]][2,], X.ita18[,j]))
  y.max <- max(c(overall.conf[[j]][1,], X.ita18[,j]))
  ylim <- c(y.min, y.max)
  png(file = paste0("Bootstrap/novel/Regressors-comparison/CB_",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
  plot(t.points, beta.j.val, type='l', col=pal[2], lwd=4, main= names.coefs[j], font.main=1, xlab='Period [s]', ylab='',
        cex.lab=2.8, cex.axis=2.5, cex.main=2.8, xaxt='n', ylim=ylim)
  points(t.points, beta.j.val, pch=19, col=pal[2])
  lines(t.points, X.ita18[,j], col=pal[1], lwd=4)
  points(t.points, X.ita18[,j], pch=19, col=pal[1])
  for(i in 1:37)
  {
    t <- t.points[i]
    c.sup <- overall.conf[[j]][1,i]
    c.inf <- overall.conf[[j]][2,i]
    lines(c(t, t), c(c.inf, c.sup), lwd=6, col=pal4.fade)
  }
  abline(h=0, col=pal[5], lty=5, lwd=2)
  axis(side=1, at=c(-2.4,xtick), labels = c(0, 10^xtick), cex.axis=2.5)
  grid()
  dev.off()
}


## Regression coefficients boxplots --------------------------------------------
library(roahd)
library(wesanderson)
pal <- wes_palette('Cavalcanti1')

pal.gb <- wes_palette('GrandBudapest2')
pal.gb.fade <- adjustcolor(pal.gb[2], alpha = 0.6)

t.points <- log10(T.period)
t.points[1] <- -2.4
xtick <- seq(-2,1,by=1)

names.coefs <- c('a','b1','b2','f1','f2','c1','c2','c3','k')
coefs.names <- c('Intercept', 'B1', 'B2', 'F1', 'F2','C1','C2','C3','K')

for(j in 1:q)
{
  beta.j.mod <- mod$betaestlist[[j]]$fd
  beta.j.val <- eval.fd(t.points, beta.j.mod)
  
  y.min <- min(overall.conf[[j]][2,])
  y.max <- max(overall.conf[[j]][1,])
  ylim <- c(y.min, y.max)
  
  png(file = paste0("Bootstrap/novel/Boxplots/reg_",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
  par(mar=c(4.5, 4, 2.5, 1)+.1)
  fda::fbplot(fit=B.list[,j,], x=t.points, method="MBD", color=pal.gb.fade, xlab='Period [s]',
              xlim=range(t.points), ylim=ylim, cex.axis=2.8, cex.lab=2.8, ylab='', xaxt='n')
  title(main=names.coefs[j], font.main=1, cex.main=2.8)
  lines(t.points, beta.j.val, lwd=4, col=pal[5], xaxt='n', yaxt='n')
  abline(h=0, col=pal[2], lty=5, lwd=3)
  axis(side=1, at=c(-2.4,xtick), labels = c(0,10^xtick), cex.axis=2.5)
  #axis(side=2, cex.axis=2.5)
  grid()
  dev.off()
}

## Regression coefficients boxplots and ITA18 comparison -----------------------
library(roahd)
library(wesanderson)
pal <- wes_palette('Cavalcanti1')

pal.gb <- wes_palette('GrandBudapest2')
pal.gb.fade <- adjustcolor(pal.gb[2], alpha = 0.6)


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

t.points <- log10(T.period)
t.points[1] <- -2.4
xtick <- seq(-2,1,by=1)

names.coefs <- c('a','b1','b2','f1','f2','Geometric attenuation','Geometric attenuation','Anelastic attenuation','k')
coefs.names <- c('Intercept', 'B1', 'B2', 'F1', 'F2','C1','C2','C3','K')

for(j in 1:q)
{
  beta.j.mod <- mod$betaestlist[[j]]$fd
  beta.j.val <- eval.fd(t.points, beta.j.mod)
  
  y.min <- min(overall.conf[[j]][2,])
  y.max <- max(overall.conf[[j]][1,])
  ylim <- c(y.min, y.max)
  
  png(file = paste0("Bootstrap/novel/Boxplots-and-comparison/reg_",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
  par(mar=c(4.5, 4, 2.5, 1)+.1)
  fda::fbplot(fit=B.list[,j,], x=t.points, method="MBD", color=pal.gb.fade, xlab='Period [s]',
              xlim=range(t.points), ylim=ylim, cex.axis=1.8, cex.lab=1.8, ylab='', xaxt='n')
  title(main=names.coefs[j], font.main=1, cex.main=1.8)
  lines(t.points, beta.j.val, lwd=4, col=pal[5], xaxt='n', yaxt='n')
  lines(t.points, X.ita18[,j], col=pal[1], lwd=4)
  points(t.points, X.ita18[,j], pch=19, col=pal[1])
  abline(h=0, col=pal[2], lty=5, lwd=3)
  axis(side=1, at=c(-2.4,xtick), labels = c(0,10^xtick), cex.axis=2.5)
  #axis(side=2, cex.axis=2.5)
  #legend(0, 0.6, legend=c("Functional", "Scalar", TeX("$M_w= 4.8$"),TeX("$M_w= 5.8$"), TeX("$M_w= 7.0$")),
  #      col=c('black','black',pal[3], pal[1],pal[5]), lwd=c(3,3,0,0,0), lty=c(1,5,0,0,0), pch=c(NA,NA,16,16,16), cex=2.2)
  
  grid()
  dev.off()
}

## Plots per la discussione ----------------------------------------------------

## C2
j <- 7

beta.j.mod <- mod$betaestlist[[j]]$fd
beta.j.val <- eval.fd(t.points, beta.j.mod)

y.min <- min(overall.conf[[j]][2,])
y.max <- max(overall.conf[[j]][1,])
ylim <- c(y.min, y.max)

png(file = paste0("Bootstrap/novel/Boxplots_and_comparison/reg_",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
par(mar=c(4.5, 4, 2.5, 1)+.1)
fda::fbplot(fit=B.list[,j,], x=t.points, method="MBD", color=pal.gb.fade, xlab='Period [s]',
            xlim=range(t.points), ylim=ylim, xaxt='n',cex.axis=1.8, cex.lab=1.8, ylab='')
title(main=names.coefs[j], font.main=1, cex.main=1.8)
lines(t.points, beta.j.val, lwd=4, col=pal[5], xaxt='n', yaxt='n')
lines(t.points, X.ita18[,j], col=pal[1], lwd=4)
points(t.points, X.ita18[,j], pch=19, col=pal[1])
abline(h=0, col=pal[2], lty=5, lwd=3)
axis(side=1, at=c(-2.4,xtick), labels = c(0,10^xtick), cex.axis=1.8)
#axis(side=2, cex.axis=2.5)
legend(-0.2, -1.62, legend=c("Functional", "Scalar"),
       col=c(pal[5],pal[1]), lwd=c(3,3), lty=c(1,1), cex=1.8)

grid()
dev.off()

## C3
j <- 8

beta.j.mod <- mod$betaestlist[[j]]$fd
beta.j.val <- eval.fd(t.points, beta.j.mod)

y.min <- min(c(overall.conf[[j]][2,],beta.j.val))
y.min <- -0.005
y.max <- max(c(overall.conf[[j]][1,], beta.j.val))
ylim <- c(y.min, y.max)

png(file = paste0("Bootstrap/novel/Boxplots-and-comparison/pres_reg_",coefs.names[j],".png"), width = 8000, height = 5000, units = "px", res=800)
par(mar=c(4.5, 4, 2.5, 1)+.1)
fda::fbplot(fit=B.list[,j,], x=t.points, method="MBD", color=pal.gb.fade, xlab='Period [s]',
            xlim=range(t.points), ylim=ylim, xaxt='n',cex.axis=1.8, cex.lab=1.8, ylab='')
title(main=names.coefs[j], font.main=1, cex.main=1.8)
lines(t.points, beta.j.val, lwd=4, col=pal[5], xaxt='n', yaxt='n')
lines(t.points, X.ita18[,j], col=pal[1], lwd=4)
points(t.points, X.ita18[,j], pch=19, col=pal[1])
abline(h=0, col=pal[5], lty=5, lwd=2)
axis(side=1, at=c(-2.4,xtick), labels = c(0,10^xtick), cex.axis=1.8)
#axis(side=2, cex.axis=2.5)
legend(-0.2, -0.0029, legend=c("Functional", "Scalar", "FBPlot"),
       col=c(pal[5],pal[1], pal.gb[2]), lwd=c(3,3, 3), lty=c(1,1,1), cex=1.8)

grid()
dev.off()

