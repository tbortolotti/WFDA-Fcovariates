setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(calculus)
library(ReconstPoFD)
library(tidyverse)
library("xtable")
library(psych)
library(progress)
library(fields)
library(ggplot2)

rm(list=ls())
graphics.off()
cat("\014")

## Utilities -------------------------------------------------------------------
load('Simulation/DATA/reg_info.RData')

B <- 100

beta0 <- reg.info$beta_estimates[[1]]$fd
beta0$coefs <- beta0$coefs/3 
beta1 <- reg.info$beta_estimates[[2]]$fd
beta2 <- reg.info$beta_estimates[[3]]$fd

grid0    <- reg.info$t.points
t.points <- seq(range(grid0)[1], range(grid0)[2], 0.25)

## PLOT RESULTS OF THE SIMULATION ----------------------------------------------

# # create a data frame
# variety=rep(LETTERS[1:7], each=40)
# treatment=rep(c("high","low"),each=20)
# note=seq(1:280)+sample(1:150, 280, replace=T)
# data=data.frame(variety, treatment, note)
# 
# # grouped boxplot
# ggplot(data, aes(x=variety, y=note, fill=treatment)) +
#   geom_boxplot()

## MSE: RECONSTRUCTION METHOD COMPARISON --------------------------------------------
# Devo creare un dataframe in cui ho:
levs <- 3 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("MSE", "Var", "Var+bias2"), each=B)
beta_info   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/Kraus-par3.RData')
# Beta 0
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_info[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))

j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
vars <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
beta_info[idxs] <- vars

j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_info[idxs] <- vars + rep(bias2, length(vars))

# Beta 1
i <- 2
j <- 1
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_info[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))

j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
vars <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
beta_info[idxs] <- vars

j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean.diff <- mean0 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_info[idxs] <- vars + rep(bias2, length(vars))

# Beta 2
i <- 3
j <- 1
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_info[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))

j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
vars <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
beta_info[idxs] <- vars

j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean.diff <- mean0 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_info[idxs] <- vars + rep(bias2, length(vars))

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_info)
method_order<- c("MSE", "Var", "Var+bias2")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_info, fill=method)) + 
  geom_boxplot() +
  #scale_y_continuous(limits=c(0,0.1)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  labs(fill = "Percentage of PO data")


