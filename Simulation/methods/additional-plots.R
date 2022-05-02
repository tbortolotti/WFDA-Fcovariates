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
library(latex2exp)

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

## BIAS: RECONSTRUCTION METHOD COMPARISON --------------------------------------
# Devo creare un dataframe in cui ho:
levs <- 6 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("Kraus:wgt", "KL-PC:wgt", "KL-AL:wgt", "Kraus", "KL-PC", "KL-AL"),
                   each=B)
beta_bias   <- numeric(3*levs)

# Riempio il beta_bias
load('Simulation/Results/repeated-simulations/Kraus-par3.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2


i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

load('Simulation/Results/repeated-simulations-new/KLNoAl-par3.RData')
i <- 1 #index of the coefficient
j <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

load('Simulation/Results/repeated-simulations/Kraus_nowgts.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

load('Simulation/Results/repeated-simulations/KLNoAl_nowgts.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

## BIAS: PO VARIATION --------------------------------------
## Pesato ----
# Devo creare un dataframe in cui ho:
levs <- 4 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs)
method      <- rep(c("10%",  "25%",  "40%",  "70%"),
                   each=1)
beta_bias   <- numeric(3*levs)

load('Simulation/Results/repeated-simulations-new/KLAl-perc1.RData')
i <- 1 #index of the coefficient
j <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

load('Simulation/Results/repeated-simulations-new/KLAl-perc2_5.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2


load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_bias)
method_order<- c("10%", "25%", "40%", "70%")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

# grouped plot
ggplot(data.box, aes(x=coefficient, y=beta_bias, colour=method)) + 
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  geom_point(aes(size = 5))+
  scale_y_continuous(limits=c(0,0.02)) +
  ylab(TeX(r'(bias$^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "PO percentage")

## Non pesato ----
# Devo creare un dataframe in cui ho:
levs <- 4 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs)
method      <- rep(c("10%", "25%", "40%",  "70%"),
                   each=1)
beta_bias   <- numeric(3*levs)

load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2


i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2


load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2


load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean.diff <- mean0 - beta0
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 2
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean1 <- mean.fd(beta1.est)
mean.diff <- mean1 - beta1
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

i <- 3
idx.inf <- (levs*(i-1) + j-1) + 1
idx.sup <- (levs*(i-1)+j)
idxs <- idx.inf:idx.sup
mean2 <- mean.fd(beta2.est)
mean.diff <- mean2 - beta2
bias2 <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
beta_bias[idxs] <- bias2

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_bias)
method_order<- c("10%", "25%", "40%", "70%")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

# grouped plot
ggplot(data.box, aes(x=coefficient, y=beta_bias, colour=method)) + 
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  geom_point(aes(size = 5))+
  #scale_y_continuous(limits=c(0,0.06)) +
  ylab(TeX(r'(bias$^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "PO percentage")


## MSE [0, 1.75] --------------------------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations-new/KLAl-perc1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations-new/KLAl-perc2_5.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
i <- 1
j <- 7
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
i <- 1
j <- 8
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.05)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[0 s, 1.75 s]")


## MSE [1.75, 3.5] --------------------------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-perc1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-perc2_5.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
i <- 1
j <- 7
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
i <- 1
j <- 8
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.05)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[1.75 s, 3.5 s]")


## VAR [0, 1.75] --------------------------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
i <- 1
j <- 1
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations-new/KLAl-perc1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations-new/KLAl-perc2_5.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
i <- 1
j <- 7
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
i <- 1
j <- 8
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.75)))

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.05)) +
  ylab(TeX(r'($||\hat{\beta}^b - \sum_{b=1}^B \hat{\beta^b} ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[0, 1.75 s]")


## VAR [1.75, 3.5] ----
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
i <- 1
j <- 1
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-perc1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-perc2_5.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
i <- 1
j <- 7
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-perc7.RData')
i <- 1
j <- 8
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.75,3.5)))

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.05)) +
  ylab(TeX(r'($||\hat{\beta}^b - \sum_{b=1}^B \hat{\beta^b} ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[1.75 s, 3.5 s]")

## NEW MSE FOR WEIGHTS DEFINITION ----------------------------------------------
## [0,1.5] -----------------------------
levs <- 6
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par2.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par5.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-parinf.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.5)))


# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.025)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[0,1.5 s]")


## [1.5, 2.5] --------------------------
levs <- 6
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par2.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par5.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-parinf.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.5,2.5)))


# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.025)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[1.5 s, 2.5 s]")

## [2.5,3.5] -------------------------------
levs <- 6
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par2.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par5.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-parinf.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(2.5,3.5)))


# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.025)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[2.5 s, 3.5 s]")

## [0, 1.75] --------------------------
levs <- 6
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations-new/KLAl-par2.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations-new/KLAl-par5.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

load('Simulation/Results/repeated-simulations-new/KLAl-parinf.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,1.75)))


# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.05)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[0, 1.75 s]")

## [1.75, 3.5] --------------------------
levs <- 6
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par2.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par5.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-parinf.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(1.75,3.5)))


# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.05)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "[1.75 s, 3.5 s]")


## MSE  COMPARISON --------------------------------------------
levs <- 6
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf"),
                   each=B)
beta_mse   <- numeric(levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
j <- 1
idx.inf <- (j-1)*B + 1
idx.sup <- j*B
idxs <- idx.inf:idx.sup
beta_mse[idxs] <- MSE*3.5

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
j <- 2
idx.inf <- (j-1)*B + 1
idx.sup <- j*B
idxs <- idx.inf:idx.sup
beta_mse[idxs] <- MSE*3.5

load('Simulation/Results/repeated-simulations-new/KLAl-par2.RData')
j <- 3
idx.inf <- (j-1)*B + 1
idx.sup <- j*B
idxs <- idx.inf:idx.sup
beta_mse[idxs] <- MSE*3.5

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
j <- 4
idx.inf <- (j-1)*B + 1
idx.sup <- j*B
idxs <- idx.inf:idx.sup
beta_mse[idxs] <- MSE*3.5

load('Simulation/Results/repeated-simulations-new/KLAl-par5.RData')
j <- 5
idx.inf <- (j-1)*B + 1
idx.sup <- j*B
idxs <- idx.inf:idx.sup
beta_mse[idxs] <- MSE*3.5

load('Simulation/Results/repeated-simulations-new/KLAl-parinf.RData')
j <- 6
idx.inf <- (j-1)*B + 1
idx.sup <- j*B
idxs <- idx.inf:idx.sup
beta_mse[idxs] <- MSE*3.5

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.025)) +
  ylab(TeX(r'($||\hat{y}^b - y ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "Weights parameter")



## PROVA ----------------------------------------------
levs <- 6
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1 #index of the coefficient
j <- 1 #index of the method
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par2.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-par5.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

load('Simulation/Results/repeated-simulations-new/KLAl-parinf.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta0.est - beta0
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta1.est - beta1
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
beta.res <- beta2.est - beta2
beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=c(0,3.5)))


# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  ylab(TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "Weights parameter")


# VARIANCE
levs <- 6
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
i <- 1
j <- 1
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
i <- 1
j <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations-new/KLAl-par2.RData')
i <- 1
j <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations-new/KLAl-par3.RData')
i <- 1
j <- 4
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations-new/KLAl-par5.RData')
i <- 1
j <- 5
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations-new/KLAl-parinf.RData')
i <- 1
j <- 6
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 2
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

i <- 3
idx.inf <- (levs*(i-1) + j-1)*B + 1
idx.sup <- (levs*(i-1)+j)*B
idxs <- idx.inf:idx.sup
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("No weights", "a=1", "a=2", "a=3", "a=5", "a=inf")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  ylab(TeX(r'($||\hat{\beta}^b - \sum_{b=1}^B \hat{\beta^b} ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  labs(fill = "Weights parameter")


