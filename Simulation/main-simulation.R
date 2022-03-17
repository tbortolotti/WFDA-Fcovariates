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

## Load Functions --------------------------------------------------------------
source('Simulation/methods/workflow_weighted_analysis.R')
source('Simulation/methods/generate_data.R')
load('Simulation/DATA/reg_info.RData')

## Simulation ------------------------------------------------------------------
fix.par <- 3
perc.po <- 4
smooth.noise <- 15

## Set the case information
case.info <- list(n.sim       = 101,
                  noise       = 0.01*smooth.noise, #all analysis with noise=0.1
                  ext.noise   = 0.5,
                  perc        = 0.1*perc.po,
                  left.bound  = 1.5,
                  right.bound = 3.5)

B <- 100
MSE<- numeric(B)

#method <- 'Kraus'
#method <- 'KLNoAl'
method <- 'KLAl'

(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  
  ## Simulate data
  seed <- 140996
  simulated_data <- generate_data(seed      = (seed+b),
                                  case      = "CASE-1",
                                  reg.info  = reg.info,
                                  case.info = case.info)
  
  ## Input data for WFDA
  t.points       <- simulated_data$t.points
  T_hp           <- simulated_data$T_hp
  curves         <- simulated_data$curves
  curves.true.fd <- simulated_data$curves.true.fd
  xlist          <- simulated_data$xlist
  
  ## WFDA
  method_evaluation <- workflow_weighted_analysis(b              = 1,
                                                  B              = case.info$n.sim,
                                                  t.points       = t.points,
                                                  breaks         = t.points,
                                                  T_hp           = T_hp,
                                                  curves         = curves,
                                                  curves.true.fd = curves.true.fd,
                                                  xlist          = xlist,
                                                  method         = method,
                                                  fix.par        = fix.par,
                                                  wgts.flag      = FALSE)
  
  MSE[b] <- method_evaluation$MSE
  
  if(b==1)
  {
    beta0.est <- method_evaluation$beta_estimates[[1]]$fd
    beta1.est <- method_evaluation$beta_estimates[[2]]$fd
    beta2.est <- method_evaluation$beta_estimates[[3]]$fd
    
    beta0.est$coefs <- beta0.est$coefs %*% rep(1,B)
    beta1.est$coefs <- beta1.est$coefs %*% rep(1,B)
    beta2.est$coefs <- beta2.est$coefs %*% rep(1,B)
  } else {
    beta0.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[1]]$fd$coefs, ncol=1)
    beta1.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[2]]$fd$coefs, ncol=1)
    beta2.est$coefs[,b] <- as.matrix(method_evaluation$beta_estimates[[3]]$fd$coefs, ncol=1)
  }
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

# name.file <- paste0('Simulation/Results/repeated-simulations-smootherr/',method,'-se',smooth.noise,'.RData')
# save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)

name.file <- paste0('Simulation/Results/repeated-simulations-smootherr/',method,'-se',smooth.noise,'_nowgts.RData')
save(MSE, beta0.est, beta1.est, beta2.est, file=name.file)

# RECONSTRUCTION METHOD COMPARISON ---------------------------------------------
beta0 <- reg.info$beta_estimates[[1]]$fd
beta0$coefs <- beta0$coefs/3 
beta1 <- reg.info$beta_estimates[[2]]$fd
beta2 <- reg.info$beta_estimates[[3]]$fd

grid0    <- reg.info$t.points
t.points <- seq(range(grid0)[1], range(grid0)[2], 0.25)

load('Simulation/Results/repeated-simulations/Kraus-par3.RData')
Kraus <- MSE
beta.res <- beta0.est - beta0
Kraus.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
Kraus.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))


load('Simulation/Results/repeated-simulations/KLNoAl-par3.RData')
KLNoAl <- MSE
beta.res <- beta0.est - beta0
KLNoAl.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLNoAl.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
KLAl <- MSE
beta.res <- beta0.est - beta0
KLAl.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/Kraus_nowgts.RData')
Kraus_nowgt <- MSE
beta.res <- beta0.est - beta0
Kraus_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
Kraus_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLNoAl_nowgts.RData')
KLNoAl_nowgt <- MSE
beta.res <- beta0.est - beta0
KLNoAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLNoAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
KLAl_nowgt <- MSE
beta.res <- beta0.est - beta0
KLAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

## MSE -------------------------------------------------------------------------
class <- c(rep('Kraus', B),
           rep('KLNoAl', B),
           rep('KLAl', B),
           rep('Kraus - unwgt', B),
           rep('KLNoAl - unwgt', B),
           rep('KLAl - unwgt', B))

data.box <- data.frame(c(Kraus,
                         KLNoAl,
                         KLAl,
                         Kraus_nowgt,
                         KLNoAl_nowgt,
                         KLAl_nowgt), as.factor(class))

names(data.box) <- c('MSE', 'class')
class_order<- c("Kraus", "KLNoAl", "KLAl", "Kraus - unwgt", "KLNoAl - unwgt", "KLAl - unwgt")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = MSE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.02)) +
  labs(y="Mean Squared Error", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")

## Beta 0 ----------------------------------------------------------------------
## BIAS BETA0
class <- c(rep('Kraus', B),
           rep('KLNoAl', B),
           rep('KLAl', B),
           rep('Kraus - unwgt', B),
           rep('KLNoAl - unwgt', B),
           rep('KLAl - unwgt', B))

data.box <- data.frame(c(Kraus.bias,
                         KLNoAl.bias,
                         KLAl.bias,
                         Kraus_nowgt.bias,
                         KLNoAl_nowgt.bias,
                         KLAl_nowgt.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c("Kraus", "KLNoAl", "KLAl", "Kraus - unwgt", "KLNoAl - unwgt", "KLAl - unwgt")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  #scale_y_continuous(limits=c(0.08,0.2)) +
  labs(y="Intercept BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA0
class <- c(rep('Kraus', B),
           rep('KLNoAl', B),
           rep('KLAl', B),
           rep('Kraus - unwgt', B),
           rep('KLNoAl - unwgt', B),
           rep('KLAl - unwgt', B))

data.box <- data.frame(c(Kraus.var,
                         KLNoAl.var,
                         KLAl.var,
                         Kraus_nowgt.var,
                         KLNoAl_nowgt.var,
                         KLAl_nowgt.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c("Kraus", "KLNoAl", "KLAl", "Kraus - unwgt", "KLNoAl - unwgt", "KLAl - unwgt")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  #scale_y_continuous(limits=c(0.08,0.2)) +
  labs(y="Intercept VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## Beta 1 ----------------------------------------------------------------------

beta1 <- reg.info$beta_estimates[[2]]$fd

load('Simulation/Results/repeated-simulations/Kraus-par3.RData')
beta.res <- beta1.est - beta1
Kraus.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
Kraus.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))


load('Simulation/Results/repeated-simulations/KLNoAl-par3.RData')
beta.res <- beta1.est - beta1
KLNoAl.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLNoAl.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
beta.res <- beta1.est - beta1
KLAl.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/Kraus_nowgts.RData')
beta.res <- beta1.est - beta1
Kraus_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
Kraus_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLNoAl_nowgts.RData')
beta.res <- beta1.est - beta1
KLNoAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLNoAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
beta.res <- beta1.est - beta1
KLAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))


## BIAS BETA 1
class <- c(rep('Kraus', B),
           rep('KLNoAl', B),
           rep('KLAl', B),
           rep('Kraus - unwgt', B),
           rep('KLNoAl - unwgt', B),
           rep('KLAl - unwgt', B))

data.box <- data.frame(c(Kraus.bias,
                         KLNoAl.bias,
                         KLAl.bias,
                         Kraus_nowgt.bias,
                         KLNoAl_nowgt.bias,
                         KLAl_nowgt.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c("Kraus", "KLNoAl", "KLAl", "Kraus - unwgt", "KLNoAl - unwgt", "KLAl - unwgt")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.1)) +
  labs(y="Beta 1 BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA 1
class <- c(rep('Kraus', B),
           rep('KLNoAl', B),
           rep('KLAl', B),
           rep('Kraus - unwgt', B),
           rep('KLNoAl - unwgt', B),
           rep('KLAl - unwgt', B))

data.box <- data.frame(c(Kraus.var,
                         KLNoAl.var,
                         KLAl.var,
                         Kraus_nowgt.var,
                         KLNoAl_nowgt.var,
                         KLAl_nowgt.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c("Kraus", "KLNoAl", "KLAl", "Kraus - unwgt", "KLNoAl - unwgt", "KLAl - unwgt")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  #scale_y_continuous(limits=c(0.08,0.2)) +
  labs(y="Beta 1 VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## Beta 2 ----------------------------------------------------------------------

beta2 <- reg.info$beta_estimates[[3]]$fd

load('Simulation/Results/repeated-simulations/Kraus-par3.RData')
beta.res <- beta2.est - beta2
Kraus.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
Kraus.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))


load('Simulation/Results/repeated-simulations/KLNoAl-par3.RData')
beta.res <- beta2.est - beta2
KLNoAl.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLNoAl.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
beta.res <- beta2.est - beta2
KLAl.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/Kraus_nowgts.RData')
beta.res <- beta2.est - beta2
Kraus_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
Kraus_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLNoAl_nowgts.RData')
beta.res <- beta2.est - beta2
KLNoAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLNoAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
beta.res <- beta2.est - beta2
KLAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))


## BIAS BETA 2
class <- c(rep('Kraus', B),
           rep('KLNoAl', B),
           rep('KLAl', B),
           rep('Kraus - unwgt', B),
           rep('KLNoAl - unwgt', B),
           rep('KLAl - unwgt', B))

data.box <- data.frame(c(Kraus.bias,
                         KLNoAl.bias,
                         KLAl.bias,
                         Kraus_nowgt.bias,
                         KLNoAl_nowgt.bias,
                         KLAl_nowgt.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c("Kraus", "KLNoAl", "KLAl", "Kraus - unwgt", "KLNoAl - unwgt", "KLAl - unwgt")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.05)) +
  labs(y="Beta 2 BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA 2
class <- c(rep('Kraus', B),
           rep('KLNoAl', B),
           rep('KLAl', B),
           rep('Kraus - unwgt', B),
           rep('KLNoAl - unwgt', B),
           rep('KLAl - unwgt', B))

data.box <- data.frame(c(Kraus.var,
                         KLNoAl.var,
                         KLAl.var,
                         Kraus_nowgt.var,
                         KLNoAl_nowgt.var,
                         KLAl_nowgt.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c("Kraus", "KLNoAl", "KLAl", "Kraus - unwgt", "KLNoAl - unwgt", "KLAl - unwgt")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  #scale_y_continuous(limits=c(0,0.005)) +
  labs(y="Beta 2 VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


# PARAMETER COMPARISON ------------------------------------------------
beta0 <- reg.info$beta_estimates[[1]]$fd
beta0$coefs <- beta0$coefs/3
beta1 <- reg.info$beta_estimates[[2]]$fd
beta2 <- reg.info$beta_estimates[[3]]$fd

grid0    <- reg.info$t.points
t.points <- seq(range(grid0)[1], range(grid0)[2], 0.25)

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
KLAl_nowgt <- MSE
beta.res <- beta0.est - beta0
KLAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
KLAl_par1 <- MSE
beta.res <- beta0.est - beta0
KLAl_par1.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_par1.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par2.RData')
KLAl_par2 <- MSE
beta.res <- beta0.est - beta0
KLAl_par2.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_par2.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
KLAl_par3 <- MSE
beta.res <- beta0.est - beta0
KLAl_par3.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_par3.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par5.RData')
KLAl_par5 <- MSE
beta.res <- beta0.est - beta0
KLAl_par5.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_par5.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

## MSE -------------------------------------------------------------------------
class <- c(rep('KLAl - nowgt', B),
           rep('KLAl - par1', B),
           rep('KLAl - par2', B),
           rep('KLAl - par3', B),
           rep('KLAl - par5', B))

data.box <- data.frame(c(KLAl_nowgt,
                         KLAl_par1,
                         KLAl_par2,
                         KLAl_par3,
                         KLAl_par5), as.factor(class))

names(data.box) <- c('MSE', 'class')
class_order<- c("KLAl - nowgt", "KLAl - par1", "KLAl - par2", "KLAl - par3", "KLAl - par5")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = MSE))+
  geom_boxplot(color="black", fill=tim.colors(5))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.01)) +
  labs(y="Mean Squared Error", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")

## Beta 0 ----------------------------------------------------------------------
## BIAS BETA0
class <- c(rep('KLAl - nowgt', B),
           rep('KLAl - par1', B),
           rep('KLAl - par2', B),
           rep('KLAl - par3', B),
           rep('KLAl - par5', B))

data.box <- data.frame(c(KLAl_nowgt.bias,
                         KLAl_par1.bias,
                         KLAl_par2.bias,
                         KLAl_par3.bias,
                         KLAl_par5.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c("KLAl - nowgt", "KLAl - par1", "KLAl - par2", "KLAl - par3", "KLAl - par5")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(5))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.1)) +
  labs(y="Intercept BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA0
class <- c(rep('KLAl - nowgt', B),
           rep('KLAl - par1', B),
           rep('KLAl - par2', B),
           rep('KLAl - par3', B),
           rep('KLAl - par5', B))

data.box <- data.frame(c(KLAl_nowgt.var,
                         KLAl_par1.var,
                         KLAl_par2.var,
                         KLAl_par3.var,
                         KLAl_par5.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c("KLAl - nowgt", "KLAl - par1", "KLAl - par2", "KLAl - par3", "KLAl - par5")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(5))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.075)) +
  labs(y="Intercept VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## Beta 1 ----------------------------------------------------------------------

beta1 <- reg.info$beta_estimates[[2]]$fd

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
beta.res <- beta1.est - beta1
KLAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
beta.res <- beta1.est - beta1
KLAl_par1.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_par1.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par2.RData')
beta.res <- beta1.est - beta1
KLAl_par2.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_par2.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
beta.res <- beta1.est - beta1
KLAl_par3.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_par3.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par5.RData')
beta.res <- beta1.est - beta1
KLAl_par5.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_par5.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

## BIAS BETA1
class <- c(rep('KLAl - nowgt', B),
           rep('KLAl - par1', B),
           rep('KLAl - par2', B),
           rep('KLAl - par3', B),
           rep('KLAl - par5', B))

data.box <- data.frame(c(KLAl_nowgt.bias,
                         KLAl_par1.bias,
                         KLAl_par2.bias,
                         KLAl_par3.bias,
                         KLAl_par5.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c("KLAl - nowgt", "KLAl - par1", "KLAl - par2", "KLAl - par3","KLAl - par5")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(5))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.1)) +
  labs(y="BETA 1 BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA1
class <- c(rep('KLAl - nowgt', B),
           rep('KLAl - par1', B),
           rep('KLAl - par2', B),
           rep('KLAl - par3', B),
           rep('KLAl - par5', B))

data.box <- data.frame(c(KLAl_nowgt.var,
                         KLAl_par1.var,
                         KLAl_par2.var,
                         KLAl_par3.var,
                         KLAl_par5.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c("KLAl - nowgt", "KLAl - par1", "KLAl - par2", "KLAl - par3", "KLAl - par5")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(5))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.025)) +
  labs(y="BETA 1 VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## Beta 2 ----------------------------------------------------------------------

beta2 <- reg.info$beta_estimates[[3]]$fd

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
beta.res <- beta2.est - beta2
KLAl_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par1.RData')
beta.res <- beta2.est - beta2
KLAl_par1.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_par1.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par2.RData')
beta.res <- beta2.est - beta2
KLAl_par2.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_par2.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
beta.res <- beta2.est - beta2
KLAl_par3.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_par3.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par5.RData')
beta.res <- beta2.est - beta2
KLAl_par5.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_par5.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))


## BIAS BETA2
class <- c(rep('KLAl - nowgt', B),
           rep('KLAl - par1', B),
           rep('KLAl - par2', B),
           rep('KLAl - par3', B),
           rep('KLAl - par5', B))

data.box <- data.frame(c(KLAl_nowgt.bias,
                         KLAl_par1.bias,
                         KLAl_par2.bias,
                         KLAl_par3.bias,
                         KLAl_par5.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c("KLAl - nowgt", "KLAl - par1", "KLAl - par2", "KLAl - par3","KLAl - par5")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(5))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.04)) +
  labs(y="BETA 2 BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA2
class <- c(rep('KLAl - nowgt', B),
           rep('KLAl - par1', B),
           rep('KLAl - par2', B),
           rep('KLAl - par3', B),
           rep('KLAl - par5', B))

data.box <- data.frame(c(KLAl_nowgt.var,
                         KLAl_par1.var,
                         KLAl_par2.var,
                         KLAl_par3.var,
                         KLAl_par5.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c("KLAl - nowgt", "KLAl - par1", "KLAl - par2", "KLAl - par3", "KLAl - par5")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(5))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.03)) +
  labs(y="BETA 2 VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


# PERCENTAGE OF PO DATA IMPACT -------------------------------------------------
beta0 <- reg.info$beta_estimates[[1]]$fd
beta0$coefs <- beta0$coefs/3
beta1 <- reg.info$beta_estimates[[2]]$fd
beta2 <- reg.info$beta_estimates[[3]]$fd

grid0    <- reg.info$t.points
t.points <- seq(range(grid0)[1], range(grid0)[2], 0.25)

load('Simulation/Results/repeated-simulations/KLAl-perc1.RData')
KLAl_perc1 <- MSE
beta.res <- beta0.est - beta0
KLAl_perc1.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_perc1.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5.RData')
KLAl_perc25 <- MSE
beta.res <- beta0.est - beta0
KLAl_perc25.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_perc25.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
KLAl_perc4 <- MSE
beta.res <- beta0.est - beta0
KLAl_perc4.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_perc4.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
KLAl_perc1_nowgt <- MSE
beta.res <- beta0.est - beta0
KLAl_perc1_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_perc1_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
KLAl_perc25_nowgt <- MSE
beta.res <- beta0.est - beta0
KLAl_perc25_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_perc25_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
KLAl_perc4_nowgt <- MSE
beta.res <- beta0.est - beta0
KLAl_perc4_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta0.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta0.est - mean0
KLAl_perc4_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

## MSE -------------------------------------------------------------------------
class <- c(rep('KLAl-perc1', B),
           rep('KLAl-perc2.5', B),
           rep('KLAl-perc4', B),
           rep('KLAl-perc1-unwgt', B),
           rep('KLAl-perc2.5-unwgt', B),
           rep('KLAl-perc4-unwgt', B))

data.box <- data.frame(c(KLAl_perc1,
                         KLAl_perc25,
                         KLAl_perc4,
                         KLAl_perc1_nowgt,
                         KLAl_perc25_nowgt,
                         KLAl_perc4_nowgt), as.factor(class))

names(data.box) <- c('MSE', 'class')
class_order<- c('KLAl-perc1', 'KLAl-perc2.5', 'KLAl-perc4', 'KLAl-perc1-unwgt', 'KLAl-perc2.5-unwgt','KLAl-perc4-unwgt')
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = MSE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.01)) +
  labs(y="Mean Squared Error", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")

## Beta 0 ----------------------------------------------------------------------
## BIAS BETA0
class <- c(rep('KLAl-perc1', B),
           rep('KLAl-perc2.5', B),
           rep('KLAl-perc4', B),
           rep('KLAl-perc1-unwgt', B),
           rep('KLAl-perc2.5-unwgt', B),
           rep('KLAl-perc4-unwgt', B))

data.box <- data.frame(c(KLAl_perc1.bias,
                         KLAl_perc25.bias,
                         KLAl_perc4.bias,
                         KLAl_perc1_nowgt.bias,
                         KLAl_perc25_nowgt.bias,
                         KLAl_perc4_nowgt.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c('KLAl-perc1', 'KLAl-perc2.5', 'KLAl-perc4', 'KLAl-perc1-unwgt', 'KLAl-perc2.5-unwgt','KLAl-perc4-unwgt')
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.1)) +
  labs(y="Intercept BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA0
class <- c(rep('KLAl-perc1', B),
           rep('KLAl-perc2.5', B),
           rep('KLAl-perc4', B),
           rep('KLAl-perc1-unwgt', B),
           rep('KLAl-perc2.5-unwgt', B),
           rep('KLAl-perc4-unwgt', B))

data.box <- data.frame(c(KLAl_perc1.var,
                         KLAl_perc25.var,
                         KLAl_perc4.var,
                         KLAl_perc1_nowgt.var,
                         KLAl_perc25_nowgt.var,
                         KLAl_perc4_nowgt.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c('KLAl-perc1', 'KLAl-perc2.5', 'KLAl-perc4', 'KLAl-perc1-unwgt', 'KLAl-perc2.5-unwgt','KLAl-perc4-unwgt')
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.075)) +
  labs(y="Intercept VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## Beta 1 ----------------------------------------------------------------------

beta1 <- reg.info$beta_estimates[[2]]$fd

load('Simulation/Results/repeated-simulations/KLAl-perc1.RData')
beta.res <- beta1.est - beta1
KLAl_perc1.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_perc1.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5.RData')
beta.res <- beta1.est - beta1
KLAl_perc25.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_perc25.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
beta.res <- beta1.est - beta1
KLAl_perc4.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_perc4.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
beta.res <- beta1.est - beta1
KLAl_perc1_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_perc1_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
beta.res <- beta1.est - beta1
KLAl_perc25_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_perc25_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
beta.res <- beta1.est - beta1
KLAl_perc4_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta1.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta1.est - mean0
KLAl_perc4_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

## BIAS BETA 1
class <- c(rep('KLAl-perc1', B),
           rep('KLAl-perc2.5', B),
           rep('KLAl-perc4', B),
           rep('KLAl-perc1-unwgt', B),
           rep('KLAl-perc2.5-unwgt', B),
           rep('KLAl-perc4-unwgt', B))

data.box <- data.frame(c(KLAl_perc1.bias,
                         KLAl_perc25.bias,
                         KLAl_perc4.bias,
                         KLAl_perc1_nowgt.bias,
                         KLAl_perc25_nowgt.bias,
                         KLAl_perc4_nowgt.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c('KLAl-perc1', 'KLAl-perc2.5', 'KLAl-perc4', 'KLAl-perc1-unwgt', 'KLAl-perc2.5-unwgt','KLAl-perc4-unwgt')
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.075)) +
  labs(y="BETA 1 BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA 1
class <- c(rep('KLAl-perc1', B),
           rep('KLAl-perc2.5', B),
           rep('KLAl-perc4', B),
           rep('KLAl-perc1-unwgt', B),
           rep('KLAl-perc2.5-unwgt', B),
           rep('KLAl-perc4-unwgt', B))

data.box <- data.frame(c(KLAl_perc1.var,
                         KLAl_perc25.var,
                         KLAl_perc4.var,
                         KLAl_perc1_nowgt.var,
                         KLAl_perc25_nowgt.var,
                         KLAl_perc4_nowgt.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c('KLAl-perc1', 'KLAl-perc2.5', 'KLAl-perc4', 'KLAl-perc1-unwgt', 'KLAl-perc2.5-unwgt','KLAl-perc4-unwgt')
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  #scale_y_continuous(limits=c(0,0.075)) +
  labs(y="BETA 1 VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")

## Beta 2 ----------------------------------------------------------------------

beta2 <- reg.info$beta_estimates[[3]]$fd

load('Simulation/Results/repeated-simulations/KLAl-perc1.RData')
beta.res <- beta2.est - beta2
KLAl_perc1.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_perc1.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5.RData')
beta.res <- beta2.est - beta2
KLAl_perc25.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_perc25.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-par3.RData')
beta.res <- beta2.est - beta2
KLAl_perc4.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_perc4.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
beta.res <- beta2.est - beta2
KLAl_perc1_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_perc1_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
beta.res <- beta2.est - beta2
KLAl_perc25_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_perc25_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
beta.res <- beta2.est - beta2
KLAl_perc4_nowgt.bias <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
mean0 <- mean.fd(beta2.est)
mean0$coefs <- mean0$coefs%*%rep(1,B)
beta.diff <- beta2.est - mean0
KLAl_perc4_nowgt.var <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))

## BIAS BETA2
class <- c(rep('KLAl-perc1', B),
           rep('KLAl-perc2.5', B),
           rep('KLAl-perc4', B),
           rep('KLAl-perc1-unwgt', B),
           rep('KLAl-perc2.5-unwgt', B),
           rep('KLAl-perc4-unwgt', B))

data.box <- data.frame(c(KLAl_perc1.bias,
                         KLAl_perc25.bias,
                         KLAl_perc4.bias,
                         KLAl_perc1_nowgt.bias,
                         KLAl_perc25_nowgt.bias,
                         KLAl_perc4_nowgt.bias), as.factor(class))

names(data.box) <- c('BIAS', 'class')
class_order<- c('KLAl-perc1', 'KLAl-perc2.5', 'KLAl-perc4', 'KLAl-perc1-unwgt', 'KLAl-perc2.5-unwgt','KLAl-perc4-unwgt')
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = BIAS))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.04)) +
  labs(y="BETA 2 BIAS", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


## VARIANCE BETA 2
class <- c(rep('KLAl-perc1', B),
           rep('KLAl-perc2.5', B),
           rep('KLAl-perc4', B),
           rep('KLAl-perc1-unwgt', B),
           rep('KLAl-perc2.5-unwgt', B),
           rep('KLAl-perc4-unwgt', B))

data.box <- data.frame(c(KLAl_perc1.var,
                         KLAl_perc25.var,
                         KLAl_perc4.var,
                         KLAl_perc1_nowgt.var,
                         KLAl_perc25_nowgt.var,
                         KLAl_perc4_nowgt.var), as.factor(class))

names(data.box) <- c('VARIANCE', 'class')
class_order<- c('KLAl-perc1', 'KLAl-perc2.5', 'KLAl-perc4', 'KLAl-perc1-unwgt', 'KLAl-perc2.5-unwgt','KLAl-perc4-unwgt')
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = VARIANCE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.035)) +
  labs(y="BETA 2 VARIANCE", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")

