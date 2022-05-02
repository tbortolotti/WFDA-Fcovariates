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
library(cowplot)
library(gridExtra)
library(ggpubr)

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

box.dir <- '/Users/teresabortolotti/Documents/R/WFDA-Fcovariates/Simulation/Boxplots-oldwgts'

## MSE: WEIGHTS DEFINITION ----------------------------------------------------------
levs <- 7
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights","a=5","a=10","a=15", "a=20", "a=inf", "0-wgts"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
{
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1 #index of the coefficient
  j <- 1 #index of the method
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par5.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par15.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par20.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par100.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-0wgts.RData')
  i <- 1
  j <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
}


# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("No weights","a=5","a=10","a=15" ,"a=20", "a=inf", "0-wgts")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
          geom_boxplot() 
pgplot <- pgplot +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() + 
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)'), fill = "Parameter: ", title="(a) Weights definition: MSE") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        legend.position="bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))
pgplot + 
  theme(legend.position="none")
ggsave(filename = "weights-parameter.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL,)

pg_legend <- cowplot::get_legend(pgplot)
as_ggplot(pg_legend)


ggsave(filename = "legend-weights-parameter.pdf",
  plot = last_plot(),
  device = NULL,
  path = box.dir,
  scale = 1,
  width = 12,
  height = 0.7,
  units = "in",
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,)

## VARIANCE: WEIGHTS DEFINITION ------------------------------------------------
# Devo creare un dataframe in cui ho:
levs <- 7
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("No weights", "a=5", "a=10", "a=15", "a=20", "a=inf", "0-wgts"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
{
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par5.RData')
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par15.RData')
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par20.RData')
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par100.RData')
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
  
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-0wgts.RData')
  i <- 1
  j <- 7
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
  
}

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("No weights", "a=5", "a=10", "a=15", "a=20", "a=inf", "0-wgts")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pg_plot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "Parameter", title="(b) Weights definition: Variance")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))
pg_plot + theme(legend.position = "none")

ggsave(filename = "weights-parameter-var.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)


## MSE: RECONSTRUCTION METHOD COMPARISON ---------------------------------------
# Devo creare un dataframe in cui ho:
levs <- 6 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("Kraus:wgt", "KL-PC:wgt", "KL-AL:wgt", "Kraus", "KL-PC", "KL-AL"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
{
  load('Simulation/Results/repeated-simulations-oldwgts/Kraus.RData')
  i <- 1 #index of the coefficient
  j <- 1 #index of the method
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLNoAl.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations/Kraus_nowgts.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations/KLNoAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
}

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order  <- c("Kraus:wgt", "KL-PC:wgt", "KL-AL:wgt", "Kraus", "KL-PC", "KL-AL")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)'),
       fill = "Method: ", title="(a) Reconstruction methods: MSE")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        legend.position="bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

pgplot + 
  theme(legend.position="none")
ggsave(filename = "Reconstruction-comparison.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)

pg_legend <- cowplot::get_legend(pgplot)
as_ggplot(pg_legend)


ggsave(filename = "legend-rec-method.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 12,
       height = 0.7,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL,)

## VARIANCE: RECONSTRUCTION METHOD COMPARISON ----------------------------------
levs <- 6 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("Kraus: wgt", "KL-PC: wgt", "KL-AL: wgt", "Kraus", "KL-PC", "KL-AL"),
                   each=B)
beta_var   <- numeric(3*levs*B)
bias0.vec   <- numeric(levs)
bias1.vec   <- numeric(levs)
bias2.vec   <- numeric(levs)

# Riempio il beta_var
{
  load('Simulation/Results/repeated-simulations-oldwgts/Kraus.RData')
  i <- 1
  j <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLNoAl.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations/Kraus_nowgts.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations/KLNoAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias0.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias1.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2.vec[j] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
}

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("Kraus: wgt", "KL-PC: wgt", "KL-AL: wgt", "Kraus", "KL-PC", "KL-AL")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "Method", title="(b) Reconstruction methods: Variance")+
theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
      axis.text.x = element_text(size = 22),
      axis.title.x = element_text(size = 22),
      axis.text.y = element_text(size = 22),
      axis.title.y = element_text(size = 22),
      legend.title = element_text(size = 22),
      legend.text = element_text(size=22))
pgplot + theme(legend.position = "none")

ggsave(filename = "Reconstruction-comparison-var.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)


## MSE: VARIATION OF THE PERCENTAGE OF PARTIALLY OBSERVES CURVES --------------------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_mse   <- numeric(3*levs*B)

# Riempio il beta_mse
{
  load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
  i <- 1 #index of the coefficient
  j <- 1 #index of the method
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc1.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc2_5.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
  i <- 1
  j <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc7.RData')
  i <- 1
  j <- 8
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta0.est - beta0
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta1.est - beta1
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  beta.res <- beta2.est - beta2
  beta_mse[idxs] <- diag(inprod(beta.res, beta.res, 0, 0, rng=range(t.points)))
  
}
# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_mse)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_mse, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \beta ||_2^2$)'),
       fill = "", title="(a) Varying fraction of PO data: MSE") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22),
        legend.position="bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

pgplot + 
  theme(legend.position="none")
ggsave(filename = "PO-percentage.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)

pg_legend <- cowplot::get_legend(pgplot)
as_ggplot(pg_legend)


ggsave(filename = "legend-PO.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 12,
       height = 0.7,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL,)



## VARIANCE: VARIATION OF THE PERCENTAGE OF PARTIALLY OBSERVES CURVES ----------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
{
  load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc1.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc2_5.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
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
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
  i <- 1
  j <- 7
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
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc7.RData')
  i <- 1
  j <- 8
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
  
}

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "PO percentage", title="(b) Varying fraction of PO data: Variance") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))

pgplot + theme(legend.position="none")

ggsave(filename = "PO-percentage-var.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)

## [0,1.75] VARIANCE ----------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
{
  load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
  i <- 1
  j <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc1.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc2_5.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
  i <- 1
  j <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc7.RData')
  i <- 1
  j <- 8
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1.5)))
  
}
# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "PO percentage", title="(a) Variation of the estimator\n in the first half of domain") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))

pgplot + theme(legend.position="none")

ggsave(filename = "0175.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)


## [1.75,3] VARIANCE: VARIATION OF THE PERCENTAGE OF PARTIALLY OBSERVES CURVES ----------
levs <- 8 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs*B)
method      <- rep(c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt"),
                   each=B)
beta_var   <- numeric(3*levs*B)

# Riempio il beta_var
{
  load('Simulation/Results/repeated-simulations/KLAl-perc1_nowgts.RData')
  i <- 1
  j <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc1.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc2_5_nowgts.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc2_5.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations/KLAl-perc7_nowgts.RData')
  i <- 1
  j <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-perc7.RData')
  i <- 1
  j <- 8
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta0.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(1.5,3.5)))
  
}

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("10%", "10%: wgt", "25%", "25%: wgt", "40%", "40%: wgt", "70%", "70%: wgt")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
pgplot <- ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  scale_y_continuous(limits=c(0,0.1)) +
  theme_bw() +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b=1}^B \hat{\beta^b} ||_2^2$)'),
       fill = "PO percentage", title="(b) Variation of the estimator \n in the second half of domain") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size=22))
pgplot + theme(legend.position="none")

ggsave(filename = "17535.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 8,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)


## BIAS2 ----------------------------------------------------------------------------

## Per ogni metodo di ricostruzione e per ogni coefficiente, mi serve l'informazione sul bias

levs <- 6 #levels of factor method
coefficient <- rep(c("Beta0", "Beta1", "Beta2"), each=levs)
method      <- rep(c("Kraus:wgt", "KL-PC:wgt", "KL-AL:wgt", "Kraus", "KL-PC", "KL-AL"),
                   each=1)
bias2   <- numeric(3*levs)

# Riempio il beta_mse
{
  load('Simulation/Results/repeated-simulations-oldwgts/Kraus.RData')
  i <- 1 #index of the coefficient
  j <- 1 #index of the method
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLNoAl.RData')
  i <- 1
  j <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations-oldwgts/KLAl-par10.RData')
  i <- 1
  j <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations/Kraus_nowgts.RData')
  i <- 1
  j <- 4
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations/KLNoAl_nowgts.RData')
  i <- 1
  j <- 5
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  load('Simulation/Results/repeated-simulations/KLAl_nowgts.RData')
  i <- 1
  j <- 6
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta0.est)
  mean.diff <- mean0 - beta0
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean.diff <- mean0 - beta1
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1) + 1
  idx.sup <- (levs*(i-1)+j)
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean.diff <- mean0 - beta2
  bias2[idxs] <- inprod(mean.diff, mean.diff, 0, 0, rng=range(t.points))
  
}
bias2


## PLOT LEGENDS SEPARATELY ----------------------------------------------------------
## Weights
library(scales)
gc.ramp <- hue_pal()(7)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:10, ylim=0:10)
legend("bottom", legend =c("No weights","a=5","a=10","a=15", "a=20", "a=inf", "0-wgts"),
       cex=1.5, bty='n',
       fill = gc.ramp, horiz=T)
mtext("Parameter:", at=0.2, cex=2)

## Reconstruction methods 


## PO Fraction

