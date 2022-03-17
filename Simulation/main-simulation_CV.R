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

## Load Data -------------------------------------------------------------------
load('Simulation/DATA/CASE1.RData')
#matplot(t.points, curves, type='l')

## Load Functions --------------------------------------------------------------
source('Simulation/methods/workflow_weighted_analysis.R')

## Utilities -------------------------------------------------------------------
n <- dim(curves)[2]
fix.par <- 2

# List of covariates
x1 <- as.matrix(regression.list$x1)
intercept <- as.matrix(rep(1,n))
xlist <- list(intercept, x1, regression.list$x2.fd)

B <- 10
MSE_cv <- numeric(B)

#method <- 'Kraus'
method <- 'KLNoAl'
method <- 'KLAl'
method <- 'KLAl-new'

(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B) # b=1
{
  method_evaluation <- workflow_weighted_analysis(b              = b,
                                                  B              = B,
                                                  t.points       = t.points,
                                                  breaks         = t.points,
                                                  T_hp           = T_hp,
                                                  curves         = curves,
                                                  curves.true    = curves.true,
                                                  curves.true.fd = curves.true.fd,
                                                  xlist          = xlist,
                                                  method         = method,
                                                  fix.par        = fix.par,
                                                  wgts.flag      = TRUE)
  
  MSE_cv[b] <- method_evaluation$MSE
  
  MSE <- mean(MSE_cv)
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

name.file <- paste0('Simulation/Results/functional-MSE/MSE_',method,'.RData')
save(MSE_cv, MSE, file=name.file)

# Store results ----------------------------------------------------------------
MSEs <- numeric(6)

load('Simulation/Results/scalar-MSE/MSE_Kraus.RData')
MSEs[1] <- MSE
Kraus <- MSE_cv

load('Simulation/Results/scalar-MSE/MSE_KLNoAl.RData')
MSEs[2] <- MSE
KLNoAl <- MSE_cv

# load('Simulation/Results/functional-MSE/MSE_KLAl-new.RData')
# MSEs[2] <- MSE
# MSEs_recon[2] <- MSE_reconstruction
# KLNoAl <- MSE_cv

load('Simulation/Results/scalar-MSE/MSE_KLAl.RData')
MSEs[3] <- MSE
KLAl <- MSE_cv

load('Simulation/Results/scalar-MSE/MSE_Kraus_nowgts.RData')
MSEs[4] <- MSE
Kraus_nowgt <- MSE_cv

load('Simulation/Results/scalar-MSE/MSE_KLNoAl_nowgts.RData')
MSEs[5] <- MSE
KLNoAl_nowgt <- MSE_cv

load('Simulation/Results/scalar-MSE/MSE_KLAl_nowgts.RData')
MSEs[6] <- MSE
KLAl_nowgt <- MSE_cv

MSEs.rel <- MSEs/min(MSEs)
MSEs.rel

MSEs
MSEs_recon

## Boxplots of the comparison ------------------------------------------------------------------

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
  #scale_y_continuous(limits=c(0.08,0.2)) +
  labs(y="10-fold Mean Squared Error", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")


