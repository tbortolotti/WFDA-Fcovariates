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
library(matrixcalc)
library(beepr)
library(latex2exp)
library(scales)

rm(list=ls())
graphics.off()
cat("\014")

## Load functions --------------------------------------------------------------
source('Simulation/methods/workflow_weighted_analysis.R')
source('Simulation/methods/generate_data.R')
load('Simulation/DATA/reg_info.RData')

## Simulate Data ---------------------------------------------------------------

# Set the case information
case.info <- list(n.sim       = 100,
                  noise       = 0.1,
                  ext.noise   = 0.5,
                  perc        = 0.4,
                  left.bound  = 1.5,
                  right.bound = 3.5)

seed <- 140996
simulated_data <- generate_data(seed      = seed,
                                case      = "CASE-1",
                                reg.info  = reg.info,
                                case.info = case.info)

t.points       <- simulated_data$t.points
T_hp           <- simulated_data$T_hp
curves         <- simulated_data$curves
curves.true.fd <- simulated_data$curves.true.fd
xlist          <- simulated_data$xlist

## Utilities -------------------------------------------------------------------
B <- 100
method <- 'KLAl'

#vec.par <- c(1,3,4,5,10,100)
#vec.par <- c(11,12,15)
#vec.par <- c(5,10,15,20)
vec.par <- c("unwgt")
MSE_cv <- matrix(data=0, nrow=B, ncol=length(vec.par))

beta0.list <- list()
beta1.list <- list()
beta2.list <- list()

(Start.Time <- Sys.time())
for(a in 1:length(vec.par)) #a=1
{
  fix.par <- vec.par[a]
  print(paste0('Parameter ', fix.par))
  pb <- progress_bar$new(total=B)
  for(b in 1:B) # b=1
  {
    method_evaluation <- workflow_weighted_analysis(b              = b,
                                                    B              = B,
                                                    t.points       = t.points,
                                                    breaks         = t.points,
                                                    T_hp           = T_hp,
                                                    curves         = curves,
                                                    curves.true.fd = curves.true.fd,
                                                    xlist          = xlist,
                                                    method         = method,
                                                    fix.par        = fix.par,
                                                    wgts.flag      = FALSE)
    
    MSE_cv[b,a] <- method_evaluation$MSE
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
  
  beta0.list[[a]] <- beta0.est
  beta1.list[[a]] <- beta1.est
  beta2.list[[a]] <- beta2.est
}

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

beep()

name.file <- paste0('Simulation/Results/LOO_',method,'_unwgt.RData')
save(vec.par, MSE_cv, beta0.list, beta1.list, beta2.list, file=name.file)

## RESULTS ----------------------------------------------------------------------
B <- 100
load('Simulation/Results/LOO_KLAl_unwgt.RData')
MSE_cv.mat <- matrix(data=0, nrow=B, ncol=1)
MSE_cv.mat[,1] <- MSE_cv

load('Simulation/Results/LOO_KLAl_oldwgts.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv)

load('Simulation/Results/LOO_KLAl_oldwgts_100.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv[,1])

#load('Simulation/Results/LOO_KLAl_0wgts.RData')
load('Simulation/Results/LOO_KLAl.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv[,6])

lamorte <- apply(MSE_cv.mat, MARGIN=2, FUN=mean)
which(lamorte==min(lamorte))

# BOXPLOTS ----------------------------------------------------------------
B <- 100
load('Simulation/Results/LOO_KLAl_unwgt.RData')
MSE_cv.mat <- matrix(data=0, nrow=B, ncol=1)
MSE_cv.mat[,1] <- MSE_cv

load('Simulation/Results/LOO_KLAl_oldwgts.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv)

load('Simulation/Results/LOO_KLAl_oldwgts_100.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv[,1])

#load('Simulation/Results/LOO_KLAl_0wgts.RData')
load('Simulation/Results/LOO_KLAl.RData')
MSE_cv.mat <- cbind(MSE_cv.mat, MSE_cv[,6])

MSE <- vec(MSE_cv.mat)

class <- c(rep("unwgt", B),
           rep('a=5', B),
           rep('a=10', B),
           rep('a=15', B),
           rep('a=20', B),
           rep('a=inf', B),
           rep('0-wgts', B))

data.box <- data.frame(MSE, as.factor(class))
names(data.box) <- c('MSE', 'class')
class_order<- c("unwgt", "a=5","a=10", "a=15", "a=20", "a=inf", "0-wgts")
data.box <- data.box %>% mutate(class = factor(x = class, levels = class_order))

box.dir <- '/Users/teresabortolotti/Desktop/flash/a-definition'
pg_plot <- ggplot(data.box, aes(x = class, y = MSE))+
  geom_boxplot(color="black", fill=hue_pal()(7))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0,0.04)) +
  theme(legend.position="None") + 
  theme_bw() +
  labs(x="", y=TeX(r'($ || \hat{y}_{i} - y_i ||_2^2$)'), title="Weights definition: LOO MSE") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

ggsave(filename = "LOO-MSE-simulation.pdf",
       plot = last_plot(),
       device = NULL,
       path = box.dir,
       scale = 1,
       width = 9,
       height = 5,
       units = "in",
       dpi = 300,
       limitsize = TRUE,
       bg = NULL)


## CASO VECCHIO ----------------------------------------------------------------
load('Simulation/Results/LOO_KLAl_111215.RData')
temp <- MSE_cv

load('Simulation/Results/LOO_KLAl.RData')
MSE_cv <- cbind(MSE_cv[,c(4,5)], temp[,c(1,2,3)])

lamorte <- apply(MSE_cv, MARGIN=2, FUN=mean)
which(lamorte==min(lamorte))
