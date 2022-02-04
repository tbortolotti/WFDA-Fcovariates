setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

library(fda)
library(devtools)
library(fastmatrix)
library(calculus)
library(ReconstPoFD)
library(snowfall)
library(progress)
library(ggplot2)
library(fields)

rm(list=ls())
graphics.off()
cat("\014")

## Load Data -------------------------------------------------------------------
load('DATA/curves.RData')
load('DATA/t_period.RData')
load('DATA/obs.RData')
load('DATA/T_hp.RData')
load('DATA/xlist-log.RData')
load('DATA/events.RData')
load('blist_options/log/blist.RData')

## Function for the evaluation of the MSE for each proposed method
source('methods/Reconstruction/methods_workflow_log.R')

# Logarithm of the period
t.points    <- log10(T.period)
t.points[1] <- -2.4

log.Thp     <- log10(T_hp)

## Cross-validation ------------------------------------------------------------

B     <- 5 #10

method <- 'extrapolation'
#method <- 'extrapolation-noweight'
#method <- 'Kraus1'
#method <- 'Kraus2'
#method <- 'KLNoAl' 
#method <- 'KLAl'


loc.vec <- c(0,1,2,3,4,5)
i       <- 6
loc     <- loc.vec[i]


MSE_cv <- numeric(B)
MSE_reconstruction_cv <- numeric(B)

(Start.Time <- Sys.time())
pb <- progress_bar$new(total=B)
for(b in 1:B)
{
  method_evaluation <- methods_workflow_log(b          = b,
                                            T.period   = T.period,
                                            t.points   = t.points,
                                            T_hp       = T_hp,
                                            log.Thp    = log.Thp,
                                            curves     = curves,
                                            events     = event.id,
                                            B          = B,
                                            xlist      = xlist,
                                            blist      = blist,
                                            method     = method,
                                            loc        = loc)
  
  MSE_cv[b] <- method_evaluation$MSE
  MSE_reconstruction_cv[b] <- method_evaluation$MSE_reconstruction
  
  MSE <- mean(MSE_cv)
  MSE_reconstruction <- mean(MSE_reconstruction_cv)
  
  pb$tick()
  
}
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

name.file <- paste0('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_',method,'_',(i-1),'.RData')
save(MSE_cv, MSE_reconstruction_cv, MSE, MSE_reconstruction, file=name.file)
mean(MSE_cv)

#name.file <- paste0('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_',method,'.RData')
#save(MSE_cv, MSE_reconstruction_cv, MSE, MSE_reconstruction, file=name.file)

load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_extrapolation-noweight.RData')
mean(MSE_cv)

# Store results ----------------------------------------------------------------
MSEs <- numeric(8)
MSEs_recon <- numeric(8)

load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_extrapolation_0.RData')
MSEs[1] <- MSE
MSEs_recon[1] <- MSE_reconstruction
extrap0 <- MSE_cv

# load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_extrapolation_1.RData')
# MSEs[2] <- MSE
# MSEs_recon[2] <- MSE_reconstruction
# extrap1 <- MSE_cv
# 
load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_extrapolation_2.RData')
MSEs[2] <- MSE
MSEs_recon[2] <- MSE_reconstruction
extrap2 <- MSE_cv

# load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_extrapolation_3.RData')
# MSEs[2] <- MSE
# MSEs_recon[2] <- MSE_reconstruction
# extrap3 <- MSE_cv

# load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_extrapolation_4.RData')
# MSEs[5] <- MSE
# MSEs_recon[5] <- MSE_reconstruction
# extrap4 <- MSE_cv

load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_extrapolation_5.RData')
MSEs[3] <- MSE
MSEs_recon[3] <- MSE_reconstruction
extrap5 <- MSE_cv

load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_extrapolation-noweight.RData')
MSEs[4] <- MSE
MSEs_recon[4] <- MSE_reconstruction
extrap_no <- MSE_cv

load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_Kraus1.RData')
MSEs[5] <- MSE
MSEs_recon[5] <- MSE_reconstruction
Kraus1 <- MSE_cv

load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_Kraus2.RData')
MSEs[6] <- MSE
MSEs_recon[6] <- MSE_reconstruction
Kraus2 <- MSE_cv

load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_KLNoAl.RData')
MSEs[7] <- MSE
MSEs_recon[7] <- MSE_reconstruction
KLNoAl <- MSE_cv

load('Results/Reconstruction-methods-comparison/5-fold/novel/MSE_KLAl.RData')
MSEs[8] <- MSE
MSEs_recon[8] <- MSE_reconstruction
KLAl <- MSE_cv

MSEs.rel <- MSEs/min(MSEs)
MSEs.rel

MSEs
MSEs_recon

## Boxplots of the comparison ------------------------------------------------------------------

# class <- c(rep('extrap 0',B),
#            rep('extrap 1', B),
#            rep('extrap 2', B),
#            rep('extrap 3', B),
#            rep('extrap 4', B),
#            rep('extrap no-wgts', B),
#            rep('Kraus 1', B),
#            rep('Kraus 2', B),
#            rep('KL - No Al', B),
#            rep('KL - Al', B))
# 
# data.box <- data.frame(c(extrap0,
#                          extrap1,
#                          extrap2,
#                          extrap3,
#                          extrap4,
#                          extrap_no,
#                          Kraus1,
#                          Kraus2,
#                          KLNoAl,
#                          KLAl), as.factor(class))

class <- c(rep('extrap 0',B),
           rep('extrap 2', B),
           rep('extrap 5', B),
           rep('extrap no-wgts', B),
           rep('Kraus 1', B),
           rep('Kraus 2', B),
           rep('KL - No Al', B),
           rep('KL - Al', B))

data.box <- data.frame(c(extrap0,
                         extrap2,
                         extrap5,
                         extrap_no,
                         Kraus1,
                         Kraus2,
                         KLNoAl,
                         KLAl), as.factor(class))
names(data.box) <- c('MSE', 'class')

pg_plot <- ggplot(data.box, aes(x = class, y = MSE))+
  geom_boxplot(color="black", fill=tim.colors(8))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  scale_y_continuous(limits=c(0.08,0.2)) +
  labs(y="5-fold Mean Squared Error", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")



