setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

library(fda)
library(devtools)
library(fastmatrix)
library(calculus)
library(ReconstPoFD)
library(snowfall)

rm(list=ls())
graphics.off()
cat("\014")

## Load Data ---------------------------------------------------------------
load('DATA/curves.RData')
load('DATA/t_period.RData')
load('DATA/obs.RData')
load('DATA/T_hp.RData')
load('DATA/xlist.RData')
load('DATA/events.RData')
load('blist_options/blist_redo.RData')

## Function for the evaluation of the MAE for each proposed method
source('methods/methods_workflow.R')

## Logarithm of the period --------------------------------------------------
t.points    <- log10(T.period)
t.points[1] <- -3

log.Thp     <- log10(T_hp)

## Cross-validation --------------------------------------------------

B     <- 5 #10

method <- 'interpolation'
#method <- 'interpolation-noweight'
#method <- 'Kraus1'
#method <- 'Kraus2'
#method <- 'KLNoAl' 
#method <- 'KLAl'

loc.vec <- seq(0, 1, by=0.2)
i       <- 1
loc     <- loc.vec[i]


### Evaluate simulation time
(Start.Time <- Sys.time())


MSE_cv <- numeric(B)
MSE_reconstruction_cv <- numeric(B)

pb <- progressBar(0, max = B, initial = 0, style = "ETA")
for(b in 1:B)
{
  method_evaluation <- method_comparison(b          = b,
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
  
  setTxtProgressBar(pb, b)
  
}


name.file <- paste0('noweight/MSE_',method,'_',loc,'.RData')

save(MSE_cv, MSE_reconstruction_cv, MSE, MSE_reconstruction, file=name.file)

##------------------------------------
End.Time <- Sys.time()
## Run-time:
round(End.Time - Start.Time, 2)
##------------------------------------



# Store results ----------------------------------------------------
MSEs <- numeric(9)
MSEs_recon <- numeric(9)

load('MSE_vals_lambdasel/MSE_interpolation_0.RData')
MSEs[1] <- MSE
MSEs_recon[1] <- MSE_reconstruction

load('MSE_vals_lambdasel/MSE_interpolation_1.RData')
MSEs[2] <- MSE
MSEs_recon[2] <- MSE_reconstruction

load('MSE_vals_lambdasel/MSE_interpolation_2.RData')
MSEs[3] <- MSE
MSEs_recon[3] <- MSE_reconstruction

load('MSE_vals_lambdasel/MSE_interpolation_3.RData')
MSEs[4] <- MSE
MSEs_recon[4] <- MSE_reconstruction

load('MSE_vals_lambdasel/MSE_interpolation_4.RData')
MSEs[5] <- MSE
MSEs_recon[5] <- MSE_reconstruction

load('MSE_vals_lambdasel/MSE_Kraus1.RData')
MSEs[6] <- MSE
MSEs_recon[6] <- MSE_reconstruction

load('MSE_vals_lambdasel/MSE_Kraus2.RData')
MSEs[7] <- MSE
MSEs_recon[7] <- MSE_reconstruction

load('MSE_vals_lambdasel/MSE_KLNoAl.RData')
MSEs[8] <- MSE
MSEs_recon[8] <- MSE_reconstruction

load('MSE_vals_lambdasel/MSE_KLAl.RData')
MSEs[9] <- MSE
MSEs_recon[9] <- MSE_reconstruction


## Analyse results ----------------------------------------
MSEs.rel <- MSEs/min(MSEs)
MSEs.rel

MSEs
MSEs_recon
