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

## Load Data ---------------------------------------------------------------
load('DATA/curves.RData')
load('DATA/t_period.RData')
load('DATA/obs.RData')
load('DATA/T_hp.RData')
load('DATA/xlist.RData')
load('DATA/data.RData')
#load('DATA/loc_1/lambda_fixed.RData')
load('DATA/events.RData')

t.points <- log10(T.period)
t.points[1] <- -3

data <- list(dJB  = dJB,
             MAG  = MAG,
             SoF  = SoF,
             VS30 = VS30)

data.f <- data.frame(dJB = dJB, MAG = MAG, SoF = SoF, VS30 = VS30)

## Load Functions ----------------------------------------------------------
source('methods/find_obs_inc.R')
source('methods/interpolation.R')
source('methods/create_weights.R')
source('methods/wt_bsplinesmoothing.R')

source('methods/Regression/weighted_fRegress.R')
source('methods/Regression/my_predict_fRegress.R')
#source('Regression/my_predict_fRegress_sing.R')
source('methods/Regression/eval_MSE.R')

## Utilities -------------------------------------------------------------------
n <- dim(curves)[2]
q <- length(xlist)
reconst_fcts  <- find_obs_inc(Y = curves)

## Interpolation ---------------------------------------------------------------
interpolate   <- interpolation(curves       = curves,
                               t.points     = T.period,
                               T_hp         = T_hp,
                               reconst_fcts = reconst_fcts)
curves.interp <- interpolate$curves.rec

#x11(width=8000, height=5000)
#matplot(T.period, curves.interp[,reconst_fcts], type='l')

## Construction of the weights
loc <- 0.3
t.points <- log10(T.period)
t.points[1] <- -3

## Construction of the weights -------------------------------------------------
wgt              <- create_weights(curves        = curves.interp,
                                   t.points      = t.points,
                                   loc.par       = loc,
                                   reconst_fcts  = reconst_fcts,
                                   T_hp          = log10(T_hp))

## Smoothing -------------------------------------------------------------------
smth             <- wt_bsplinesmoothing(curves   = curves.interp,
                                        wgts.obs = wgt$wgts.obs,
                                        t.points = t.points,
                                        lambda   = 1e-5,
                                        set.cb   = FALSE)
curves.interp.fd <- smth$curves.fd

#y2cmaps <- smth$y2cmaps
#save(y2cmaps, file='DATA/y2cmaps.RData')

## Sphi for confidence bands ---------------------------------------------------
# source('methods/build_Sphi.R')
# build_Sphi(y2cmaps = smth$y2cmaps, curves = curves.interp.fd)

## B-list ----------------------------------------------------------------------
load('blist_options/blist_redo.RData')

## Regression and beta estimation ----------------------------------------------
mod <- weighted_fRegress(y            = curves.interp.fd,
                         xfdlist      = xlist,
                         betalist     = blist,
                         wgts         = wgt$wgts.fd)

## y_hat   ---------------------------------------------------------------------
curves.interp.hat   <- my_predict_fRegress(mod          = mod,
                                           xlist        = xlist,
                                           t.points     = t.points)
curves.interp.hat.v <- eval.fd(t.points, curves.interp.hat)

## MSE evaluation - through an event-wise crossvalidation ----------------------
source('methods/Regression/pwMSE.R')
pwMSE.val <- pwMSE(curves    = curves,
                   curves.fd = curves.interp.fd,
                   xlist     = xlist,
                   t.points  = t.points,
                   events    = event.id,
                   blist     = blist,
                   B         = 10,
                   wgts.fd   = wgt$wgts.fd,
                   set.ITA18 = TRUE,
                   data      = data.f)

#### DIAGNOSTIC ----------------------------------------------------------------

## Coefficients standard errors ------------------------------------------------
library(Matrix)

source("methods/stderrors.R")
load('DATA/Sphi.RData')

res       <- curves.interp.hat - curves.interp.fd
E         <- t(eval.fd(T.period, res))
SigmaE    <- 1/(n-q)*t(E)%*%E

st.err <- stderrors(Sphi         = Sphi,
                    obs.inc      = reconst_fcts,
                    SigmaE       = SigmaE,
                    returnMatrix = FALSE,
                    fRegressList = mod)

save(st.err, file='DATA/st_err.RData')

## PLOTS -----------------------------------------------------------------------
name_dir     <- paste0("fitITA18")
input        <- list(xlist    = xlist,
                     t.points = t.points,
                     blist    = blist)

## Save plots
source("methods/plots/save_plots.R")
save_plots(mod.fit=mod, input=input, name_dir=name_dir)

## Goodness of fit
source('methods/plots/goodnessoffit.R')
goodnessoffit(mod.fit=mod, input=input, data=data, name_dir=name_dir)

## Model comparison ------------------------------------------------------------
# source('methods/fit_ITA18.R')
# mod.ITA18 <- fit_ITA18(data.f, curves)
# 
# save(mod.ITA18, file='DATA/mod_ITA18.RData')
load('DATA/mod_ITA18.RData')

coefs.ITA18 <- mod.ITA18$coefficients
sigma.ITA18 <- mod.ITA18$sigma.ITA18
y.hat.ITA18 <- mod.ITA18$y.hat.ITA18

## Plot sigma comparison
source('methods/plots/plot_sigma.R')

res       <- curves.interp.hat - curves.interp.fd
E         <- t(eval.fd(t.points, res))
SigmaE    <- 1/(n-q-1)*t(E)%*%E

plot_sigma(SigmaE, sigma.ITA18, t.points, name_dir)


## Plot MSE comparison
source('methods/plots/plotMSE_pw.R')
plotMSE_pw(MSE.vec   = pwMSE.val$MSE_t,
           MSE.ita18 = pwMSE.val$MSE_t18,
           t.points  = t.points,
           name_dir  = name_dir)

# source('plotMSE_pw_discuss.R')
# plotMSE_pw_discuss(MSE.vec   = pwMSE.val$MSE_t,
#                    MSE.ita18 = pwMSE.val$MSE_t18,
#                    t.points  = T.period,
#                    name_dir  = name_dir,
#                    set.log   = FALSE,
#                    plot.name = "MSE")

## Prediction comparison
source("methods/plots/model_comparison.R")
t.idxs <- c(2,7,21)
model_comparison(mod.fit    = mod,
                 t.points   = t.points,
                 name_dir   = name_dir,
                 t.idxs     = t.idxs,
                 data       = data.f,
                 curves     = curves.interp,
                 corrective = FALSE)

# source("model_comparison_discuss.R")
# t.idxs <- c(2,21)
# model_comparison_discuss(mod.fit    = mod,
#                          t.points   = T.period,
#                          name_dir   = name_dir,
#                          t.idxs     = t.idxs,
#                          data       = data,
#                          curves     = curves.interp,
#                          set.log    = FALSE)

## Prediction comparison with ITA18 corrected
source("methods/plots/model_comparison_corrected.R")
t.idxs <- c(2,7,21)
model_comparison_corrected(mod.fit    = mod,
                           t.points   = t.points,
                           name_dir   = name_dir,
                           t.idxs     = t.idxs,
                           data       = data,
                           curves     = curves.interp)

## Near source comparison with ITA18
source("methods/plots/plot_nearsource_comparison.R")
t.idxs <- c(2,7,21)
plot_nearsource_comparison(mod.fit    = mod,
                           t.points   = t.points,
                           name_dir   = name_dir,
                           t.idxs     = t.idxs,
                           data       = data,
                           curves     = curves.interp)

# source("plot_nearsource_comparison_discuss.R")
# t.idxs <- c(2,7)
# plot_nearsource_comparison_discuss(mod.fit    = mod,
#                                    t.points   = T.period,
#                                    name_dir   = name_dir,
#                                    t.idxs     = t.idxs,
#                                    data       = data,
#                                    curves     = curves.interp,
#                                    set.log    = FALSE)

## Magnitude check --> a check that the spectral acceleration does not diminish
#                      as magnitude increases
source("methods/plots/magnitude_check.R")
t.idxs <- c(2,7,21)
magnitudes <- c(4.0, 5.0, 6.0, 7.0, 8.0)
magnitude_check(mod.fit    = mod,
                t.points   = t.points,
                name_dir   = name_dir,
                t.idxs     = t.idxs,
                magnitudes = magnitudes, 
                data       = data.f,
                curves     = curves.interp)

## Confidence bands
load('WFA/DATA/st_err_loc1.RData')
source('plot_CB.R')
plot_CB(st.err = st.err, name_dir, alpha=0.05, conf = FALSE)

source('plot_regressors.R')
plot_regressors(st.err=st.err, name_dir, intercept.idxs = c(7,21,30,37))


## Source ----------------------------------------------------------------------
source('methods/plots/Source.R')
dir.current <- getwd()
my.dir <- paste0(dir.current,"/Results/",name_dir,"/Source")
dir.create(my.dir)

for(t.idx in 1:37)
{
  Source(my.dir=my.dir, name_dir=name_dir, data=data, t.idx=t.idx)
}

## Path ------------------------------------------------------------------------
source('methods/plots/Path.R')
dir.current <- getwd()
my.dir <- paste0(dir.current,"/Results/",name_dir,"/Path")
dir.create(my.dir)
for(t.idx in 1:37)
{
  Path(my.dir=my.dir, name_dir=name_dir, data=data, t.idx=t.idx)
}

source('methods/plots/Path_unified.R')
dir.current <- getwd()
my.dir <- paste0(dir.current,"/Results/",name_dir,"/Path_unified")
dir.create(my.dir)
for(t.idx in 1:37)
{
  Path_unified(my.dir=my.dir, name_dir=name_dir, data=data, t.idx=t.idx)
}

## Site ------------------------------------------------------------------------
source('methods/plots/Site.R')
dir.current <- getwd()
my.dir <- paste0(dir.current,"/Results/",name_dir,"/Site")
dir.create(my.dir)

for(t.idx in 1:37)
{
  Site(my.dir=my.dir, name_dir=name_dir, data=data, t.idx=t.idx)
}

