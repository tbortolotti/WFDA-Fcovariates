setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

rm(list=ls())
graphics.off()
cat("\014")

library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(fastmatrix)

## DATA EXTRACTION -------------------------------------------------------------

# Joyner-Boore distance
data <- readMat('DATA/dist.mat')
dJB <- data[[1]]

# Moment Magnitude
data <- readMat('DATA/MAG.mat')
MAG <- data[[1]]

# Style of Faulting
data <- readMat('DATA/sof.mat')
SoF <- c()
for(i in 1:length(data[[1]]))
{
  SoF[i] <- as.character(data[[1]][[i]])
}

# Share-wave velocity V_S30
data <- readMat('DATA/vs30.mat')
VS30 <- data[[1]]

# Curves
data <- readMat('DATA/DATA.mat')
U_hp <- c()
V_hp <- c()
for(i in 1:length(data[[1]][[1]]))
{
  U_hp[i] <- as.numeric(data[[1]][[1]][[i]])
  V_hp[i] <- as.numeric(data[[1]][[2]][[i]])
}


rotD50 <- list()
for(i in 3:length(data[[1]]))
{
  rotD50[[i-2]] <- data[[1]][[i]]
}

## Save curves values in matrix
rotD50.m <- matrix(data=0, nrow=5778, ncol=37)
for(i in 1:length(rotD50))
{
  for(j in 1:length(rotD50[[1]]))
  {
    rotD50.m[j,i] <- as.numeric(rotD50[[i]][[j]])
  }
}

T.period <- c(0,0.010,0.025,0.040,0.050,0.070,0.100,
          0.150,0.200,0.250,0.300,0.350,0.400,0.450,
          0.500,0.600,0.700,0.750,0.800,0.900,1.000,
          1.200,1.400,1.600,1.800,2.000,2.500,3.000,
          3.500,4.000,4.500,5.000,6.000,7.000,8.000,
          9.000,10.000)

#filtered
rotD50.f <- rotD50.m
# period filtering
for(i in 1:dim(rotD50.m)[1])
{
  f_max <- max(U_hp[i],V_hp[i])
  T_max <- 1/f_max
  for(j in 1:length(T.period))
  {
    if(T.period[j]>T_max)
    {
      rotD50.f[i,j] <- NA
    }
  }
}


## EVENT
# load event characteristics
data <- readMat('DATA/event.mat')
event.id <- list()
for(i in 1:length(data[[1]][[1]]))
{
  event.id[[i]] <- as.character(data[[1]][[1]][[i]][[1]])
}

uni <- unique(event.id)
aux <- numeric(length(event.id))
for(i in 1:length(uni))
{
  aux[which(event.id==uni[[i]])] <- i
}
event.id <- aux

event.lat <- c()
for(i in 1:length(data[[1]][[1]]))
{
  event.lat[i] <- as.numeric(data[[1]][[2]][[i]][[1]])
}

event.long <- c()
for(i in 1:length(data[[1]][[1]]))
{
  event.long[i] <- as.numeric(data[[1]][[3]][[i]][[1]])
}


## Remove missing values
id.na <- which(is.na(dJB) | is.na(MAG) | is.na(SoF) | is.na(VS30) | is.na(U_hp) |
                 is.na(V_hp) | is.na(event.id))

dJB <- dJB[-id.na]
MAG <- MAG[-id.na]
SoF <- SoF[-id.na]
VS30 <- VS30[-id.na]
U_hp <- U_hp[-id.na]
V_hp <- V_hp[-id.na]
event.id <- event.id[-id.na]
event.lat <- event.lat[-id.na]
event.long <- event.long[-id.na]
rotD50.f <- rotD50.f[-id.na,]
rotD50.m <- rotD50.m[-id.na,]

## Incomplete records analysis -------------------------------------------------
# count NA's per observation of unfiltered matrix
na.m <- c()
for(i in 1:dim(rotD50.m)[1])
{
  na.m[i] <- sum(is.na(rotD50.m[i,]))
}
plot(na.m, lwd=1, col='darkblue',type='h', main="NAs before filtering")

# Remove curves that are completely unobserved
dJB <- dJB[-which(na.m==37)]
MAG <- MAG[-which(na.m==37)]
VS30 <- VS30[-which(na.m==37)]
SoF <- SoF[-which(na.m==37)]
U_hp <- U_hp[-which(na.m==37)]
V_hp <- V_hp[-which(na.m==37)]
event.id <- event.id[-which(na.m==37)]
event.lat <- event.lat[-which(na.m==37)]
event.long <- event.long[-which(na.m==37)]
rotD50.m <- rotD50.m[-which(na.m==37),]
rotD50.f <- rotD50.f[-which(na.m==37),]

x11()
par(mfrow=c(2,1))
matplot(T.period, t(rotD50.m), type='l', main="Unfiltered observations")
matplot(T.period, t(rotD50.f), type='l', main="Filtered observations")

# Incomplete records analysis
F_hp <- c()
T_hp <- c()
for(i in 1:length(U_hp))
{
  F_hp[i] <- max(V_hp[i],U_hp[i])
  T_hp[i] <- min(1/F_hp[i],10)
}

length(which(T_hp < 10))
x11()
par(mfrow=c(1,2))
hist(T_hp,ylim=c(0,length(T_hp)), xlab="T_max", ylab="", main="Occurrences of T_max in the observations")

# Percentage of observations that have T_hp<T.period, as function of T.period
prop <- c()
for(t in 1:length(T.period))
{
  prop[t] <- length(which(T_hp>=T.period[t]))/length(T_hp) 
}
plot(T.period, prop, ylim=c(0.5,1), xlab="", ylab="", type='o-', main="% complete records per period")

# which observations have incomplete domain
na.f <- c()
for(i in 1:dim(rotD50.f)[1])
{
  na.f[i] <- sum(is.na(rotD50.f[i,]))
}
obs.inc <- which(na.f>0)

x11()
par(mfrow=c(2,1))
matplot(T.period, t(rotD50.f), ylim=c(0,5000), type='l', main="Complete records")
matplot(T.period, t(rotD50.f[obs.inc,]), ylim=c(0,5000), type='l', main="Incomplete records")


curves <- log10(t(rotD50.f))
save(curves, file='DATA/curves.RData')
save(T.period, file='DATA/t_period.RData')
save(T_hp, file='DATA/T_hp.RData')
save(event.id, event.lat, event.long, file='DATA/events.RData')

n <- dim(curves)[2]
obs.comp <- (1:n)[-obs.inc]

save(obs.comp, obs.inc, file='DATA/obs.RData')
save(rotD50.f, dJB, MAG, VS30, SoF, U_hp, V_hp, T.period, file='DATA/data.RData')

## Regressors construction ------------------------------------------------------
load('DATA/ITA18_regressors.RData')
load('DATA/t_period.RData')

#scaletta
#1. rendi Mh, Mref e h degli fd objects. Mi serve per capire quale sia l'fPar
#   ottimo per smoothare i regressori
#2. Costruisci gli fd objects dei regressori.

## 1. Build Mh, Mref and h functional objects
names(ITA18.regressors)

par(mfrow=c(3,1))
plot(T.period, ITA18.regressors$Mh.vec, type='l', xlab='Period', ylab='Mh')
plot(T.period, ITA18.regressors$Mref.vec, type='l', xlab='Period', ylab='Mref')
plot(T.period, ITA18.regressors$h.vec, type='l', xlab='Period', ylab='h')

t.points <- log10(T.period)
t.points[1] <- -2.4

par(mfrow=c(3,1))
plot(t.points, ITA18.regressors$Mh.vec, type='l', xlab='log10(T)', ylab='Mh')
plot(t.points, ITA18.regressors$Mref.vec, type='l', xlab='log10(T)', ylab='Mref')
plot(t.points, ITA18.regressors$h.vec, type='l', xlab='log10(T)', ylab='h')

## Functional covariates on T ## -------------------------------------------
t.points <- T.period
breaks <- c(seq(0,1,0.1), seq(2,10,0.5))

# Mh
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=1)
Mh.fd <- smooth.basis(t.points, ITA18.regressors$Mh.vec, fPar)
plot(Mh.fd)
lines(t.points, ITA18.regressors$Mh.vec, type='l', xlab='Period', ylab='Mh', col='red')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=1)
Mh.fd <- smooth.basis(t.points, ITA18.regressors$Mh.vec, fPar)
plot(Mh.fd)
lines(t.points, ITA18.regressors$Mh.vec, type='l', xlab='Period', ylab='Mh', col='red')

# Mref
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,10, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, ITA18.regressors$Mref.vec, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt #0.01

plot(log10(lambda.vec), gcv.vec, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.opt)
Mref.fd <- smooth.basis(t.points, ITA18.regressors$Mref.vec, fPar)
plot(Mref.fd)
lines(t.points, ITA18.regressors$Mref.vec, type='l', xlab='Period', ylab='Mref', col='red')

# h
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, ITA18.regressors$h.vec, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt #0.01

plot(log10(lambda.vec)[2:5], gcv.vec[2:5], type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
h.fd <- smooth.basis(t.points, ITA18.regressors$h.vec, fPar)
plot(h.fd)
lines(t.points, ITA18.regressors$h.vec, type='l', xlab='Period', ylab='h', col='red')


## 2. Create the list of functional regressors
load('DATA/data.RData')

n <- dim(rotD50.f)[1]
N <- length(t.points)

# Create the matrix of values assumed by the regressors at the sampling instants
reg.Mlow <- matrix(data=0, nrow=N, ncol=n)
reg.Mhigh <- matrix(data=0, nrow=N, ncol=n)
reg.D1 <- matrix(data=0, nrow=N, ncol=n)
reg.D2 <- matrix(data=0, nrow=N, ncol=n)
reg.D3 <- matrix(data=0, nrow=N, ncol=n)

R.mat <- matrix(data=0, nrow=N, ncol=n)

for(t in 1:N)
{
  Mh   <- ITA18.regressors$Mh.vec[t]
  Mref <- ITA18.regressors$Mref.vec[t]
  h    <- ITA18.regressors$h.vec[t]
  R    <- sqrt(dJB^2 + h^2)
  R.mat[t,] <- R
  
  reg.Mlow[t,]  <- ifelse(MAG<=Mh, MAG - Mh, 0)
  reg.Mhigh[t,] <- ifelse(MAG>Mh, MAG - Mh, 0)
  reg.D1[t,]    <- (MAG - Mref)*log10(R)
  reg.D2[t,]    <- log10(R)
  reg.D3[t,]    <- R
}

reg.SS      <- ifelse(SoF=="SS", 1, 0)
reg.TF      <- ifelse(SoF=="TF", 1, 0)
reg.S       <- ifelse(VS30<=1500, log10(VS30/800), log10(1500/800))

# Create the fd objects by projecting the regressors over an appropriate functional basis

# reg.Mlow
matplot(t.points, reg.Mlow, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, R.mat, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt #0.01

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.opt)
Mlow.fd <- smooth.basis(t.points, reg.Mlow, fPar)$fd
plot(Mlow.fd)

# I decide to force it to be smoother
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=1) 
Mlow.fd <- smooth.basis(t.points, reg.Mlow, fPar)$fd
plot(Mlow.fd)

# reg.Mhigh
matplot(t.points, reg.Mhigh, type='l')

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=1)
Mhigh.fd <- smooth.basis(t.points, reg.Mhigh, fPar)$fd
plot(Mhigh.fd)

# R
matplot(t.points, R.mat, type='l', main='True R')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, R.mat, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.opt)
R.fd <- smooth.basis(t.points, R.mat, fPar)$fd
plot(R.fd, main='Smoothed R')

# Look what happens for short and long distances: indexes 21 - 36
r <- R.mat[,21]
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, r, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.opt)
r.fd <- smooth.basis(t.points, r, fPar)
plot(r.fd)
lines(t.points, r, type='l', xlab='Period', ylab='R - idx 21', col='red')

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
r.fd <- smooth.basis(t.points, r, fPar)
plot(r.fd)
lines(t.points, r, type='l', xlab='Period', ylab='R - idx 21', col='red')

r <- R.mat[,36]
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, r, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.opt)
r.fd <- smooth.basis(t.points, r, fPar)
plot(r.fd)
lines(t.points, r, type='l', xlab='Period', ylab='R - idx 21', col='red')

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
r.fd <- smooth.basis(t.points, r, fPar)
plot(r.fd)
lines(t.points, r, type='l', xlab='Period', ylab='R - idx 21', col='red')

# reg.D1
matplot(t.points, reg.D1, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D1, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt #0.01

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.01)
D1.fd <- smooth.basis(t.points, reg.D1, fPar)$fd
plot(D1.fd)

# reg.D2
matplot(t.points, reg.D2, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D2, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt #0.01

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.01) #
D2.fd <- smooth.basis(t.points, reg.D2, fPar)$fd
plot(D2.fd)

# reg.D3
matplot(t.points, reg.D3, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D3, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt #0.01

d3 <- reg.D3[,21]
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
d3.fd <- smooth.basis(t.points, d3, fPar)$fd
plot(d3.fd, t.points, type='l')
lines(t.points, d3, type='l', xlab='Period', ylab='Reg D3 - idx 21', col='red')

d3 <- reg.D3[,21]
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.01)
d3.fd <- smooth.basis(t.points, d3, fPar)$fd
plot(d3.fd, t.points, type='l')
lines(t.points, d3, type='l', xlab='Period', ylab='Reg D3 - idx 36', col='red')


fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.01)
D3.fd <- smooth.basis(t.points, reg.D3, fPar)$fd
plot(D3.fd)

# Project constant coefficients over a constant basis
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)

reg.SS <- t(replicate(N, reg.SS))
SS.fd <- smooth.basis(t.points, reg.SS, fPar)$fd
par(mfrow=c(2,1))
matplot(t.points, reg.SS, type='l')
plot(SS.fd)

reg.TF <- t(replicate(N, reg.TF))
TF.fd <- smooth.basis(t.points, reg.TF, fPar)$fd

reg.S <- t(replicate(N, reg.S))
S.fd <- smooth.basis(t.points, reg.S, fPar)$fd


# Create the intercept
intercept <- t(replicate(N, rep(1,n)))
intercept.fd <- smooth.basis(t.points, intercept, fPar)$fd

## Finally build the list of functional regressors
xlist <- list(intercept.fd, Mlow.fd, Mhigh.fd, SS.fd, TF.fd, D1.fd, D2.fd, D3.fd, S.fd)
save(xlist, file='DATA/xlist.RData')

## Functional covariates on log(T) ## --------------------------------------------------
t.points <- log10(T.period)
t.points[1] <- -2.4
breaks <- t.points
# Mh
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=1)
Mh.fd <- smooth.basis(t.points, ITA18.regressors$Mh.vec, basis)
plot(Mh.fd)
lines(t.points, ITA18.regressors$Mh.vec, type='l', xlab='Period', ylab='Mh', col='red')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
Mh.fd <- smooth.basis(t.points, ITA18.regressors$Mh.vec, fPar)
plot(Mh.fd)
lines(t.points, ITA18.regressors$Mh.vec, type='l', xlab='Period', ylab='Mh', col='red')


# Mref
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,10, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, ITA18.regressors$Mref.vec, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

plot(log10(lambda.vec), gcv.vec, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
Mref.fd <- smooth.basis(t.points, ITA18.regressors$Mref.vec, fPar)
plot(Mref.fd)
lines(t.points, ITA18.regressors$Mref.vec, type='l', xlab='Period', ylab='Mref', col='red')

# h
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, ITA18.regressors$h.vec, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

plot(log10(lambda.vec)[2:5], gcv.vec[2:5], type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
h.fd <- smooth.basis(t.points, ITA18.regressors$h.vec, fPar)
plot(h.fd)
lines(t.points, ITA18.regressors$h.vec, type='l', xlab='Period', ylab='h', col='red')


## 2. Create the list of functional regressors
load('DATA/data.RData')

n <- dim(rotD50.f)[1]
N <- length(t.points)

# Create the matrix of values assumed by the regressors at the sampling instants
reg.Mlow <- matrix(data=0, nrow=N, ncol=n)
reg.Mhigh <- matrix(data=0, nrow=N, ncol=n)
reg.D1 <- matrix(data=0, nrow=N, ncol=n)
reg.D2 <- matrix(data=0, nrow=N, ncol=n)
reg.D3 <- matrix(data=0, nrow=N, ncol=n)

R.mat <- matrix(data=0, nrow=N, ncol=n)

for(t in 1:N)
{
  Mh   <- ITA18.regressors$Mh.vec[t]
  Mref <- ITA18.regressors$Mref.vec[t]
  h    <- ITA18.regressors$h.vec[t]
  R    <- sqrt(dJB^2 + h^2)
  R.mat[t,] <- R
  
  reg.Mlow[t,]  <- ifelse(MAG<=Mh, MAG - Mh, 0)
  reg.Mhigh[t,] <- ifelse(MAG>Mh, MAG - Mh, 0)
  reg.D1[t,]    <- (MAG - Mref)*log10(R)
  reg.D2[t,]    <- log10(R)
  reg.D3[t,]    <- R
}

reg.SS      <- ifelse(SoF=="SS", 1, 0)
reg.TF      <- ifelse(SoF=="TF", 1, 0)
reg.S       <- ifelse(VS30<=1500, log10(VS30/800), log10(1500/800))

# Create the fd objects by projecting the regressors over an appropriate functional basis

# reg.Mlow
matplot(t.points, reg.Mlow, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=3)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.Mlow, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
Mlow.fd <- smooth.basis(t.points, reg.Mlow, fPar)$fd
plot(Mlow.fd)

# reg.Mhigh
matplot(t.points, reg.Mhigh, type='l')

fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
Mhigh.fd <- smooth.basis(t.points, reg.Mhigh, fPar)$fd
plot(Mhigh.fd)

# R
matplot(t.points, R.mat, type='l', main='True R')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, R.mat, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
R.fd <- smooth.basis(t.points, R.mat, fPar)$fd
plot(R.fd, main='Smoothed R')

# Look what happens for short and long distances: indexes 21 - 36
r <- R.mat[,21]
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, r, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
r.fd <- smooth.basis(t.points, r, fPar)
plot(r.fd)
lines(t.points, r, type='l', xlab='Period', ylab='R - idx 21', col='red')

r <- R.mat[,36]
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- smooth.basis(t.points, r, fPar)$gcv
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
r.fd <- smooth.basis(t.points, r, fPar)
plot(r.fd)
lines(t.points, r, type='l', xlab='Period', ylab='R - idx 21', col='red')

# reg.D1
matplot(t.points, reg.D1, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D1, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
D1.fd <- smooth.basis(t.points, reg.D1, fPar)$fd
plot(D1.fd)

# reg.D2
matplot(t.points, reg.D2, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D2, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
D2.fd <- smooth.basis(t.points, reg.D2, fPar)$fd
plot(D2.fd)

# reg.D3
matplot(t.points, reg.D3, type='l')

basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
esp        <- seq(-1,12, by=1)
lambda.vec <- sort(10^-esp)
gcv.vec    <- numeric(length(lambda.vec))
for(j in 1:length(lambda.vec))
{
  fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
  gcv.vec[j] <- sum(smooth.basis(t.points, reg.D3, fPar)$gcv)/n
}
lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
lambda.opt

fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
D3.fd <- smooth.basis(t.points, reg.D3, fPar)$fd
plot(D3.fd)

# Project constant coefficients over a constant basis
basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=1)

reg.SS <- t(replicate(N, reg.SS))
SS.fd <- smooth.basis(t.points, reg.SS, fPar)$fd
matplot(t.points, reg.SS, type='l')
plot(SS.fd)

reg.TF <- t(replicate(N, reg.TF))
TF.fd <- smooth.basis(t.points, reg.TF, fPar)$fd

reg.S <- t(replicate(N, reg.S))
S.fd <- smooth.basis(t.points, reg.S, fPar)$fd
matplot(t.points, reg.S, type='l')
plot(S.fd)

# Create the intercept
intercept <- t(replicate(N, rep(1,n)))
intercept.fd <- smooth.basis(t.points, intercept, basis)$fd

## Finally build the list of functional regressors
xlist <- list(intercept.fd, Mlow.fd, Mhigh.fd, SS.fd, TF.fd, D1.fd, D2.fd, D3.fd, S.fd)
save(xlist, file='DATA/xlist-log.RData')

