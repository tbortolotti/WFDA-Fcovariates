## Plot utili risultati della scelta del parametro che definisce il peso

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
library(beepr)
library(tidyverse)
library(latex2exp)

rm(list=ls())
graphics.off()
cat("\014")

## Utilities -----------------------------------------------------------------
# Logarithm of the period
load('DATA/t_period.RData')
t.points    <- log10(T.period)
t.points[1] <- -2.5
N <- length(t.points)
B <- 10

## MSE POINTWISE: metodo non pesato, pesato, con peso che cade a 0 ------------
load('Results/a-definition/MSE_10fold_nowgts.RData')
MSE.mat <- matrix(data=0, nrow=N, ncol=1)
MSE.mat[,1] <- apply(MSE_cv[,,1], MARGIN=1, FUN=mean)

load('Results/a-definition/MSE_10fold_oldwgts.RData')
for(a in 2:length(vec.par))
{
  MSE.mat <- cbind(MSE.mat,apply(MSE_cv[,,a], MARGIN=1, FUN=mean))
}

load('Results/a-definition/MSE_10fold_oldwgts_hugevals.RData')
MSE.mat <- cbind(MSE.mat,apply(MSE_cv[,,2], MARGIN=1, FUN=mean))


load('Results/a-definition/MSE_10fold_0wgts.RData')
temp <- MSE.mat
MSE.mat <- cbind(temp, apply(MSE_cv[,,1], MARGIN=1, FUN=mean))

## Plot
y.limits <- c(min(MSE.mat), max(MSE.mat))
plot(t.points, MSE.mat[,1], type='l', lwd=2, col=1, ylim=y.limits, xlab='log-Period', ylab="Pointwise MSE")
for(a in 2:dim(MSE.mat)[2])
{
  lines(t.points, MSE.mat[,a], type='l', lwd=2, col=a)
}
legend(x=0, y=0.18, legend=c("no-wgts", "a=5", "a=10","a=15", "a=20", "a=inf" ,"0-wgts"), col=c(1,2,3,4,5,6,7), lty=1, lwd=2)

## Zoom the segment [0,1]
idxs <- seq(21,37,by=1)
y.limits <- c(min(MSE.mat[idxs,]), max(MSE.mat[idxs,]))
plot(t.points[idxs], MSE.mat[idxs,1], type='l', lwd=3, col=1, ylim=y.limits, xlab='log-Period', ylab="Pointwise MSE")
for(a in 2:dim(MSE.mat)[2])
{
  lines(t.points[idxs], MSE.mat[idxs,a], type='l', lwd=3, col=a)
}
legend(x=0, y=0.098, legend=c("no-wgts", "a=5", "a=10","a=15", "a=20", "a=inf" ,"0-wgts"), col=c(1,2,3,4,5,6,7), lty=1, lwd=2)


y.limits <- c(min(MSE.mat), max(MSE.mat))
plot(T.period, MSE.mat[,1], type='l', lwd=3, col=1,  ylim=y.limits, xlab='Period', ylab="Pointwise MSE")
for(a in 2:dim(MSE.mat)[2])
{
  lines(T.period, MSE.mat[,a], type='l', lwd=3, col=a)
}
legend(x=8, y=0.18, legend=c("no-wgts", "a=5", "a=10","a=15", "a=20", "a=inf" ,"0-wgts"), col=c(1,2,3,4,5,6,7), lty=1, lwd=2)

## Zoom the segment [6,10]
idxs <- seq(21,37,by=1)
y.limits <- c(min(MSE.mat[idxs,]), max(MSE.mat[idxs,]))
plot(T.period[idxs], MSE.mat[idxs,1], type='l', lwd=3, col=1, ylim=y.limits, xlab='Period', ylab="Pointwise MSE")
for(a in 2:dim(MSE.mat)[2])
{
  lines(T.period[idxs], MSE.mat[idxs,a], type='l', lwd=3, col=a)
}
legend(x=2, y=0.098, legend=c("no-wgts", "a=5", "a=10","a=15", "a=20", "a=inf" ,"0-wgts"), col=c(1,2,3,4,5,6,7), lty=1, lwd=2)


## BOXPLOT degli MSE valutati curva per curva ----------------------------------
MSE.glob.mat <- matrix(data=0, nrow=B, ncol=1)

load('Results/a-definition/MSE_10fold_nowgts.RData')
MSE.glob.mat[,1] <- unlist(lapply(MSE_glob_bis.list[[1]], FUN=mean))

load('Results/a-definition/MSE_10fold_oldwgts.RData')
for(a in 2:length(vec.par))
{
  MSE.glob.mat <- cbind(MSE.glob.mat, unlist(lapply(MSE_glob_bis.list[[a]], FUN=mean)))
}

load('Results/a-definition/MSE_10fold_oldwgts_hugevals.RData')
MSE.glob.mat <- cbind(MSE.glob.mat, unlist(lapply(MSE_glob_bis.list[[2]], FUN=mean)))

load('Results/a-definition/MSE_10fold_0wgts.RData')
MSE.glob.mat <- cbind(MSE.glob.mat, unlist(lapply(MSE_glob_bis.list[[1]], FUN=mean)))

## Minimo brutale
MSE.mean <- apply(MSE.glob.mat , MARGIN=2, FUN=mean)
MSE.mean
MSE.sd <- apply(MSE.glob.mat , MARGIN=2, FUN=sd)
MSE.sd
which(MSE.mean==min(MSE.mean))

MSE.glob.vec <- vec(MSE.glob.mat)
## Plot
B <- 10
class <- c(rep('no-wgts', B),
           rep('a=5', B),
           rep('a=10', B),
           rep('a=15', B),
           rep('a=20', B),
           rep('a=inf', B),
           rep('0-wgts', B))

data.box <- data.frame(MSE.glob.vec, as.factor(class))
names(data.box) <- c('MSE', 'class')
class_order<- c("no-wgts", "a=5", "a=10","a=15", "a=20", "a=inf", "0-wgts")
data.box <- data.box %>% mutate(class=factor(x=class, levels=class_order))

box.dir <- '/Users/teresabortolotti/Desktop/flash/a-definition'
pg_plot <- ggplot(data.box, aes(x = class, y = MSE))+
  geom_boxplot(color="black", fill=tim.colors(length(class)/B))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  #scale_y_continuous(limits=c(0.08,0.2)) +
  theme(legend.position="None") + 
  theme_bw() +
  labs(x="", y=TeX(r'($\frac{1}{n^b} \; \sum_{i=1}^{n^b} err_i^2$)'), title="Weights definition: 10-fold MSE") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

ggsave(filename = "global-MSE.pdf",
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

## Boxplot della varianza per ciascun regressore separatamente -------------------------------------
levs <- 7
coefficient <- rep(c("k"), each=levs*B)
name.coef <- names(table(coefficient))[1]
method      <- rep(c("no-wgts", "a=5", "a=10", "a=15" ,"a=20","a=inf", "0-wgts"), each=B)
beta_var   <- numeric(1*levs*B)

# choose among c(1,2,3,6,7,8,9)
load('Results/a-definition/MSE_10fold_oldwgts.RData')
beta.list <- beta9.list

load('Results/a-definition/MSE_10fold_oldwgts_hugevals.RData')
betalarge.list <- beta9.list

load('Results/a-definition/MSE_10fold_nowgts.RData')
betano.list <- beta9.list

load('Results/a-definition/MSE_10fold_0wgts.RData')
betazero.list<- beta9.list


# Riempio il beta_var
{
  ## j=1
  j <- 1
  beta1.est <- betano.list[[1]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  ## j=2
  j <- 2
  beta1.est <- beta.list[[2]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  ## j=3
  j <- 3
  beta1.est <- beta.list[[3]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  ## j=4
  j <- 4
  beta1.est <- beta.list[[4]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  ## j=5
  j <- 5
  beta1.est <- beta.list[[5]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  ## j=6
  j <- 6
  beta1.est <- betalarge.list[[2]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
  
  ## j=7
  j <- 7
  beta1.est <- betazero.list[[1]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=range(t.points)))
}

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("no-wgts", "a=5", "a=10", "a=15" ,"a=20","a=inf", "0-wgts")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
box.dir <- '/Users/teresabortolotti/Desktop/flash/a-definition'

ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  #scale_y_continuous(limits=c(0,0.1)) +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b} \hat{\beta^b} ||_2^2$)'),
       title=paste0("Coefficient ", name.coef, " : Variance"), fill="Parameter") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 22)) +
  theme_bw()
  
ggsave(filename = paste0("Var-global-",name.coef,".pdf"),
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


## Boxplot della varianza per ciascun regressore separatamente nell'ultima parte del dominio -------------------------------------
levs <- 7
coefficient <- rep(c("k"), each=levs*B)
name.coef <- names(table(coefficient))[1]
method      <- rep(c("no-wgts", "a=5", "a=10", "a=15" , "a=20", "a=inf", "0-wgts"), each=B)
beta_var   <- numeric(1*levs*B)

# choose among c(1,2,3,6,7,8,9)
load('Results/a-definition/MSE_10fold_oldwgts.RData')
beta.list <- beta9.list

load('Results/a-definition/MSE_10fold_oldwgts_hugevals.RData')
betalarge.list <- beta9.list

load('Results/a-definition/MSE_10fold_nowgts.RData')
betano.list <- beta9.list

load('Results/a-definition/MSE_10fold_0wgts.RData')
betazero.list<- beta9.list

# Riempio il beta_var
{
  ## j=1
  j <- 1
  beta1.est <- betano.list[[1]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0.77,1)))
  
  ## j=2
  j <- 2
  beta1.est <- beta.list[[2]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0.77,1)))
  
  ## j=3
  j <- 3
  beta1.est <- beta.list[[3]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0.77,1)))
  
  ## j=4
  j <- 4
  beta1.est <- beta.list[[4]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0.77,1)))
  
  ## j=5
  j <- 5
  beta1.est <- beta.list[[5]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0.77,1)))
  
  ## j=6
  j <- 6
  beta1.est <- betalarge.list[[2]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0.77,1)))
  
  ## j=7
  j <- 7
  beta1.est <- betazero.list[[1]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0.77,1)))
  
}

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("no-wgts", "a=5", "a=10", "a=15" , "a=20", "a=inf", "0-wgts")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
box.dir <- '/Users/teresabortolotti/Desktop/flash/a-definition'

ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  #scale_y_continuous(limits=c(0,0.1)) +
  labs(x="", y=TeX(r'($||\hat{\beta}^b - \frac{1}{B} \; \sum_{b} \hat{\beta^b} ||_2^2$)'),
       title=paste0("Coefficient ", name.coef, " : Partial Variance"), fill="Parameter") + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 22)) +
  theme_bw()

ggsave(filename = paste0("Var-partial-",name.coef,".pdf"),
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

## RECONSTRUCTION METHODS COMPARISON: EXTRAP vs KLAL ---------------------------

load('Results/rec-method-comparison/MSE_10fold_oldwgts_extrap.RData')
MSE.glob.vec <- unlist(lapply(MSE_glob.list[[1]], FUN=mean))

load('Results/a-definition/MSE_10fold_oldwgts_hugevals.RData')
MSE.glob.vec <- c(MSE.glob.vec, unlist(lapply(MSE_glob.list[[2]], FUN=mean)))

load('Results/rec-method-comparison/MSE_10fold_oldwgts_extrap_nowgts.RData')
MSE.glob.vec <- c(MSE.glob.vec, unlist(lapply(MSE_glob.list[[1]], FUN=mean)))

load('Results/a-definition/MSE_10fold_nowgts.RData')
MSE.glob.vec <- c(MSE.glob.vec, unlist(lapply(MSE_glob.list[[1]], FUN=mean)))


## Plot
B <- 10
class <- c(rep('extrap:wgt', B),
           rep('KL-AL:wgt', B),
           rep('extrap', B),
           rep('KL-AL', B))

data.box <- data.frame(MSE.glob.vec, as.factor(class))
names(data.box) <- c('MSE', 'class')
class_order<- c('extrap:wgt', 'KL-AL:wgt', 'extrap', 'KL-AL')
data.box <- data.box %>% mutate(class=factor(x=class, levels=class_order))

pg_plot <- ggplot(data.box, aes(x = class, y = MSE))+
  geom_boxplot(color="black", fill=tim.colors(length(class)/B))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  #scale_y_continuous(limits=c(0.08,0.2)) +
  labs(y="10-fold Mean Squared Error", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")

## RESULT SECCO ------------------------------------------------------------------
load('Results/rec-method-comparison/MSE_10fold_oldwgts_extrap.RData')
MSE.comp <- unlist(lapply(MSE_glob.list[[1]], FUN=mean))

load('Results/a-definition/MSE_10fold_oldwgts_hugevals.RData')
MSE.comp <- cbind(MSE.comp, unlist(lapply(MSE_glob.list[[2]], FUN=mean)))

load('Results/rec-method-comparison/MSE_10fold_oldwgts_extrap_nowgts.RData')
MSE.comp <- cbind(MSE.comp, unlist(lapply(MSE_glob.list[[1]], FUN=mean)))

load('Results/a-definition/MSE_10fold_nowgts.RData')
MSE.comp <- cbind(MSE.comp, unlist(lapply(MSE_glob.list[[1]], FUN=mean)))

lamorte.mean <- apply(MSE.comp, MARGIN=2, FUN=mean)
lamorte.mean
lamorte.sd <- apply(MSE.comp, MARGIN=2, FUN=sd)
lamorte.sd
which(lamorte.mean==min(lamorte.mean))

## RESULT SECCO BIS ------------------------------------------------------------------
load('Results/rec-method-comparison/MSE_10fold_oldwgts_extrap.RData')
MSE.comp <- unlist(lapply(MSE_glob_bis.list[[1]], FUN=mean))

load('Results/a-definition/MSE_10fold_oldwgts_hugevals.RData')
MSE.comp <- cbind(MSE.comp, unlist(lapply(MSE_glob_bis.list[[2]], FUN=mean)))


lamorte.mean <- apply(MSE.comp, MARGIN=2, FUN=mean)
lamorte.mean
lamorte.sd <- apply(MSE.comp, MARGIN=2, FUN=sd)
lamorte.sd
which(lamorte.mean==min(lamorte.mean))



