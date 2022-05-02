## Plot utili risultati della scelta del parametro che definisce il peso

## MSE point-wise: metodo pesato, non pesato, con peso che cade a 0 -------------

load('Results/a-definition/MSE_5fold_fullstats_nowgts.RData')
N <- 37
MSE.mat <- matrix(data=0, nrow=N, ncol=1)
MSE.mat[,1] <- apply(MSE_cv[,,1], MARGIN=1, FUN=mean)

load('Results/a-definition/MSE_5fold_fullstats_lowvals.RData')
temp <- MSE.mat
MSE.mat <- matrix(data=0, nrow=N, ncol=4)
MSE.mat[,1] <- temp
for(a in 1:3)
{
  MSE.mat[,(1+a)] <- apply(MSE_cv[,,a+2], MARGIN=1, FUN=mean)
}

load('Results/a-definition/MSE_5fold_fullstats_newwgts.RData')
temp <- MSE.mat
curr.dim <- dim(temp)[2]
MSE.mat <- matrix(data=0, nrow=N, ncol=(curr.dim+1))
MSE.mat[,1:curr.dim] <- temp
MSE.mat[,(curr.dim+1)] <- apply(MSE_cv[,,3], MARGIN=1, FUN=mean)

load('Results/a-definition/MSE_5fold_fullstats_newwgts_largevals.RData')
temp <- MSE.mat
curr.dim <- dim(temp)[2]
MSE.mat <- matrix(data=0, nrow=N, ncol=(curr.dim+1))
MSE.mat[,1:curr.dim] <- temp
MSE.mat[,(curr.dim+1)] <- apply(MSE_cv[,,2], MARGIN=1, FUN=mean)


## Plot
y.limits <- c(min(MSE.mat), max(MSE.mat))
plot(t.points, MSE.mat[,1], type='l', lwd=2, col=1, ylim=y.limits, xlab='log-Period', ylab="Pointwise MSE")
for(a in 2:dim(MSE.mat)[2])
{
  lines(t.points, MSE.mat[,a], type='l', lwd=2, col=a)
}
legend(x=0, y=0.18, legend=c("no-wgts", "a=1", "a=2", "a=4","a=10", "0-wgts"), col=c(1,2,3,4,5,6), lty=1, lwd=2)

## Zoom the segment [0,1]
idxs <- seq(21,37,by=1)
y.limits <- c(min(MSE.mat[idxs,]), max(MSE.mat[idxs,]))
plot(t.points[idxs], MSE.mat[idxs,1], type='l', lwd=3, col=1, ylim=y.limits, xlab='log-Period', ylab="Pointwise MSE")
for(a in 2:dim(MSE.mat)[2])
{
  lines(t.points[idxs], MSE.mat[idxs,a], type='l', lwd=3, col=a)
}
legend(x=0, y=0.095, legend=c("no-wgts","a=1", "a=2", "a=4","a=10", "0-wgts"), col=c(1,2,3,4,5,6), lty=1, lwd=2)


y.limits <- c(min(MSE.mat), max(MSE.mat))
plot(T.period, MSE.mat[,1], type='l', lwd=3, col=1,  ylim=y.limits, xlab='Period', ylab="Pointwise MSE")
for(a in 2:dim(MSE.mat)[2])
{
  lines(T.period, MSE.mat[,a], type='l', lwd=3, col=a)
}
legend(x=8, y=0.18, legend=c("no-wgts","a=1", "a=2", "a=4","a=10", "0-wgts"), col=c(1,2,3,4,5,6), lty=1, lwd=2)

## Zoom the segment [6,10]
idxs <- seq(21,37,by=1)
y.limits <- c(min(MSE.mat[idxs,]), max(MSE.mat[idxs,]))
plot(T.period[idxs], MSE.mat[idxs,1], type='l', lwd=3, col=1, ylim=y.limits, xlab='Period', ylab="Pointwise MSE")
for(a in 2:dim(MSE.mat)[2])
{
  lines(T.period[idxs], MSE.mat[idxs,a], type='l', lwd=3, col=a)
}
legend(x=2, y=0.095, legend=c("no-wgts", "a=1", "a=2", "a=4","a=10", "0-wgts"), col=c(1,2,3,4,5,6), lty=1, lwd=2)



## Boxplots of the MSE functional ------------------------------------------------
load('Results/a-definition/MSE_5fold_fullstats_lowvals.RData')
B <- dim(MSE_fun)[1]

MSE_vals <- c(MSE_fun[,3], MSE_fun[,4], MSE_fun[,5])

load('Results/a-definition/MSE_5fold_fullstats_nowgts.RData')
MSE_vals <- c(MSE_vals, MSE_fun[,1])

load('Results/a-definition/MSE_5fold_fullstats_0wgts.RData')
MSE_vals <- c(MSE_vals, MSE_fun[,1])

load('Results/a-definition/MSE_5fold_fullstats_newwgts.RData')
MSE_vals <- c(MSE_vals, MSE_fun[,3])

class <- c(rep('a=1', B),
           rep('a=2', B),
           rep('a=4', B),
           rep('no-wgts', B),
           rep('0-wgts', B),
           rep('a=10', B))

data.box <- data.frame(MSE_vals, as.factor(class))
names(data.box) <- c('MSE', 'class')

pg_plot <- ggplot(data.box, aes(x = class, y = MSE))+
  geom_boxplot(color="black", fill=tim.colors(6))
pg_plot +
  scale_x_discrete(limits = levels(data.box$class)) +
  #scale_y_continuous(limits=c(0.08,0.2)) +
  labs(y="5-fold Mean Squared Error", x="") +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position="None")

## Boxplots for the variance of the functional coefficients estimators ---------
levs <- 5
coefficient <- rep(c("a", "b1", "b2","c1","c2","c3","k"), each=levs*B)
method      <- rep(c("no-wgts", "a=1", "a=2", "a=4", "0-wgts"), each=B)
beta_var   <- numeric(7*levs*B)

# Riempio il beta_var
{
  load('Results/a-definition/MSE_5fold_fullstats_nowgts.RData')
  ## j=1
  j <- 1
  beta1.est <- beta1.list[[1]]
  beta2.est <- beta2.list[[1]]
  beta3.est <- beta3.list[[1]]
  beta6.est <- beta6.list[[1]]
  beta7.est <- beta7.list[[1]]
  beta8.est <- beta8.list[[1]]
  beta9.est <- beta9.list[[1]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta3.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta3.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta6.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta6.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta7.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta7.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta8.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta8.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta9.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta9.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  
  load('Results/a-definition/MSE_5fold_fullstats_lowvals.RData')
  ## j=2
  j <- 2
  beta1.est <- beta1.list[[3]]
  beta2.est <- beta2.list[[3]]
  beta3.est <- beta3.list[[3]]
  beta6.est <- beta6.list[[3]]
  beta7.est <- beta7.list[[3]]
  beta8.est <- beta8.list[[3]]
  beta9.est <- beta9.list[[3]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta3.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta3.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta6.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta6.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta7.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta7.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta8.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta8.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta9.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta9.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  
  ## j=3
  j <- 3
  beta1.est <- beta1.list[[4]]
  beta2.est <- beta2.list[[4]]
  beta3.est <- beta3.list[[4]]
  beta6.est <- beta6.list[[4]]
  beta7.est <- beta7.list[[4]]
  beta8.est <- beta8.list[[4]]
  beta9.est <- beta9.list[[4]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta3.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta3.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta6.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta6.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta7.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta7.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta8.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta8.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta9.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta9.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  ## j=4
  j <- 4
  beta1.est <- beta1.list[[5]]
  beta2.est <- beta2.list[[5]]
  beta3.est <- beta3.list[[5]]
  beta6.est <- beta6.list[[5]]
  beta7.est <- beta7.list[[5]]
  beta8.est <- beta8.list[[5]]
  beta9.est <- beta9.list[[5]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta3.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta3.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta6.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta6.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta7.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta7.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta8.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta8.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta9.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta9.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  
  ## j=5
  load('Results/a-definition/MSE_5fold_fullstats_0wgts.RData')
  j <- 5
  beta1.est <- beta1.list[[1]]
  beta2.est <- beta2.list[[1]]
  beta3.est <- beta3.list[[1]]
  beta6.est <- beta6.list[[1]]
  beta7.est <- beta7.list[[1]]
  beta8.est <- beta8.list[[1]]
  beta9.est <- beta9.list[[1]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 2
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta2.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta2.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 3
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta3.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta3.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 4
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta6.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta6.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 5
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta7.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta7.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 6
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta8.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta8.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  i <- 7
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta9.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta9.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
}
  
# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("no-wgts", "a=1", "a=2", "a=4", "0-wgts")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  #scale_y_continuous(limits=c(0,0.1)) +
  ylab(TeX(r'($||\hat{\beta}^b - \frac{1}{B} \sum_{b=1}^B \hat{\beta^b} ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  theme_bw() +
  labs(fill = "Weights parameter")

## Boxplot per ciascun regressore separatamente -------------------------------------
levs <- 6
coefficient <- rep(c("c3"), each=levs*B)
method      <- rep(c("no-wgts", "a=1", "a=2", "a=4", "a=10" ,"0-wgts"), each=B)
beta_var   <- numeric(1*levs*B)

# choose among c(1,2,3,6,7,8,9)
load('Results/a-definition/MSE_5fold_fullstats_lowvals.RData')
beta.list <- beta8.list

load('Results/a-definition/MSE_5fold_fullstats_nowgts.RData')
betano.list <- beta8.list

load('Results/a-definition/MSE_5fold_fullstats_newwgts_largevals.RData')
beta0.list <- beta8.list

load('Results/a-definition/MSE_5fold_fullstats_newwgts.RData')
beta10.list <- beta8.list

# load('Results/a-definition/MSE_5fold_fullstats_0wgts.RData')
# beta0.list <- beta1.list

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
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  ## j=2
  j <- 2
  beta1.est <- beta.list[[3]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  ## j=3
  j <- 3
  beta1.est <- beta.list[[4]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  ## j=4
  j <- 4
  beta1.est <- beta.list[[5]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  ## j=5
  j <- 5
  beta1.est <- beta10.list[[3]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
  
  ## j=6
  j <- 6
  beta1.est <- beta0.list[[2]]
  
  i <- 1
  idx.inf <- (levs*(i-1) + j-1)*B + 1
  idx.sup <- (levs*(i-1)+j)*B
  idxs <- idx.inf:idx.sup
  mean0 <- mean.fd(beta1.est)
  mean0$coefs <- mean0$coefs%*%rep(1,B)
  beta.diff <- beta1.est - mean0
  beta_var[idxs] <- diag(inprod(beta.diff, beta.diff, 0, 0, rng=c(0,1)))
}

# Costruisco il dataframe per il plot
data <- data.frame(coefficient, method, beta_var)
method_order<- c("no-wgts", "a=1", "a=2", "a=4", "a=10" ,"0-wgts")
data.box <- data %>% mutate(method=factor(x=method, levels=method_order))

##grouped boxplot
ggplot(data.box, aes(x=coefficient, y=beta_var, fill=method)) + 
  geom_boxplot() +
  #scale_y_continuous(limits=c(0,0.1)) +
  ylab(TeX(r'($||\hat{\beta}^b - \frac{1}{B} \sum_{b} \hat{\beta^b} ||_2^2$)')) + 
  xlab("") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size=18),
        legend.text=element_text(size=19)) +
  theme_bw() +
  labs(fill = "Weights parameter")

