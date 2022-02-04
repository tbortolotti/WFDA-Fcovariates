

save_plots <- function(mod.fit, input, name_dir)
{

  library(ggplot2)
  library(ggpubr)
  
  load("DATA/ITA18_regressors.RData")
  
  # create directory to save plots
  dir.current <- getwd()
  new.dir <- paste0(dir.current,"/Results/",name_dir)
  dir.create(new.dir)
  
  # save RData ----------------------------------------------------------------------
  save(mod.fit, input, file = paste0(new.dir,"/Regression.RData"))
  
  xlist <- input$xlist
  blist <- input$blist
  t.points <- input$t.points
  basis <- mod.fit$yhatfdobj$basis
  theta.t <- eval.basis(t.points, basis)
  theta <- t(theta.t)
  n <- dim(xlist[[1]]$coefs)[2]
  
  # Fitted values  --------------------------------------------
  y.hat     <- mod.fit$yhatfdobj
  y         <- mod.fit$yfdobj
  rangeval  <- y$basis$rangeval
  res       <- y - y.hat
  y.vals    <- theta.t %*% y$coefs
  yhat.vals <- theta.t %*% y.hat$coefs
  
  a0.est <- mod.fit$betaestlist[[1]]$fd
  b1.est <- mod.fit$betaestlist[[2]]$fd
  b2.est <- mod.fit$betaestlist[[3]]$fd
  f1.est <- mod.fit$betaestlist[[4]]$fd
  f2.est <- mod.fit$betaestlist[[5]]$fd
  c1.est <- mod.fit$betaestlist[[6]]$fd
  c2.est <- mod.fit$betaestlist[[7]]$fd
  c3.est <- mod.fit$betaestlist[[8]]$fd
  k0.est <- mod.fit$betaestlist[[9]]$fd
  
  # Plots ----------------------------------------------------------------------
  
  # True curves
  png(file = paste0(new.dir,"/true_curves.png"), width = 8000, height = 5000, units = "px", res=800)
  matplot(t.points, y.vals, type='l')
  title(main="True curves")
  dev.off()
  
  # Fitted curves
  png(file = paste0(new.dir,"/fitted_curves.png"), width = 8000, height = 5000, units = "px", res=800)
  matplot(t.points, yhat.vals, type='l')
  title(main="Fitted curves")
  dev.off()
  
  # Estimated Regressors
  png(file = paste0(new.dir,"/regressors.png"), width = 8000, height = 5000, units = "px", res=800)
  par(mfrow=c(3,3))
  plot(a0.est, t.points, type='l', lwd=3, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated A")
  plot(b1.est, t.points, type='l', lwd=3, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated B1")
  plot(b2.est, t.points, type='l', lwd=3, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated B2")
  plot(f1.est, t.points, type='l', lwd=3, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated F1")
  plot(f2.est, t.points, type='l', lwd=3, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated F2")
  plot(c1.est, t.points, type='l', lwd=3, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated C1")
  plot(c2.est, t.points, type='l', lwd=3, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated C2")
  plot(c3.est, t.points, type='l', lwd=3, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated C3")
  plot(k0.est, t.points, type='l', lwd=2, col='forestgreen')
  abline(h=0, lty=2, lwd=1, col="grey")
  title(main="Estimated K")
  dev.off()
  
  # ITA18 Comparison
  
  # a0 <- ITA18.regressors[,1]
  # b1 <- ITA18.regressors[,2]
  # b2 <- ITA18.regressors[,3]
  # c1 <- ITA18.regressors[,4]
  # c2 <- ITA18.regressors[,5]
  # c3 <- ITA18.regressors[,6]
  # k0 <- ITA18.regressors[,7]
  # f1 <- ITA18.regressors[,8]
  # f2 <- ITA18.regressors[,9]
  load('DATA/ITA18_coefs.RData')
  a0 <- coefs.ITA18$a0
  b1 <- coefs.ITA18$b1
  b2 <- coefs.ITA18$b2
  c1 <- coefs.ITA18$c1
  c2 <- coefs.ITA18$c2
  c3 <- coefs.ITA18$c3
  k0 <- coefs.ITA18$k0
  f1 <- coefs.ITA18$f1
  f2 <- coefs.ITA18$f2
  
  a0.est.vals <- eval.basis(t.points, a0.est$basis) %*% a0.est$coefs
  b1.est.vals <- eval.basis(t.points, b1.est$basis) %*% b1.est$coefs
  b2.est.vals <- eval.basis(t.points, b2.est$basis) %*% b2.est$coefs
  c1.est.vals <- eval.basis(t.points, c1.est$basis) %*% c1.est$coefs
  c2.est.vals <- eval.basis(t.points, c2.est$basis) %*% c2.est$coefs
  c3.est.vals <- eval.basis(t.points, c3.est$basis) %*% c3.est$coefs
  k0.est.vals <- eval.basis(t.points, k0.est$basis) %*% k0.est$coefs
  f1.est.vals <- eval.basis(t.points, f1.est$basis) %*% f1.est$coefs
  f2.est.vals <- eval.basis(t.points, f2.est$basis) %*% f2.est$coefs
  
  png(file = paste0(new.dir,"/ITA18comparison.png"), width = 8000, height = 5000, units = "px", res=800)
  par(mfrow=c(3,3))
  matplot(t.points, a0.est.vals, type='l', lwd=2, col='forestgreen', main="A")
  points(t.points, a0, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  matplot(t.points, b1.est.vals, type='l', lwd=2, ylim=c(0,1), col='forestgreen', main="B1")
  points(t.points, b1, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  matplot(t.points, b2.est.vals, type='l', lwd=2, ylim=c(-0.2,0.5), col='forestgreen', main="B2")
  points(t.points, b2, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  matplot(t.points, f1.est.vals, type='l', lwd=2, ylim=c(-0.1,0.3), col='forestgreen', main="F1")
  points(t.points, f1, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  matplot(t.points, f2.est.vals, type='l', lwd=2, col='forestgreen', main="F2")
  points(t.points, f2, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  matplot(t.points, c1.est.vals, type='l', lwd=2, ylim=c(0.1,0.5), col='forestgreen', main="C1")
  points(t.points, c1, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  matplot(t.points, c2.est.vals, type='l', lwd=2, ylim=c(-1.7,-1.1), col='forestgreen', main="C2")
  points(t.points, c2, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  matplot(t.points, c3.est.vals, type='l', lwd=2, ylim=c(-0.01,0.01), col='forestgreen', main="C3")
  points(t.points, c3, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  matplot(t.points, k0.est.vals, type='l', lwd=2, ylim=c(-0.9,-0.2), col='forestgreen', main="K")
  points(t.points, k0, type='l', lwd=2, col='darkolivegreen3')
  abline(h=0, lty=2, lwd=1, col="grey")
  dev.off()

  # Residuals
  png(file = paste0(new.dir,"/residuals.png"), width = 8000, height = 5000, units = "px", res=800)
  plot(res, t.points)
  title(main="Residuals")
  dev.off()
  
  
  # Squared multiple correlation function -------------------------------
  # Squared residuals
  c.hat <- y.hat$coefs
  c.y <- y$coefs
  c.res2 <- (c.y-c.hat) * (c.y-c.hat)
  
  # SSE - sum of squared residuals
  SSE <- rowSums(c.res2)
  SSE.manual <- theta.t%*%SSE
  
  # SSY
  c.a0 <- a0.est$coefs
  c.A0 <- c.a0%*%t(rep(1,n))
  ssy <- (c.y - c.A0) * (c.y - c.A0)
  SSY <- rowSums(ssy)
  SSY.manual <- theta.t%*%SSY
  
  RSQ.t <- (SSY.manual - SSE.manual)/SSY.manual
  png(file = paste0(new.dir,"/RSQ.png"), width = 8000, height = 5000, units = "px", res=800)
  matplot(t.points, RSQ.t, lwd=2, col="darkblue", ylim=c(0.8,1),
          type='l', main="Squared multiple correlation function")
  abline(h=1, col="lightblue")
  dev.off()
  
  # F-ratio ------------------------------------------------------------
  F.stat <- Fstat.fd(y,y.hat,argvals=NULL)$F

  p <- 9
  alpha <- 0.05
  F.alpha <- qf(1-alpha, p-1, n-p)
  
  png(file = paste0(new.dir,"/F_statistic.png"), width = 8000, height = 5000, units = "px", res=800)
  plot(F.stat, type='l', lwd=2, ylim=c(0,16), col="darkblue", main="Point-wise F statistic")
  abline(h=F.alpha, lty='dashed', col="grey")
  dev.off()

}