Site <- function(my.dir, name_dir, data, t.idx)
{
  library(latex2exp)
  library(wesanderson)
  
  pal <- wes_palette('Cavalcanti1')
  # Directory where to save plots ---------------------------------------
  dir.current <- getwd()
  
  # DATA and Utilities ----------------------------------------------------------
  load(paste0("Results/",name_dir,"/Regression.RData"))
  
  xlist     <- input$xlist
  t.points  <- input$t.points
  n         <- dim(xlist[[1]]$coefs)[2]
  q         <- length(xlist)
  y         <- mod.fit$yfdobj
  y.hat     <- mod.fit$yhatfdobj

  load('DATA/ITA18_regressors.RData')
  Mh.vec    <- ITA18.regressors$Mh.vec
  h.vec     <- ITA18.regressors$h.vec
  Mref.vec  <- ITA18.regressors$Mref.vec
  M.h.T     <- Mh.vec[t.idx]
  h.T       <- h.vec[t.idx]
  Mref.T    <- Mref.vec[t.idx]
  
  dJB       <- data$dJB
  MAG       <- data$MAG
  SoF       <- data$SoF
  VS30      <- data$VS30
  
  y.T       <- eval.fd(t.points[t.idx], y)
  
  # Analysis of the source term ------------------------------------------------
  a0.est <- mod.fit$betaestlist[[1]]$fd
  b1.est <- mod.fit$betaestlist[[2]]$fd
  b2.est <- mod.fit$betaestlist[[3]]$fd
  f1.est <- mod.fit$betaestlist[[4]]$fd
  f2.est <- mod.fit$betaestlist[[5]]$fd
  c1.est <- mod.fit$betaestlist[[6]]$fd
  c2.est <- mod.fit$betaestlist[[7]]$fd
  c3.est <- mod.fit$betaestlist[[8]]$fd
  k0.est <- mod.fit$betaestlist[[9]]$fd
  
  a0.est.T <- eval.fd(t.points[t.idx], a0.est)
  b1.est.T <- eval.fd(t.points[t.idx], b1.est)
  b2.est.T <- eval.fd(t.points[t.idx], b2.est)
  f1.est.T <- eval.fd(t.points[t.idx], f1.est)
  f2.est.T <- eval.fd(t.points[t.idx], f2.est)
  c1.est.T <- eval.fd(t.points[t.idx], c1.est)
  c2.est.T <- eval.fd(t.points[t.idx], c2.est)
  c3.est.T <- eval.fd(t.points[t.idx], c3.est)
  k0.est.T <- eval.fd(t.points[t.idx], k0.est)
  beta.T   <- as.matrix(c(a0.est.T, b1.est.T, b2.est.T, f1.est.T, f2.est.T, c1.est.T, c2.est.T, c3.est.T, k0.est.T))
  
  ## Construct design matrix at time t -----------------------------------------
  X <- matrix(data=0, nrow=n, ncol=q)
  for(j in 1:q)
  {
    X[,j] <- eval.fd(t.points[t.idx],xlist[[j]])
  }
  
  ## FVS30 ---------------------------------------------------------------------
  idxs <- c(1,2,3,4,5,6,7,8)
  X_red <- as.matrix(X[,idxs])
  B_red <- as.matrix(beta.T[idxs,])
  FS.points <- t(y.T) - X_red%*%B_red
  
  VVS <- seq(range(X_red)[1], range(X_red)[2], length.out=100)

  t <- t.points[t.idx]
  t.plot <- round(10^t, digits=2)
  
  k0.T <- eval.fd(t,k0.est)
  
  FS <- k0.T*VVS
  X.ax <- VVS + log10(800)
  
  png(file = paste0(my.dir,"/site_t",t.idx,".png"), width = 5000, height = 5000, units = "px", res=800)
  par(mar=c(4.2,5.4,1,1)+.1)
  plot(log10(VS30), FS.points, pch=1, col="grey", xaxt='n', ylim=c(-1.5,2), xlab=TeX("$V_{S30}$"),ylab=TeX("$F_{VS30}$"),
       cex.lab=2.2, cex.axis=2.2)
  title(main=paste0("T = ",t.plot,"s"), font.main=1, cex.main=2, line=-2)
  lines(X.ax, FS, lwd=3, col=pal[2])
  xtick<-seq(2, 3.5, by=0.5)
  axis(side=1, at=xtick, labels = round(10^xtick, digits=1), cex.axis=1.8)
  grid()
  dev.off()
  
}