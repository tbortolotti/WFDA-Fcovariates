Path_unified <- function(my.dir, name_dir, data, t.idx, set.log=FALSE)
{
  library(latex2exp)
  library(wesanderson)
  
  pal <- wes_palette('Cavalcanti1')
  
  # DATA and Utilities ----------------------------------------------------------
  load(paste0("Results/",name_dir,"/Regression.RData"))
  
  xlist     <- input$xlist
  t.points  <- input$t.points
  n         <- dim(xlist[[1]]$coefs)[2]
  q         <- length(xlist)
  y         <- mod.fit$yfdobj
  y.hat     <- mod.fit$yhatfdobj
  theta.t   <- eval.basis(t.points, y$basis)
  theta     <- t(theta.t)
  
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
  
  ## FD --------------------------------------------------------
  idxs <- c(1,2,3,4,5,9)
  X_red <- as.matrix(X[,idxs])
  B_red <- as.matrix(beta.T[idxs,])
  
  FD.points <- t(y.T) - X_red%*%B_red
  
  ## MAG = 4
  in.mag <- which(MAG>=3.60 & MAG<=4.40)
  
  dJB.ord <- seq(range(dJB)[1], range(dJB)[2], length.out=100)
  
  R           <- sqrt(dJB.ord^2 + h.T^2)
  
  reg.D_1     <- (4 - M.ref.T)*log10(R)
  reg.D_2     <- log10(R)
  reg.D_3     <- R
  XX          <- cbind(reg.D_1, reg.D_2, reg.D_3)
  
  t <- t.points[t.idx]
  if(set.log)
  {
    t.plot <- round(10^t, digits=2)
  } else {
    t.plot <- t
  }
  
  c1.T <- eval.fd(t,c1.est)
  c2.T <- eval.fd(t,c2.est)
  c3.T <- eval.fd(t,c3.est)
  
  FD <- c1.T*reg.D_1 + c2.T*reg.D_2 + c3.T*reg.D_3
  
  log10dJB    <- log10(dJB)
  log10dJB[1] <- 0
  log10dJB.ord    <- log10(dJB.ord)
  log10dJB.ord[1] <- 0
  
  png(file = paste0(my.dir,"/path_t",t.idx,".png"), width = 5000, height = 5000, units = "px", res=800)
  par(mar=c(4.2,5.4,1,1)+.1)
  plot(log10dJB[in.mag],FD.points[in.mag], pch=1, col="gray80",
       xlim=c(0,max(log10dJB[in.mag])), ylim=c(-5.5,0), xaxt='n',
       xlab="Distance [km]", ylab=TeX("$F_D$"), cex.lab=2.2, cex.axis=2.2)
  title(main=paste0("T = ",t.plot,"s"), font.main=1, cex.main=2, line=-2)
  lines(log10dJB.ord, FD, lwd=3, col=pal[2])
  xtick<-seq(0, 2, by=1)
  axis(side=1, at=xtick, labels = 10^xtick, cex.axis=2.2)
  grid()
  
  ## MAG = 6
  in.mag <- which(MAG>=5.60 & MAG<=6.43)
  
  dJB.ord <- seq(range(dJB)[1], range(dJB)[2], length.out=100)
  
  R           <- sqrt(dJB.ord^2 + h.T^2)
  
  reg.D_1     <- (6 - M.ref.T)*log10(R)
  reg.D_2     <- log10(R)
  reg.D_3     <- R
  XX          <- cbind(reg.D_1, reg.D_2, reg.D_3)
  
  c1.T <- eval.fd(t,c1.est)
  c2.T <- eval.fd(t,c2.est)
  c3.T <- eval.fd(t,c3.est)
  
  FD <- c1.T*reg.D_1 + c2.T*reg.D_2 + c3.T*reg.D_3
  
  points(log10dJB[in.mag],FD.points[in.mag], pch=1, col="gray50")
  lines(log10dJB.ord, FD, lwd=3, col=pal[1])
  legend(0, -3.8, legend=c("SA, M = 4.0","SA, M = 6.0", TeX("Data, $M_w= 4.0 \\pm 0.3$"), TeX("Data, $M_w= 6.0 \\pm 0.3$")),
         col=c(pal[2],pal[1], 'gray80', 'gray50'), lty=c(1,1,0,0), pch=c(NA,NA,16,16), cex=1.4)
  
  dev.off()
  
}