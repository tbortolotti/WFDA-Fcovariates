
model_comparison <- function(mod.fit, t.points, name_dir, t.idxs, data, curves, corrective=FALSE,
                             set.log=FALSE)
{
  library(fda)
  library(coda)
  library(devtools)
  library(wesanderson)
  
  pal <- wes_palette('Cavalcanti1')
  fondo1 <- 'grey80'
  fondo2 <- 'grey50'
  
  dir.current <- getwd()
  new.dir <- paste0(dir.current,"/Results/",name_dir,"/ITA18_comparison")
  dir.create(new.dir)

  # Fitted regressors --------------------------------------------
  a0.est <- mod.fit$betaestlist[[1]]$fd
  b1.est <- mod.fit$betaestlist[[2]]$fd
  b2.est <- mod.fit$betaestlist[[3]]$fd
  f1.est <- mod.fit$betaestlist[[4]]$fd
  f2.est <- mod.fit$betaestlist[[5]]$fd
  c1.est <- mod.fit$betaestlist[[6]]$fd
  c2.est <- mod.fit$betaestlist[[7]]$fd
  c3.est <- mod.fit$betaestlist[[8]]$fd
  k0.est <- mod.fit$betaestlist[[9]]$fd
  
  # DATA ------------------------------------------------------------------
  dJB       <- data$dJB
  MAG       <- data$MAG
  SoF       <- data$SoF
  VS30      <- data$VS30
  
  # Fit ITA18

  # ITA18 regressors
  load("DATA/ITA18_coefs.RData")
  a0       <- coefs.ITA18$a0
  b1       <- coefs.ITA18$b1
  b2       <- coefs.ITA18$b2
  c1       <- coefs.ITA18$c1
  c2       <- coefs.ITA18$c2
  c3       <- coefs.ITA18$c3
  k0       <- coefs.ITA18$k0
  f1       <- coefs.ITA18$f1
  f2       <- coefs.ITA18$f2
  
  load('DATA/ITA18_regressors.RData')
  Mh.vec   <- ITA18.regressors$Mh.vec
  Mref.vec <- ITA18.regressors$Mref.vec
  h.vec    <- ITA18.regressors$h.vec
  
  # ITA18 corrective terms -------------------------------------------
  if(corrective)
  {
    load('DATA/ITA18_corrected.RData')
    a.corr    <- ITA18.corrected$a
    b.corr    <- ITA18.corrected$b1
    c.corr    <- ITA18.corrected$c
    f1.corr   <- ITA18.corrected$f1
    f2.corr   <- ITA18.corrected$f2
    h.corr    <- ITA18.corrected$h
    Mref.corr <- ITA18.corrected$Mref
  }

  # Repeat the plot for values ----------------------------------------
  #vs.vec <- c(600,450,300)
  vs.vec <- 450
  sof.vec <- c('NF','SS','TF')
  
  # for every period, repeat
  for(i in 1:length(t.idxs)) #i=1
  {

    t.idx <- t.idxs[i]
    t <- t.points[t.idx]
    if(set.log)
    {
      t.plot <- 10^t
    } else {
      t.plot <- t
    }
    
    
    a0.T <- eval.fd(t,a0.est)
    b1.T <- eval.fd(t,b1.est)
    b2.T <- eval.fd(t,b2.est)
    f1.T <- eval.fd(t,f1.est)
    f2.T <- eval.fd(t,f2.est)
    c1.T <- eval.fd(t,c1.est)
    c2.T <- eval.fd(t,c2.est)
    c3.T <- eval.fd(t,c3.est)
    k0.T <- eval.fd(t,k0.est)
    
    beta.T <- as.matrix(c(a0.T, b1.T, b2.T, f1.T, f2.T, c1.T, c2.T, c3.T, k0.T))
    
    a0.ITA18 <- a0[t.idx]
    b1.ITA18 <- b1[t.idx]
    b2.ITA18 <- b2[t.idx]
    c1.ITA18 <- c1[t.idx]
    c2.ITA18 <- c2[t.idx]
    c3.ITA18 <- c3[t.idx]
    k0.ITA18 <- k0[t.idx]
    f1.ITA18 <- f1[t.idx]
    f2.ITA18 <- f2[t.idx]
    
    beta.T.ITA18 <- as.matrix(c(a0.ITA18, b1.ITA18, b2.ITA18, f1.ITA18,
                                f2.ITA18, c1.ITA18, c2.ITA18, c3.ITA18, k0.ITA18))
    
    M.h.T    <- Mh.vec[t.idx]
    M.ref.T  <- Mref.vec[t.idx]
    h.T      <- h.vec[t.idx]
    
    if(corrective)
    {
      a.corr.T    <- a.corr[t.idx]
      b.corr.T    <- b.corr[t.idx]
      c.corr.T    <- c.corr[t.idx]
      f1.corr.T   <- f1.corr[t.idx]
      f2.corr.T   <- f2.corr[t.idx]
      h.corr.T    <- h.corr[t.idx]
      Mref.corr.T <- Mref.corr[t.idx]
    }
    
    ## VS -------------------------------------------------------------
    for(j in 1:length(vs.vec)) #j=1
    {

      vs <- vs.vec[j]
      
      for(k in 1:length(sof.vec)) #k=1
      {
      
        # FS = NF -------------------------------------------------------------------
        sof <- sof.vec[k]
        
        # 1. Mw = 4.0
        Mw <- 4
        include <- which(MAG>=3.70 & MAG<=4.30 & SoF==sof)
        SA.points <- curves[,include]
        
        # generate dataset for the functional model
        dJB.seq  <- seq(range(dJB)[1],range(dJB)[2],length.out=100)
        nn       <- length(dJB.seq)
        MAG.seq  <- rep(Mw,nn)
        SoF.seq  <- rep(sof,nn)
        VS30.seq <- rep(vs,nn)
        
        R           <- sqrt(dJB.seq^2 + h.T^2)
        intercept   <- rep(1,nn)
        reg.M_l     <- ifelse(MAG.seq<=M.h.T, MAG.seq - M.h.T, 0)
        reg.M_h     <- ifelse(MAG.seq>M.h.T, MAG.seq - M.h.T, 0)
        reg.SS      <- ifelse(SoF.seq=="SS", 1, 0)
        reg.TF      <- ifelse(SoF.seq=="TF", 1, 0)
        reg.D_1     <- (Mw - M.ref.T)*log10(R)
        reg.D_2     <- log10(R)
        reg.D_3     <- R
        reg.S       <- ifelse(VS30.seq<=1500, log10(VS30.seq/800), log10(1500/800))
        XX          <- cbind(intercept, reg.M_l, reg.M_h, reg.SS, reg.TF, reg.D_1, reg.D_2, reg.D_3, reg.S)
        
        SA <- XX%*%beta.T
        SA.ITA18 <- XX%*%beta.T.ITA18
        
        # corrected prediction
        if(corrective == TRUE)
        {
         reg.M <- ifelse(MAG.seq<=Mref.corr.T, MAG.seq - M.h.T, 0)
         delta <- a.corr.T + b.corr.T*reg.M + f1.corr.T*reg.SS + f2.corr.T*reg.TF +
                  c.corr.T*log10(R^2 + h.corr.T^2)/2
         SA.ITA18.corr <- SA.ITA18 + delta
        }
        
        # Plot Mw = 4
        log10dJB        <- log10(dJB)
        log10dJB[1]     <- 0
        log10dJB.seq    <- log10(dJB.seq)
        log10dJB.seq[1] <- 0
        
        png(file = paste0(new.dir,"/VS",vs,"_",sof,"_T",t.plot,".png"), width = 5000, height = 5000, units = "px", res=800)
        par(mar=c(4.2,5.4,1.2,1)+.1)
        plot(log10dJB[include], SA.points[t.idx,], pch=16, col=fondo1,
             xlim=c(0,max(log10dJB[include])), ylim=c(-2,4), xlab="Distance [km]",
             ylab=paste0("SA [cm/s2]"), xaxt='n', yaxt='n', cex.lab=2)
        title(main=paste0("T=",t.plot,"s, SOF= ",sof), line=-2, font.main=1, cex.main=2)
        lines(log10dJB.seq, SA, lwd=4, col=pal[2])
        lines(log10dJB.seq, SA.ITA18, lty=5, lwd=4, col=pal[1])
        if(corrective)
        {
          lines(log10dJB.seq, SA.ITA18.corr, lty=3, lwd=3, col=pal[3])
        }
        grid()
  
        # 2. Mw = 6.8
        Mw = 6.8
        
        include <- which(MAG>=6.50 & MAG<=7.10 & SoF==sof)
        SA.points <- curves[,include]
        
        # generate dataset for the functional model
        dJB.seq  <- seq(range(dJB)[1],range(dJB)[2],length.out=100)
        nn       <- length(dJB.seq)
        MAG.seq  <- rep(Mw,nn)
        SoF.seq  <- rep(sof,nn)
        VS30.seq <- rep(vs,nn)
        
        R           <- sqrt(dJB.seq^2 + h.T^2)
        intercept   <- rep(1,nn)
        reg.M_l     <- ifelse(MAG.seq<=M.h.T, MAG.seq - M.h.T, 0)
        reg.M_h     <- ifelse(MAG.seq>M.h.T, MAG.seq - M.h.T, 0)
        reg.SS      <- ifelse(SoF.seq=="SS", 1, 0)
        reg.TF      <- ifelse(SoF.seq=="TF", 1, 0)
        reg.D_1     <- (Mw - M.ref.T)*log10(R)
        reg.D_2     <- log10(R)
        reg.D_3     <- R
        reg.S       <- ifelse(VS30.seq<=1500, log10(VS30.seq/800), log10(1500/800))
        XX          <- cbind(intercept, reg.M_l, reg.M_h, reg.SS, reg.TF, reg.D_1, reg.D_2, reg.D_3, reg.S)
        
        SA <- XX%*%beta.T
        SA.ITA18 <- XX%*%beta.T.ITA18
        
        # corrected prediction
        if(corrective)
        {
          delta <- a.corr.T + f1.corr.T*reg.SS + f2.corr.T*reg.TF +
                   c.corr.T*log10(R^2 + h.corr.T^2)/2
          SA.ITA18.corr <- SA.ITA18 + delta
        }
        
        
        # Plots Mw = 6.8
        
        points(log10(dJB[include]), SA.points[t.idx,], pch=16, col=fondo2)
        lines(log10dJB.seq, SA, lwd=4, col=pal[2])
        lines(log10dJB.seq, SA.ITA18, lty=5, lwd=4, col=pal[1])
        if(corrective)
        {
          lines(log10dJB.seq, SA.ITA18.corr, lty=3, lwd=4, col=pal[3])
          legend(0, -0.3, legend=c("Functional", "Scalar","Corrected", TeX("Data, $M_w= 4.0 \\pm 0.3$"), TeX("Data, $M_w= 6.8 \\pm 0.3$")),
                 col=c(pal[2], pal[5], pal[5], fondo1, fondo2), lty=c(1,6,3,0,0), pch=c(NA,NA,NA,16,16), cex=1.4)
        } else {
          legend(0, -0.3, legend=c("Functional", "Scalar", TeX("Data, $M_w= 4.0 \\pm 0.3$"), TeX("Data, $M_w= 6.8 \\pm 0.3$")),
                 col=c(pal[2], pal[1], fondo1, fondo2), lwd=c(3,3,0,0,0), lty=c(1,5,0,0), pch=c(NA,NA,16,16), cex=1.2)
        }
        
        xtick <- seq(0, 2, by=1)
        ytick <- seq(-2,4)
        axis(side=1, at=xtick, labels = 10^xtick,cex.axis=2)
        axis(side=2, at=ytick, labels = 10^ytick,cex.axis=2) 
        dev.off()
 
      }

    }
  
  }
  
}
