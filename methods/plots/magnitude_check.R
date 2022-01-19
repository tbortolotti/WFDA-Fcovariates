# Magnitude saturation analysis -----------------------------------------

magnitude_check <- function(mod.fit, t.points, name_dir, t.idxs, magnitudes,
                            data, curves)
{
  library(fda)
  library(coda)
  library(devtools)
  library(pbmcapply)
  library(wesanderson)
  
  ## Plot utilities that depend on magnitudes -------------------------
  
  pal    <- wes_palette('Cavalcanti1')
  #cols   <- c(length(magnitudes))
  fondo  <- c(length(magnitudes))
  leg    <- c(length(magnitudes))
  for(m in 1:length(magnitudes))
  {
    val      <- magnitudes[m]
    #cols[m] <- paste0('gray',10*val)
    fondo[m] <- adjustcolor(pal[m], alpha=0.3)
    leg[m]   <- paste0("Data, Mag = ", val)
  }
  pchs   <- rep(16, length(magnitudes))

  
  ## Directory specification
  dir.current <- getwd()
  new.dir <- paste0(dir.current,"/Results/",name_dir,"/magnitude_check")
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
  
  load('DATA/ITA18_regressors.RData')
  Mh.vec   <- ITA18.regressors$Mh.vec
  Mref.vec <- ITA18.regressors$Mref.vec
  h.vec    <- ITA18.regressors$h.vec
  
  # DATA ------------------------------------------------------------------
  dJB       <- data$dJB
  MAG       <- data$MAG
  SoF       <- data$SoF
  VS30      <- data$VS30

  # Repeat the plot for values ----------------------------------------
  vs.vec <- c(450)
  sof.vec <- c('NF','SS','TF')
  
  # for every period, repeat
  for(i in 1:length(t.idxs)) #i=1
  {
    t.idx  <- t.idxs[i]
    t      <- t.points[t.idx]
    t.plot <- 10^t
    
    M.h.T   <- Mh.vec[i]
    M.ref.T <- Mref.vec[i]
    h.T    <- h.vec[i]
    
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
    
    ## VS -------------------------------------------------------------
    for(j in 1:length(vs.vec)) #j=1
    {
      
      vs <- vs.vec[j]
      
      for(k in 1:length(sof.vec)) #k=1
      {
        
        # FS = NF -------------------------------------------------------------------
        sof <- sof.vec[k]
        
        for(m in 1:length(magnitudes)) #m=1
        {
          Mw <- magnitudes[m]
          include <- which(MAG>=(Mw - 0.3) & MAG<=(Mw + 0.3) & SoF==sof)
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
          
          # Plot
          log10dJB        <- log10(dJB)
          log10dJB[1]     <- 0
          log10dJB.seq    <- log10(dJB.seq)
          log10dJB.seq[1] <- 0
          
          if(m==1)
          {
            # open plot
            setwd(new.dir)
            png(file = paste0("VS",vs,"_",sof,"_T",t.plot,".png"), width = 5000, height = 5000, units = "px", res=800)
            par(mar=c(4.2,5.4,1,1)+.1)
            plot(log10dJB[include], SA.points[t.idx,], pch=16, col=fondo[m],
                 xlim=c(0,max(log10dJB[include])), ylim=c(-2,4), xlab="Distance [km]",
                 ylab=paste0("SA [cm/s2]"), xaxt='n', yaxt='n',cex.lab=2.2)
            title(main=paste0('T = ', t.plot," s"), font.main=1, cex.main=2, line=-2)
            
          } else {
            points(log10(dJB[include]), SA.points[t.idx,], pch=16, col=fondo[m])
          }
          
          if(!(length(include)==0))
          {
            lines(log10dJB.seq, SA, lwd=3, col=pal[m])
          }
          
          if(m==length(magnitudes))
          {
            # close plot
            grid()
            legend(0, 0, legend=leg, col=pal, pch=pchs, cex=1.4)
            xtick <- seq(0, 2, by=1)
            ytick <- seq(-2,4)
            axis(side=1, at=xtick, labels = 10^xtick, cex.axis=2.2)
            axis(side=2, at=ytick, labels = 10^ytick, cex.axis=2.2) 
            dev.off()
            setwd(dir.current)
          }
          
        }

      }
      
    }
    
  }
  
}
