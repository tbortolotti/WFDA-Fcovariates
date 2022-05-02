create_old_weights <- function(curves.rec, t.points, breaks, fix.par, reconst_fcts, Thp)
{
  ## Utilities -----------------------------------------------------------------
  n <- dim(curves.rec)[2]
  
  ## Building the weights ------------------------------------------------------
  #tt <- seq(range(t.points)[1], range(t.points)[2], by=0.01)
  tt <- t.points
  
  T.last <- range(t.points)[2]
  wgts.obs <- matrix(data=1, nrow=length(tt), ncol=n)
  
  for(i in 1:length(reconst_fcts)) # i = 1
  {
    t.max    <- Thp[reconst_fcts[i]]
    imax     <- tail(which(t.points<=t.max),n=1)
    imax.new <- tail(which(tt<=t.max),n=1)
    T.max    <- t.points[imax]
    
    scale    <- fix.par*sd(curves.rec[imax,-reconst_fcts])
    loc      <- T.max + (T.last - T.max)/2
    
    # if(set.log)
    # {
    #   loc    <- 10^T.max + loc.par
    #   x      <- 10^tt[(imax.new+1):length(tt)]
    # } else {
    #   loc    <- T.max + loc.par
    #   x      <- tt[(imax.new+1):length(tt)]
    # }
    
    x      <- tt[(imax.new+1):length(tt)]
    
    y.x      <- 1/(1+exp((x - loc)*scale)) + (1 - 1/(1+exp((x[1]-loc)*scale)))
    y        <- rep(1,imax.new)
    y        <- c(y,y.x)
    
    wgts.obs[,reconst_fcts[i]] <- y
  }
  
  # five <- wgts.obs[,reconst_fcts[i]]
  # ten <- wgts.obs[,reconst_fcts[i]]
  # fifteen <- wgts.obs[,reconst_fcts[i]]
  # twenty <- wgts.obs[,reconst_fcts[i]]
  # inf <- wgts.obs[,reconst_fcts[i]]
  # wgts0 <- wgts.obs[,reconst_fcts[i]]
  
  # temp <- which(tt %in% t.points)
  # pdf(file = paste0("useful-pics/weight.pdf"), width = 8, height = 5)
  # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
  # plot(tt, wgts.obs[,reconst_fcts[1]], type='l', lwd=3, col='darkorange',
  #      main="(b) Weight", ylim=c(0,1), xlab="t", ylab="w(t)",
  #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
  # points(tt[temp], wgts.obs[temp,reconst_fcts[1]],pch=16, col='darkorange')
  # grid()
  # dev.off()
  
  # save(five, ten, fifteen, twenty, inf, wgts0, file='useful-pics/info_plot_wgts.RData')
  # load('useful-pics/info_plot_wgts.RData')
  # library(ggplot2)
  # gc.ramp <- hue_pal()(7)
  # temp <- which(tt %in% t.points)
  # pdf(file = paste0("useful-pics/weight.pdf"), width = 8, height = 5)
  # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
  # plot(tt, ten, type='l', lwd=3, lty=1, col=gc.ramp[3],
  #      main="(b) Weights", ylim=c(0,1), xlab="t", ylab="w(t)",
  #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
  # points(tt[temp], ten[temp],pch=16, col=gc.ramp[3])
  # lines(tt, five, lwd=3, lty=2, col=gc.ramp[2])
  # lines(tt, fifteen, lwd=3, lty=3, col=gc.ramp[4])
  # lines(tt, twenty, lwd=3, lty=4, col=gc.ramp[5])
  # lines(tt, inf, lwd=3, lty=5, col=gc.ramp[6])
  # lines(tt, wgts0, lwd=3, lty=6, col=gc.ramp[7])
  # legend(0,0.7, legend=c('a=5', 'a=10', 'a=15', 'a=20', 'a=inf', '0-wgts'),
  #        col=c(gc.ramp[2],gc.ramp[3],gc.ramp[4],gc.ramp[5],gc.ramp[6],gc.ramp[7]),
  #        lty=c(1,2,3,4,5,6), lwd=3, cex=1.6)
  # grid()
  # dev.off()
  
  # 
  ## Weights as functional data -----------------------------------------
  ## Smooth the weights on a B-spline basis of degree 2
  #basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=4)
  
  basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=2)
  wgts.fd  <- smooth.basis(tt, wgts.obs, basis)$fd
  
  ## Output -------------------------------------------------------------
  
  out        <- list(wgts.obs, wgts.fd)
  names(out) <- c('wgts.obs', 'wgts.fd')
  
  return(out)
  
}

