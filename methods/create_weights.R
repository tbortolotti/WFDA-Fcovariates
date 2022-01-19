


create_weights <- function(curves, t.points, loc.par=0.3, reconst_fcts, T_hp)
{
  ## Utilities -----------------------------------------------------------------
  n <- dim(curves)[2]
  
  ## Building the weights ------------------------------------------------------
  ## For each observation, evaluate mu so as to be loc.par after
  ## the last one valid, and the scale so as to be T_max/fix.par, i.e. the
  ## last valid period divided by the standard deviation of the complete
  ## observations in that period.
  step <- round(min(t.points[2:length(t.points)] - t.points[1:(length(t.points)-1)]), digits=2)
  tt <- seq(range(t.points)[1], range(t.points)[2], by=step)

  wgts.obs <- matrix(data=1, nrow=length(tt), ncol=n)
  
  for(i in 1:length(reconst_fcts))
  {
    t.max    <- T_hp[reconst_fcts[i]]
    imax     <- tail(which(t.points<=t.max),n=1)
    imax.new <- tail(which(tt<=t.max),n=1)
    T.max    <- t.points[imax]
    
    fix.par  <- sd(curves[imax,-reconst_fcts])
    
    loc      <- T.max + loc.par
    scale    <- 1/fix.par
    
    y        <- rep(1,imax.new)
    x        <- tt[(imax.new+1):length(tt)]
    y.x      <- 1/(1+exp((x - loc)/scale)) + (1 - 1/(1+exp((x[1]-loc)/scale)))
    y        <- c(y,y.x)
    
    wgts.obs[,reconst_fcts[i]] <- y
  }
  
  ## Weights as functional data -----------------------------------------
  ## Smooth the weights on a B-spline basis of order 2
  
  breaks <- t.points
  
  basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=3)
  
  # lambda.vec <- seq(0.01, 0.1, by=0.01)
  # gcv.vec    <- numeric(length(lambda.vec))
  # for(j in 1:length(lambda.vec))
  # {
  #   fPar <- fdPar(fdobj=basis, Lfdobj=0, lambda=lambda.vec[j])
  #   gcv.vec[j] <- sum(smooth.basis(tt, wgts.obs, fPar)$gcv)/n
  # }
  # lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
  
  fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
  wgts.s  <- smooth.basis(tt, wgts.obs, fPar)
  wgts.fd <- wgts.s$fd
  
  wgts.ev <- eval.fd(t.points, wgts.fd)
  
  ## Output -------------------------------------------------------------
  
  out        <- list(wgts.ev, wgts.fd)
  names(out) <- c('wgts.obs', 'wgts.fd')
  
  return(out)
  
}

