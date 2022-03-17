create_weights_simulation <- function(curves.rec, t.points, breaks, fix.par=2, reconst_fcts, Thp)
{
  ## Utilities -----------------------------------------------------------------
  n <- dim(curves.rec)[2]
  
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
    t.max    <- Thp[reconst_fcts[i]]
    imax     <- tail(which(t.points<=t.max),n=1)
    imax.new <- tail(which(tt<=t.max),n=1)
    T.max    <- t.points[imax]
    
    scale    <- fix.par*sd(curves.rec[imax,-reconst_fcts])
    loc      <- T.max
    
    x        <- tt[(imax.new+1):length(tt)]
    
    y.x      <- 2/(1+exp((x - loc)*scale))
    y        <- rep(1,imax.new)
    y        <- c(y,y.x)
    
    wgts.obs[,reconst_fcts[i]] <- y
  }
  
  ## Weights as functional data -----------------------------------------
  ## Smooth the weights on a B-spline basis of degree 2
  basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=2)
  wgts.fd  <- smooth.basis(tt, wgts.obs, basis)$fd
  wgts.ev <- eval.fd(t.points, wgts.fd)
  
  ## Output -------------------------------------------------------------
  
  out        <- list(wgts.ev, wgts.fd)
  names(out) <- c('wgts.obs', 'wgts.fd')
  
  return(out)
  
}

