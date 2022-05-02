create_weights_new <- function(curves.rec, t.points, breaks, fix.par=2, reconst_fcts, Thp)
{
  ## Utilities -----------------------------------------------------------------
  n <- dim(curves.rec)[2]
  
  ## Building the weights ------------------------------------------------------
  tt <- t.points

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
  wgts.fd <- smooth.basis(t.points, wgts.obs, basis)$fd
  
  ## Output -------------------------------------------------------------
  
  out        <- list(wgts.obs, wgts.fd)
  names(out) <- c('wgts.obs', 'wgts.fd')
  
  return(out)
  
}