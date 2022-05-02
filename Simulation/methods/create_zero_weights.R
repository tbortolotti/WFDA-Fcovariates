create_zero_weights <- function(curves.rec, t.points, breaks, reconst_fcts, Thp, log.flag=FALSE)
{
  ## Utilities -----------------------------------------------------------------
  n <- dim(curves.rec)[2]
  
  ## Building the weights ------------------------------------------------------
  if(log.flag)
  {
    tt <- t.points
  } else {
    step <- round(min(t.points[2:length(t.points)] - t.points[1:(length(t.points)-1)]), digits=2)
    tt <- seq(range(t.points)[1], range(t.points)[2], by=step)
  }
  
  
  wgts.obs <- matrix(data=1, nrow=length(tt), ncol=n)
  
  for(i in 1:length(reconst_fcts))
  {
    t.max    <- Thp[reconst_fcts[i]]
    imax     <- tail(which(t.points<=t.max),n=1)
    imax.new <- tail(which(tt<=t.max),n=1)
    
    x        <- tt[(imax.new+1):length(tt)]
    
    y.x      <- rep(1e-6, length(x))
    y        <- rep(1,imax.new)
    y        <- c(y,y.x)
    
    wgts.obs[,reconst_fcts[i]] <- y
  }
  
  ## Weights as functional data -----------------------------------------
  ## Smooth the weights on a B-spline basis of degree 2
  basis   <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=1)
  wgts.fd <- smooth.basis(tt, wgts.obs, basis)$fd
  
  ## Output -------------------------------------------------------------
  
  out        <- list(wgts.obs, wgts.fd)
  names(out) <- c('wgts.obs', 'wgts.fd')
  
  return(out)
  
}