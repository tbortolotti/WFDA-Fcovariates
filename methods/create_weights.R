


create_weights <- function(curves.rec, t.points, breaks, loc.par=1, reconst_fcts, Thp, set.log=FALSE)
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
    
    fix.par  <- sd(curves.rec[imax,-reconst_fcts])
    scale    <- 1/fix.par
    loc      <- T.max + loc.par
    
    if(set.log)
    {
      loc    <- 10^T.max + loc.par
      x      <- 10^tt[(imax.new+1):length(tt)]
    } else {
      loc    <- T.max + loc.par
      x      <- tt[(imax.new+1):length(tt)]
    }
    
    y.x      <- 1/(1+exp((x - loc)/scale)) + (1 - 1/(1+exp((x[1]-loc)/scale)))
    y        <- rep(1,imax.new)
    y        <- c(y,y.x)
    
    wgts.obs[,reconst_fcts[i]] <- y
  }
  
  ## Weights as functional data -----------------------------------------
  ## Smooth the weights on a B-spline basis of degree 2
  #basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=4)
  if(set.log)
  {
    basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=3)
    # exp.vec <- seq(-7,2,by=1)
    # lambda.vec <- 10^exp.vec
    # gcv.vec    <- numeric(length(lambda.vec))
    # for(j in 1:length(lambda.vec))
    # {
    #   fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=lambda.vec[j])
    #   gcv.vec[j] <- sum(smooth.basis(tt, wgts.obs, fPar)$gcv)/n
    # }
    # lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
    fPar <- fdPar(fdobj=basis, Lfdobj=1, lambda=0.1)
    wgts.fd  <- smooth.basis(tt, wgts.obs, fPar)$fd
    wgts.ev <- eval.fd(t.points, wgts.fd)
  } else {
    basis  <- create.bspline.basis(rangeval=range(breaks), breaks=breaks, norder=2)
    wgts.fd  <- smooth.basis(tt, wgts.obs, basis)$fd
    wgts.ev <- eval.fd(t.points, wgts.fd)
  }

  ## Output -------------------------------------------------------------
  
  out        <- list(wgts.ev, wgts.fd)
  names(out) <- c('wgts.obs', 'wgts.fd')
  
  return(out)
  
}

