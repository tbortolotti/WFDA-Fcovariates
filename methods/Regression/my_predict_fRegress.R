my_predict_fRegress <- function(mod, xlist, t.points)
{
  
  ## Utilities -----------------------------------------------------------------
  y         <- mod$yfdobj
  p         <- length(xlist)
  n         <- dim(xlist[[1]]$coefs)[2]
  rangeval  <- range(t.points)
  ybasisobj <- y$basis
  ynbasis   <- ybasisobj$nbasis
  nfine     <- max(501,10*ynbasis+1)
  tfine     <- seq(rangeval[1], rangeval[2], len=nfine)
  
  ## Evaluate it component per component and function by function --------------
  yhatmat  <- matrix(0, nrow=nfine, ncol=n)
  for (j in 1:p) {
    xfdj       <- eval.fd(tfine, xlist[[j]], 0, FALSE)
    betafdParj <- mod$betaestlist[[j]]
    betafdj    <- betafdParj$fd
    betavecj   <- eval.fd(tfine, betafdj, 0, FALSE)
    yhatmat    <- yhatmat + as.vector(betavecj)*xfdj
  }
  yhatfdobj <- smooth.basis(tfine, yhatmat, ybasisobj)$fd
  
  ## Output --------------------------------------------------------------------
  return(yhatfdobj)
  
}