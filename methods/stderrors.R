  

stderrors <- function(Sphi = NULL,
                      obs.inc,
                      SigmaE = NULL,
                      returnMatrix = FALSE,
                      fRegressList)
{
  
  # Utilities -------------------------------------------------
  
  yfdobj        <- fRegressList$yfdobj
  xfdlist       <- fRegressList$xlist
  betalist      <- fRegressList$betalist
  betaestlist   <- fRegressList$betaestlist
  Cmat          <- fRegressList$Cmat
  
  ycoef         <- yfdobj$coefs
  ycoefdim      <- dim(ycoef)
  TT            <- ycoefdim[1]
  N             <- ycoefdim[2]
  ybasisobj     <- yfdobj$basis
  rangeval      <- ybasisobj$rangeval
  ynbasis       <- ybasisobj$nbasis
  
  nfine         <- max(501,10*ynbasis+1)
  tfine         <- seq(rangeval[1], rangeval[2], len=nfine)
  
  ybasismat     <- eval.basis(tfine, ybasisobj, 0, returnMatrix)
  
  p             <- length(betalist)
  ncoef         <- dim(Cmat)[1]
  
  #  -----------------------------------------------------------------------
  #        Compute pointwise standard errors of regression coefficients
  #               if both y2cMap and SigmaE are supplied.
  #        y2cMap is supplied by the smoothing of the data that defined
  #        the dependent variable.
  #        SigmaE has to be computed from a previous analysis of the data.
  #  -----------------------------------------------------------------------
  
  if (!(is.null(Sphi) || is.null(SigmaE))) {
    
    ybasismat = eval.basis(tfine, ybasisobj, 0, returnMatrix)
    
    deltat    = tfine[2] - tfine[1]
    
    #  compute BASISPRODMAT
    
    basisprodmat = matrix(0,ncoef,ynbasis*N) # pL x Ln --> questa ? la c2bMap
    
    mj2 = 0
    for (j in 1:p) {
      betafdParj = betalist[[j]]
      betabasisj = betafdParj$fd$basis
      ncoefj     = betabasisj$nbasis
      bbasismatj = eval.basis(tfine, betabasisj, 0, returnMatrix)
      xfdj       = xfdlist[[j]]
      onemat     = as.matrix(rep(1,length(tfine)))
      tempj      = onemat %*% xfdj
      #  row indices of BASISPRODMAT to fill
      mj1    = mj2 + 1
      mj2    = mj2 + ncoefj
      indexj = mj1:mj2
      #  inner products of beta basis and response basis
      #    weighted by covariate basis functions
      mk2 = 0
      for (k in 1:ynbasis) {
        #  col indices of BASISPRODMAT to fill
        mk1    = mk2 + 1
        mk2    = mk2 + N
        indexk = mk1:mk2
        tempk  = bbasismatj*ybasismat[,k]
        basisprodmat[indexj,indexk] =
          deltat*crossprod(tempk,tempj)
      }
    }
    
    #  compute variances of regression coefficient function values
    
    Sb = solve(Cmat,basisprodmat) #c2bMap
    
    A  = Sb %*% Sphi
    # Compute A x (SigmaE \otimes In)
    In          <- Diagonal(N, x=1)
    SigmaE.I    <- base::kronecker(SigmaE, In)
    bvar        <- A %*% SigmaE.I %*% t(A)

    betastderrlist = vector("list", p)
    mj2 = 0
    for (j in 1:p) {
      betafdParj = betalist[[j]]
      betabasisj = betafdParj$fd$basis
      ncoefj     = betabasisj$nbasis
      mj1 	     = mj2 + 1
      mj2 	     = mj2 + ncoefj
      indexj     = mj1:mj2
      bbasismat  = eval.basis(tfine, betabasisj, 0, returnMatrix)
      bvarj      = bvar[indexj,indexj]
      bstderrj   = sqrt(diag(as.matrix(bbasismat %*% bvarj %*% t(bbasismat))))
      bstderrfdj = smooth.basis(tfine, bstderrj, betabasisj)$fd
      betastderrlist[[j]] = bstderrfdj
    }
  } else {
    betastderrlist = NULL
    bvar           = NULL
    Sb             = NULL
  }
  
  #  -------------------------------------------------------------------
  #                       Set up output list object
  #  -------------------------------------------------------------------
  
  outputList <- list(yfdobj         = yfdobj,
                     xfdlist        = xfdlist,
                     betalist       = betalist,
                     betaestlist    = betaestlist,
                     Cmat           = Cmat,
                     SigmaE         = SigmaE,
                     betastderrlist = betastderrlist,
                     bvar           = bvar,
                     c2bMap         = Sb)
  
  
  return(outputList)
}


  