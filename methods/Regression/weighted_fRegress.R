weighted_fRegress <- function(y, xfdlist, betalist, wgts=NULL, returnMatrix=FALSE) {
  
  #  FREGRESS  Fits a functional linear model using multiple
  #  functional independent variables with the dependency being
  #  pointwise or concurrent.
  #  The case of a scalar independent variable is included by treating
  #  it as a functional independent variable with a constant basis
  #  and a unit coefficient.
  #  
  #  Arguments:
  #  Y            ... an object for the dependent variable,
  #                   which may be:
  #                       a functional data object or a numerical vector
  #  XFDLIST      ... a list object of length p with each list
  #                   containing an object for an independent variable.
  #                   the object may be:
  #                       a functional data object or
  #                       a vector
  #                   if XFDLIST is a functional data object or a vector,
  #                   it is converted to a list of length 1.
  #  BETALIST     ... a list object of length p with each list
  #                   containing a functional parameter object for
  #                   the corresponding regression function.  If any of
  #                   these objects is a functional data object, it is
  #                   converted to the default functional parameter object.
  #                   if BETALIST is a functional parameter object
  #                   it is converted to a list of length 1.
  #  WGTS         ... a functional object containing the functional weight for each observation.
  #                   If NULL, the function performs an unweighted function-on-function regression.
  #  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
  #                   from a call to function BsplineS.  See this function for
  #                   enabling this option.
  #
  #  Returns FREGRESSLIST  ... A list containing seven members with names:
  #    yfdobj      ... first  argument of FREGRESS
  #    xfdlist     ... second argument of FREGRESS
  #    betalist    ... third  argument of FREGRESS
  #    betaestlist ... estimated regression functions
  #    yhatfdobj   ... functional data object containing fitted functions
  #    Cmat        ... coefficient matrix for the linear system defining
  #                    the regression coefficient basis vector
  #    Dmat        ... right side vector for the linear system defining
  #                    the regression coefficient basis vector
  #    Cmatinv     ... inverse of the coefficient matrix, needed for
  #                    function FREGRESS.STDERR that computes standard errors
  #    df          ... degrees of freedom for fit
  
  if (is.fdPar(y)) y <- y$fd

  arglist <- my_fRegressArgCheck(y, xfdlist, betalist)
  
  yfdobj   <- arglist$yfd
  xfdlist  <- arglist$xfdlist
  betalist <- arglist$betalist
  
  p <- length(xfdlist)

  #  ----------------------------------------------------------------
  #           YFDOBJ is a functional data object
  #  ----------------------------------------------------------------
  
  #  extract dependent variable information
  ycoef     <- yfdobj$coefs
  ycoefdim  <- dim(ycoef)
  N         <- ycoefdim[2]
  ybasisobj <- yfdobj$basis
  rangeval  <- ybasisobj$rangeval
  ynbasis   <- ybasisobj$nbasis
  onesbasis <- create.constant.basis(rangeval)
  onesfd    <- fd(1,onesbasis)
  
  if (length(ycoefdim) > 2) stop("YFDOBJ from YFD is not univariate.")
  
  #  --------  set up the linear equations for the solution  -----------
  
  #  compute the total number of coefficients to be estimated
  
  ncoef <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (betafdParj$estimate) {
      ncoefj     <- betafdParj$fd$basis$nbasis
      ncoef      <- ncoef + ncoefj
    }
  }
  
  Cmat <- matrix(0,ncoef,ncoef)
  Dmat <- rep(0,ncoef)
  
  #  ------------------------------------------------------------------------
  #  Compute the symmetric positive definite matrix CMAT and
  #  the column matrix DMAT.  CMAT contains weighted inner products of
  #  bases for each pair of terms plus, for lambda > 0, a roughness penalty
  #  matrix to ensure that the estimated coefficient functions will be smooth
  #  The weight vector is the point-wise product of the associated functional
  #  covariates.  
  #  Dmat contains for each covariate the weighted integral of the basis 
  #  functions, with the weight function being the covariate function
  #  pointwise multiplied the dependent variate yobj.
  #  The estimated coefficients functions are defined by the solution
  #  CMAT %*% COEF = DMAT.
  #  ------------------------------------------------------------------------
  
  #  loop through rows of CMAT
  
  mj2 <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (betafdParj$estimate) {
      #  get jth beta basis
      betafdj    <- betafdParj$fd
      betabasisj <- betafdj$basis
      ncoefj     <- betabasisj$nbasis
      #  row indices of CMAT and DMAT to fill
      mj1    <- mj2 + 1
      mj2    <- mj2 + ncoefj
      indexj <- mj1:mj2
      #  compute right side of equation DMAT
      #  compute weight function for DMAT
      xfdj <- xfdlist[[j]]
      if (is.null(wgts)) {
        xyfdj <- xfdj*yfdobj
      } else {
        xyfdj <- (xfdj*wgts)*yfdobj
      }
      wtfdj <- sum(xyfdj)
      #  Compute jth component of DMAT
      Dmatj <- inprod(betabasisj,onesfd,0,0,rangeval,wtfdj)
      Dmat[indexj] <- Dmatj
      #  loop through columns of CMAT
      mk2 <- 0
      for (k in 1:j) {
        betafdPark <- betalist[[k]]
        if (betafdPark$estimate) {
          #  get the kth basis
          betafdk    <- betafdPark$fd
          betabasisk <- betafdk$basis
          ncoefk     <- betabasisk$nbasis
          #  column indices of CMAT to fill
          mk1 <- mk2 + 1
          mk2 <- mk2 + ncoefk
          indexk <- mk1:mk2
          #  set up weight function for CMAT component
          xfdk <- xfdlist[[k]]
          if (is.null(wgts)) {
            xxfdjk <- xfdj*xfdk
          } else {
            xxfdjk <- (xfdj*wgts)*xfdk
          }
          wtfdjk <- sum(xxfdjk)
          #  compute the inner product
          Cmatjk <- inprod(betabasisj, betabasisk, 0, 0,
                           rangeval, wtfdjk)
          Cmat[indexj,indexk] <- Cmatjk
          Cmat[indexk,indexj] <- t(Cmatjk)
        }
      }
      #  attach penalty term to diagonal block if required
      lambdaj <- betafdParj$lambda
      if (lambdaj > 0) {
        Rmatj <- betafdParj$penmat
        if (is.null(Rmatj)) {
          Lfdj  <- betafdParj$Lfd
          Rmatj <- eval.penalty(betabasisj, Lfdj)
        }
        Cmat[indexj,indexj] <- Cmat[indexj,indexj] +
          lambdaj*Rmatj
      }
    }
  }
  
  #  ensure symmetry
  
  Cmat <- (Cmat+t(Cmat))/2
  
  Cmatinv <- solve(Cmat)
  
  betacoef <- Cmatinv %*% Dmat
  
  #  set up fdPar objects for reg. fns. in BETAESTLIST
  
  betaestlist <- betalist
  mj2 <- 0
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (betafdParj$estimate) {
      betafdj <- betafdParj$fd
      ncoefj  <- betafdj$basis$nbasis
      mj1     <- mj2 + 1
      mj2     <- mj2 + ncoefj
      indexj  <- mj1:mj2
      coefj   <- betacoef[indexj]
      betafdj$coefs <- as.matrix(coefj)
      betafdParj$fd <- betafdj
    }
    betaestlist[[j]] <- betafdParj
  }
  
  #  set up fd objects for predicted values in YHATFDOBJ
  
  nfine   <- max(501,10*ynbasis+1)
  tfine   <- seq(rangeval[1], rangeval[2], len=nfine)
  yhatmat <- matrix(0,nfine,N)
  for (j in 1:p) {
    xfdj       <- xfdlist[[j]]
    xmatj      <- eval.fd(tfine, xfdj, 0, returnMatrix)
    betafdParj <- betaestlist[[j]]
    betafdj    <- betafdParj$fd
    betavecj   <- eval.fd(tfine, betafdj, 0, returnMatrix)
    yhatmat    <- yhatmat + xmatj*as.vector(betavecj)
  }
  yhatfdobj <- smooth.basis(tfine, yhatmat, ybasisobj)$fd
  
  df <- NA
  
  #  -------------------------------------------------------------------
  #                       Set up output list object
  #  -------------------------------------------------------------------
  
  fRegressList <- list(yfdobj         = yfdobj,
                       xfdlist        = xfdlist,
                       betalist       = betalist,
                       betaestlist    = betaestlist,
                       yhatfdobj      = yhatfdobj,
                       Cmat           = Cmat,
                       Dmat           = Dmat,
                       Cmatinv        = Cmatinv,
                       df             = df)
  
  return(fRegressList)
  
}

## Auxiliary functions ---------------------------------------------------------

eigchk <- function(Cmat) {
  
  #  Last modified 25 August 2020 by Jim Ramsay
  
  #  Cmat for NA's
  
  if (any(is.na(Cmat))) stop("Cmat has NA values.")
  
  #  check Cmat for Cmatmetry
  
  if (max(abs(Cmat-t(Cmat)))/max(abs(Cmat)) > 1e-10) {
    stop('CMAT is not symmetric.')
  } else {
    Cmat <- (Cmat + t(Cmat))/2
  }
  
  #  check Cmat for singularity
  
  eigval <- eigen(Cmat)$values
  ncoef  <- length(eigval)
  if (eigval[ncoef] < 0) {
    neig <- min(length(eigval),10)
    cat("\nSmallest eigenvalues:\n")
    print(eigval[(ncoef-neig+1):ncoef])
    cat("\nLargest  eigenvalues:\n")
    print(eigval[1:neig])
    stop("Negative eigenvalue of coefficient matrix.")
  }
  if (eigval[ncoef] == 0) stop("Zero eigenvalue of coefficient matrix.")
  logcondition <- log10(eigval[1]) - log10(eigval[ncoef])
  if (logcondition > 12) {
    warning("Near singularity in coefficient matrix.")
    cat(paste("\nLog10 Eigenvalues range from\n",
              log10(eigval[ncoef])," to ",log10(eigval[1]),"\n"))
  }
}



my_fRegressArgCheck <- function(yfd, xfdlist, betalist) 
{
  #  FREGRESS_ARGCHECK checks the first three arguments for the functions
  #  for function regression, including FREGRESS.
  
  #  --------------------  Check classes of arguments  --------------------
  
  #  check that YFD is of class either 'fd' or 'numeric' and compute sample size N
  
  if (!(is.fdPar(yfd) || is.fd(yfd) || is.numeric(yfd) || is.matrix(yfd))) stop(
    "First argument is not of class 'fdPar', 'fd', 'numeric' or 'matrix'.")
  
  if (is.fdPar(yfd)) yfd <- yfd$fd
  
  if (inherits(yfd, "fd")) {
    ycoef <- yfd$coefs
    N     <- dim(ycoef)[2]
  } else {
    N <- length(yfd)
  } 
  
  #  check that xfdlist is a list object and compute number of covariates p
  
  #  check XFDLIST
  
  if (inherits(xfdlist, "fd") || inherits(xfdlist, "numeric")) 
    xfdlist <- list(xfdlist)
  
  if (!inherits(xfdlist, "list")) stop(
    "Argument XFDLIST is not a list object.")
  
  #  get number of independent variables p
  
  p <- length(xfdlist)
  
  #  check BETALIST
  
  if (inherits(betalist, "fd")) betalist <- list(betalist)
  
  if (!inherits(betalist, "list")) stop(
    "Argument BETALIST is not a list object.")
  
  if (length(betalist) != p)  {
    cat(paste("\nNumber of regression coefficients does not match\n",
              "number of independent variables."))
    stop("")
  }
  
  #  extract the range if YFD is functional
  
  if (inherits(yfd, "fd")) {
    rangeval <- yfd$basis$rangeval
  } else {
    rangeval = c(0,1)
  }
  
  #  --------------------  check contents of XFDLIST  -------------------
  
  #  If the object is a vector of length N,
  #  it is converted to a functional data object with a
  #  constant basis
  
  onebasis <- create.constant.basis(rangeval)
  onesfd   <- fd(1,onebasis)
  
  xerror <- FALSE
  for (j in 1:p) {
    xfdj <- xfdlist[[j]]
    if (inherits(xfdj, "fd")) {
      xcoef <- xfdj$coefs
      if (length(dim(xcoef)) > 2) stop(
        paste("Covariate",j,"is not univariate."))
      #  check size of coefficient array
      Nj <- dim(xcoef)[2]
      if (Nj != N) {
        print(
          paste("Incorrect number of replications in XFDLIST",
                "for covariate",j))
        xerror = TRUE
      }
    } 
    if (inherits(xfdj, "numeric")) {
      if (!is.matrix(xfdj)) xfdj = as.matrix(xfdj)
      Zdimj <- dim(xfdj)
      if (Zdimj[1] != N && Zdimj != 1) {
        print(paste("Vector in XFDLIST[[",j,"]] has wrong length."))
        xerror = TRUE 
      } 
      if (Zdimj[2] != 1) {
        print(paste("Matrix in XFDLIST[[",j,"]] has more than one column."))
        xerror = TRUE 
      } 
      xfdlist[[j]] <- fd(matrix(xfdj,1,N), onebasis)
    } 
    if (!(inherits(xfdlist[[j]], "fd"     ) || 
          inherits(xfdlist[[j]], "numeric") ||
          inherits(xfdlist[[j]], "matrix" ))) {
      print(paste("XFDLIST[[",j,"]] is not an FD or numeric or matrix object."))
      xerror = TRUE
    }
  }
  
  #  --------------------  check contents of BETALIST  -------------------
  
  berror <- FALSE
  for (j in 1:p) {
    betafdParj <- betalist[[j]]
    if (inherits(betafdParj, "fd") || inherits(betafdParj, "basisfd")) {
      betafdParj    <- fdPar(betafdParj)
      betalist[[j]] <- betafdParj
    }
    if (!inherits(betafdParj, "fdPar")) {
      print(paste("BETALIST[[",j,"]] is not a FDPAR object."))
      berror <- TRUE
    }
  }
  
  if (xerror || berror) stop(
    "An error has been found in either XFDLIST or BETALIST.")
  
  #  ---------------------  return the argument list  --------------------
  
  # The older versions of fda package used yfdPar as the name for the first member.
  
  return(list(yfd=yfd, xfdlist=xfdlist, betalist=betalist))
  
}
