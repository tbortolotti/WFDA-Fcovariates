my_fpca <- function(Ly, Lu, reconst_fcts = NULL, pev = 0.99, CEscores = FALSE, PACE = FALSE, PACE_E0 = FALSE, center = TRUE, maxbins = NULL){
  
  ## Useful method ------------------------------------------------
  source('KL/my_irreg2mat.R')
  source('KL/my_quadWeights.R')
  
  ##
  n      <- length(Ly)
  id_vec <- NULL
  for(i in 1:n){id_vec <- c(id_vec, rep(i, length(Ly[[i]])))}
  ##
  ydata  <-  data.frame(".id"    = id_vec, 
                        ".index" = unname(unlist(Lu)), 
                        ".value" = unname(unlist(Ly)))
  ##
  nbasis         <- 10
  maxbins        <- ifelse(is.null(maxbins), 1000, maxbins)
  useSymm        <- FALSE # if true, there will be no smoothing accross the diagonal
  makePD         <- FALSE # if true, the NP-estimate of cov is forced to be positive definite
  ##
  # Y (nobs x length(argvals)) with NA's if no data at a argvalue
  Y        <- my_irreg2mat(ydata, binning = TRUE, maxbins = maxbins) 
  argvals  <- as.numeric(colnames(Y))
  ##
  # functions to be reconstructed
  if(is.null(reconst_fcts)){reconst_fcts <- 1:n}
  Y.pred   <- Y[reconst_fcts,,drop=FALSE]
  # argvals of observed fragments to be reconstructed
  argvalsO <- vector("list", length = length(reconst_fcts))
  for(i in seq_len(length(reconst_fcts))){
    minmax        <- range(as.numeric(names(c(stats::na.omit(Y.pred[i,])))))
    argvalsO[[i]] <- argvals[argvals >= minmax[1] & argvals <= minmax[2]]
  }
  ##
  D      <- NCOL(Y)      # nobs per function
  I      <- NROW(Y)      # number of functions
  I.pred <- NROW(Y.pred) # number of functions to be reconstructed
  ##
  d.vec <- rep(argvals, each = I)
  id    <- rep(1:I, rep(D, I))
  
  ## MEAN ##############################################################
  if(center){
    ## mean
    gam0    <- mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu      <- mgcv::predict.gam(gam0, newdata = data.frame(d.vec = argvals))
    Y.tilde <- Y - matrix(mu, I, D, byrow = TRUE)
  }else{
    Y.tilde <- Y
    mu      <- rep(0, D)
  }
  ## plot(x=argvals,y=mu, type="l")
  ##
  
  ## COV ###############################################################
  ## 1. pointwise (at argvalues) sample covariance matrix (=naive cov-estimator)
  ## 2. smooth this matrix
  cov.sum = cov.count = cov.mean = matrix(0, D, D)
  for(i in 1:I){
    obs.points = which(!is.na(Y[i, ]))
    cov.count[obs.points, obs.points] <- cov.count[obs.points, obs.points] + 1
    cov.sum[  obs.points, obs.points] <- cov.sum[  obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
  }
  G.0       <- ifelse(cov.count == 0, NA, cov.sum/cov.count)
  diag.G0   <- diag(G.0)
  diag(G.0) <- NA
  if(!useSymm){
    row.vec <- rep(argvals, each = D)
    col.vec <- rep(argvals, D)
    npc.0   <- matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis),
                                                  weights = as.vector(cov.count)), 
                                        newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
    npc.0 = (npc.0 + t(npc.0))/2
  }else{
    use          <- upper.tri(G.0, diag = TRUE)
    use[2, 1]    <- use[ncol(G.0), ncol(G.0) - 1] <- TRUE
    usecov.count <- cov.count
    usecov.count[2, 1] <- usecov.count[ncol(G.0), ncol(G.0) - 1] <- 0
    usecov.count <- as.vector(usecov.count)[use]
    use          <- as.vector(use)
    vG.0         <- as.vector(G.0)[use]
    row.vec      <- rep(argvals, each = D)[use]
    col.vec      <- rep(argvals, times = D)[use]
    mCov         <- mgcv::gam(vG.0 ~ te(row.vec, col.vec, k = nbasis), weights = usecov.count)
    npc.0        <- matrix(NA, D, D)
    spred        <- rep(argvals, each = D)[upper.tri(npc.0, diag = TRUE)]
    tpred        <- rep(argvals, times = D)[upper.tri(npc.0, diag = TRUE)]
    # Estimated covariance function:
    smVCov       <- mgcv::predict.gam(mCov, newdata = data.frame(row.vec = spred, col.vec = tpred))
    npc.0[upper.tri(npc.0, diag = TRUE)] <- smVCov
    npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
    # slct <- seq.int(1,length(argvals),len=25)
    # persp(z=npc.0[slct,slct],x=argvals[slct],y=argvals[slct])
  }
  if(makePD){
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, trace = TRUE)
      as.matrix(tmp$mat)
    }
  }
  # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch. 8)
  w          <- my_quadWeights(argvals, method = "trapezoidal")
  Wsqrt      <- diag(sqrt(w))
  Winvsqrt   <- diag(1/(sqrt(w)))
  V          <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  evalues    <- replace(evalues, which(evalues <= 0), 0)
  npc        <- length(evalues[evalues>0])
  npc        <- ifelse(is.null(pev), npc, which(cumsum(evalues[evalues>0])/sum(evalues[evalues>0])>=pev)[1])
  efunctions <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
  evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
  # Estimated covariance function
  cov        <- efunctions %*% tcrossprod(diag(evalues), efunctions)
  # Numerical integration for estimation of sigma2
  T.len      <- argvals[D] - argvals[1]  # total interval length
  T1.min     <- min(which(argvals >= argvals[1] + 0.25 * T.len))  # left bound of narrower interval T1
  T1.max     <- max(which(argvals <= argvals[D] - 0.25 * T.len))  # right bound of narrower interval T1
  DIAG       <- (diag.G0 - diag(cov))[T1.min:T1.max]  # function values
  w2         <- my_quadWeights(argvals[T1.min:T1.max], method = "trapezoidal")
  sigma2     <- max(stats::weighted.mean(DIAG, w = w2, na.rm = TRUE), 0)
  ##
  
  ## PACE ################################################################################
  if(PACE | PACE_E0){
    scoresP           <- vector(mode = "list", length(reconst_fcts))
    # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
    wP                <- my_quadWeights(argvals, method = "trapezoidal")
    WsqrtP            <- diag(sqrt(wP))
    WinvsqrtP         <- diag(1/(sqrt(wP)))
    # Cov
    VP                <- WsqrtP %*% cov %*% WsqrtP
    evalP             <- eigen(VP, symmetric = TRUE, only.values = TRUE)$values
    evalP             <- replace(evalP, which(evalP <= 0), 0)
    npcP              <- length(evalP[evalP>0])  
    npcP              <- ifelse(is.null(pev), npcP, which(cumsum(evalP[evalP>0])/sum(evalP[evalP>0])>=pev)[1])
    efunctionsP       <- matrix(WinvsqrtP %*% eigen(VP, symmetric = TRUE)$vectors[, seq(len = npcP)], nrow = nrow(VP), ncol = npcP)
    evaluesP          <- evalP[1:npcP] 
    ##
    D.invP            <- diag(1/evaluesP, nrow = npcP, ncol = npcP)
    ##
    for(i in seq_len(length(reconst_fcts))){# i <- 1
      obs_locP        <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvals))
      Y.centP         <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, D))
      ##
      if(sigma2 == 0){sigma2 <- 1e-6}
      if(PACE_E0){sigma2  <- 0}
      if(length(obs_locP) < npcP){npcP <- length(obs_locP)}
      ZcurP           <- efunctionsP[obs_locP, 1:npcP, drop=FALSE]
      ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcP, 1:npcP]), silent = TRUE)
      # while(is.error(ZtZ_sD.invP)){# for preventing singularities
      #   sigma2          <- sigma2 + .Machine$double.eps
      #   ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcP, 1:npcP]), silent = TRUE)
      # }
      if(is.error(ZtZ_sD.invP)){
        ZtZ_sD.invP     <- matrix(NA, npcP, npcP)
      }
      scoresP[[i]]    <- c(ZtZ_sD.invP %*% t(ZcurP) %*% c(stats::na.omit(Y.centP)))
    }
  } else {
    efunctionsP <- NA
    evaluesP    <- NA
    scoresP     <- NA  
  }
  ## End PACE ############################################################################
  
  
  ## computations for observed fragments
  muO          <- vector("list", length(reconst_fcts))
  scoresO      <- vector("list", length(reconst_fcts))
  CE_scoresO   <- vector("list", length(reconst_fcts))
  evaluesO     <- vector("list", length(reconst_fcts))
  efunctionsO  <- vector("list", length(reconst_fcts))
  efun_reconst <- vector("list", length(reconst_fcts))
  ##
  obs_argvalsO <- vector("list", length(reconst_fcts))
  ##
  for(i in seq_len(length(reconst_fcts))){# i <- 1
    # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
    w                 <- my_quadWeights(argvalsO[[i]], method = "trapezoidal")
    Wsqrt             <- diag(sqrt(w))
    Winvsqrt          <- diag(1/(sqrt(w)))
    locO              <- match(argvalsO[[i]],argvals)
    # CovOO
    VO                <- Wsqrt %*% cov[locO,locO] %*% Wsqrt
    evalO             <- eigen(VO, symmetric = TRUE, only.values = TRUE)$values
    evalO             <- replace(evalO, which(evalO <= 0), 0)
    npcO              <- length(evalO[evalO>0])  
    npcO              <- ifelse(is.null(pev), npcO, which(cumsum(evalO[evalO>0])/sum(evalO[evalO>0])>=pev)[1])
    efunctionsO[[i]]  <- matrix(Winvsqrt %*% eigen(VO, symmetric = TRUE)$vectors[, seq(len = npcO)], nrow = nrow(VO), ncol = npcO)
    evaluesO[[i]]     <- evalO[1:npcO]  
    ##
    D.inv             <- diag(1/evaluesO[[i]])
    Z                 <- efunctionsO[[i]]
    Y.cent            <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, D))
    obs_locO          <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvalsO[[i]]))
    obs_argvalsO[[i]] <- argvalsO[[i]][obs_locO]
    ## 
    if(CEscores){
      ## CEScores (i.e., PACE-Scores)
      if(sigma2           ==   0){sigma2 <- 1e-6}
      if(length(obs_locO) < npcO){npcO <- length(obs_locO)}
      ##
      Zcur           <- Z[obs_locO,1:npcO,drop=FALSE]
      ZtZ_sD.inv     <- solve(crossprod(Zcur) + sigma2 * D.inv[1:npcO,1:npcO])
      CE_scoresO[[i]] <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
    } else {
      CE_scoresO[[i]] <- NA
    }
    ## Classical scores (via intergral approximation)
    scoresO[[i]] <- apply(X      = efunctionsO[[i]][obs_locO,,drop=FALSE], 
                          MARGIN = 2, 
                          FUN    = function(ef){pracma::trapz(y=ef*c(stats::na.omit(Y.cent)),x=obs_argvalsO[[i]])})
    ##
    muO[[i]]     <- mu[locO]
    ## ##################################################################
    ## Reconstructive eigenfunctions
    efun_reconst[[i]]  <- matrix(NA, nrow=length(argvals), ncol=npcO)
    ##
    for(k in seq_len(npcO)){
      efun_reconst[[i]][,k] <- apply(X   = cov[locO,,drop=FALSE], MARGIN = 2,
                                     FUN = function(x){pracma::trapz(x=argvalsO[[i]], efunctionsO[[i]][,k] * x)})
      efun_reconst[[i]][,k] <- efun_reconst[[i]][,k] / evaluesO[[i]][k]
    }
    ## ##################################################################
  }
  
  ## Return results ##################################################
  ret.objects <- c("Y", "mu", "muO", "cov", "sigma2",
                   "argvals",    "argvalsO", "obs_argvalsO",
                   "CE_scoresO",  "scoresO", "scoresP",
                   "efunctions", "efunctionsO", "efun_reconst", "efunctionsP", 
                   "evalues",    "evaluesO", "evaluesP")
  ret         <- lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret)  <- ret.objects
  return(ret)
}
