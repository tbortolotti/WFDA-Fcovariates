my_gcvKneipLiebl <- function(fpca_obj, argvalsO, method, pev = 0.99, progrbar = FALSE){
  
  ## Useful methods --------------------------------------------------
  source('methods/Reconstruction/KL/is_error.R')
  source('methods/Reconstruction/KL/my_quadWeights.R')
  source('methods/Reconstruction/KL/my_reconstKneipLiebl_fun.R')
  
  ##
  Y           <- fpca_obj$Y # data pre-processed by irreg2mat()
  mu          <- fpca_obj$mu
  cov         <- fpca_obj$cov
  sigma2      <- fpca_obj$sigma2
  argvals     <- fpca_obj$argvals
  efunctionsP <- fpca_obj$efunctionsP
  evaluesP    <- fpca_obj$evaluesP
  npcP        <- length(evaluesP)
  a           <- min(argvals)
  b           <- max(argvals)
  ##
  ## a function is considered as complete if: 
  ## 1. its first and its last observation are non-missing
  ## 2. not considered so far
  compl_fcts <- apply(Y, 1, function(x){
    !is.na(utils::head(x, 1)) & !is.na(utils::tail(x, 1)) 
    # & length(c(stats::na.omit(x)))>=floor(.8*length(argvals))
  })
  ##
  # functions to be reconstructed
  Y.compl       <- Y[compl_fcts,,drop=FALSE]
  # make artifical fragements (all over argvalsO)
  locO          <- match(argvalsO, argvals)
  locM          <- c(1:length(argvals))[-locO]
  Y.pred        <- Y.compl
  Y.pred[,locM] <- NA
  ## the artifical fragments must have sufficiently many observations (>=5) for the pre-smoothing
  slct          <- apply(Y.pred, 1, function(x){length(c(stats::na.omit(x)))>=5})
  Y.pred        <- Y.pred[slct,,drop=FALSE]
  n_compl       <- nrow(Y.pred)
  ##
  if(any(method == c(1,2,3,4,5))){ # all, but PACE and PACE_E0
    muO               <- mu[locO]
    w                 <- my_quadWeights(argvalsO, method = "trapezoidal")
    Wsqrt             <- diag(sqrt(w))
    Winvsqrt          <- diag(1/(sqrt(w)))
    # CovOO
    V                 <- Wsqrt %*% cov[locO,locO] %*% Wsqrt
    evalues           <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
    evalues           <- replace(evalues, which(evalues <= 0), 0)
    ##
    npc               <- length(evalues[evalues>0])
    npc               <- ifelse(is.null(pev), npc, which(cumsum(evalues[evalues>0])/sum(evalues[evalues>0])>=pev)[1])
    ##
    efunctionsO       <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
    evaluesO          <- evalues[1:npc]  
    D.inv             <- diag(1/evaluesO)
    Z                 <- efunctionsO
    ## ##################################################################
    ## Reconstructive eigenfunctions
    efun_reconst      <- matrix(NA, nrow=length(argvals), ncol=npc)
    ##
    for(k in seq_len(npc)){
      efun_reconst[,k] <- apply(X   = cov[locO,,drop=FALSE], MARGIN = 2,
                                FUN = function(x){pracma::trapz(x=argvalsO, efunctionsO[,k] * x)})
      efun_reconst[,k] <- efun_reconst[,k] / evaluesO[k]
    }
    ## ##################################################################
  } else {
    npc <- npcP  
  }
  ##
  if(n_compl <= 1){
    # warning("Too few complete functions; we use the pve=0.9 criterium")
    if(any(method == c(6,7))){
      K.pve <- which(cumsum(evaluesP[evaluesP>0])/sum(evaluesP[evaluesP>0])>=.9)[1]
    } else {
      K.pve <- which(cumsum(evalues[evalues>0])/sum(evalues[evalues>0])>=.9)[1]
    }
    return(K.pve)
  }
  ##
  rss_mat           <- matrix(NA, n_compl, npc)
  if(progrbar){
    ## Progress Bar
    cat("Select K via GCV:\n")
    pb                <- utils::txtProgressBar(min = 0, max = n_compl*npc, style = 3)
    counter           <- 0
  }
  ##
  for(i in seq_len(n_compl)){# i <- 3
    ## Needed for all following ###########################################
    Y.cent            <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, ncol(Y)))
    obs_locO          <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvalsO))
    ## ####################################################################
    ## 
    if(any(method==c(1,2,3))){
      ## Classical scores (via intergration)
      scoresO      <- apply(X      = efunctionsO[obs_locO,,drop=FALSE], 
                            MARGIN = 2, 
                            FUN    = function(ef){pracma::trapz(y=ef*c(stats::na.omit(Y.cent)),x=argvalsO[obs_locO])})
    }
    ##
    if(any(method==c(2,4))){
      ## Pre-smoothing of fragmO
      smooth.fit        <- suppressMessages(stats::smooth.spline(y=c(stats::na.omit((Y.pred[i,]))), x=argvalsO[obs_locO]))
      fragmO_presmooth  <- stats::predict(smooth.fit, argvalsO)$y
    }
    ##
    if(any(method==c(4,5))){
      ## CEscores (PACE scores with respect to efunctionsO)
      if(sigma2           ==  0){sigma2 <- 1e-6}
      if(length(obs_locO) < npc){
        npcO        <- length(obs_locO)
        Zcur        <- Z[obs_locO,1:npcO,drop=FALSE]
        ZtZ_sD.inv  <- solve(crossprod(Zcur) + sigma2 * D.inv[1:npcO,1:npcO])
        CE_scoresO  <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
        CE_scoresO  <- c(CE_scoresO, rep(0, npc-npcO))
      }else{
        Zcur        <- Z[obs_locO,,drop=FALSE]
        ZtZ_sD.inv  <- solve(crossprod(Zcur) + sigma2 * D.inv)
        CE_scoresO  <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
      }
    }
    if(any(method == c(6,7))){
      ## PACE (PACE scores with respect to efunctionsP)
      obs_locP        <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvals))
      Y.centP         <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, ncol(Y)))
      D.invP          <- diag(1/evaluesP, nrow = npc, ncol = npc)
      ##
      if(sigma2 == 0){sigma2 <- 1e-6}
      if(method == 7){sigma2 <- 0}
      if(length(obs_locP) < npc){
        npcO            <- length(obs_locP)
        ZcurP           <- efunctionsP[obs_locP, 1:npcO, drop=FALSE]
        ##
        ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcO, 1:npcO]), silent = TRUE)
        # while(is_error(ZtZ_sD.invP)){# for preventing singularities
        #   sigma2          <- sigma2 + .Machine$double.eps
        #   ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcO, 1:npcO]), silent = TRUE)
        # }
        if(is_error(ZtZ_sD.invP)){
          ZtZ_sD.invP     <- matrix(NA, npcO, npcO)
        }
        scoresP         <- c(ZtZ_sD.invP %*% t(ZcurP) %*% c(stats::na.omit(Y.centP)))
        scoresP         <- c(scoresP, rep(0, npc-npcO))
      }else{
        ZcurP           <- efunctionsP[obs_locP, , drop=FALSE]
        ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npc, 1:npc]), silent = TRUE)
        # while(is_error(ZtZ_sD.invP)){# for preventing singularities
        #   sigma2          <- sigma2 + .Machine$double.eps
        #   ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npc, 1:npc]), silent = TRUE)
        # }
        if(is_error(ZtZ_sD.invP)){
          ZtZ_sD.invP     <- matrix(NA, npc, npc)
        }
        scoresP         <- c(ZtZ_sD.invP %*% t(ZcurP) %*% c(stats::na.omit(Y.centP)))
      }
    }
    # if(method == 8){
    #   if(sigma2 == 0){sigma2 <- 1e-6}
    #   ##
    #   efun_reconst_orth <- pracma::gramSchmidt(efun_reconst)$Q
    #   Zcur              <- efun_reconst_orth
    #   Zcur              <- Zcur[obs_locO,1:npc,drop=FALSE]
    #   ZtZ_sD.inv        <- solve(crossprod(Zcur) + sigma2 * D.inv[1:npc,1:npc])
    #   CE_scores_orth    <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
    # }
    ##
    for(k in seq_len(npc)){# k <- 5
      ##
      if(method == 1){# Classical scores, with alignment of reconstr and fully observed fragmO
        result_tmp <- my_reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=scoresO,efun_reconst=efun_reconst,fragmO=Y.pred[i,locO], K=k)
      }
      if(method == 2){# Classical scores, with alignment of reconstr and fragmO
        result_tmp <- my_reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=scoresO,efun_reconst=efun_reconst,fragmO=fragmO_presmooth, K=k)
      }
      if(method == 3){# Classical scores, without alignment of reconstr and pre-smoothed fragmO
        result_tmp <- my_reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=scoresO,efun_reconst=efun_reconst,fragmO=NULL,K=k)
      }
      if(method == 4){# CEscores, with alignment of reconstr and pre-smoothed fragmO
        result_tmp <- my_reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=CE_scoresO,efun_reconst=efun_reconst,fragmO=fragmO_presmooth,K=k)
      }
      if(method == 5){# CEscores, without alignment of reconstr and fragmO
        result_tmp <- my_reconstKneipLiebl_fun(mu=mu,argvals=argvals,argvalsO=argvalsO,scoresO=CE_scoresO,efun_reconst=efun_reconst,fragmO=NULL,  K=k)
      }
      if(any(method == c(6,7))){# PACE
        result_tmp <- my_PACE_fun(mu=mu, argvals=argvals, scoresP=scoresP, efunctionsP=efunctionsP, K=k)
      }
      # if(method == 8){# Orth
      #   result_tmp <- reconstKneipLiebl_orth_fun(mu=mu,argvals=argvals,CE_scores_orth=CE_scores_orth,efun_reconst_orth=efun_reconst_orth,K=k)
      # }
      ##
      rss_mat[i,k] <- sum((result_tmp[['y_reconst']][locM] - Y.compl[i,locM])^2, na.rm = TRUE)
      # cat("i=",i," k=",k,"\n")
      if(progrbar){
        counter <- counter+1
        utils::setTxtProgressBar(pb, counter)
      }
    }
  }
  if(progrbar){close(pb)}
  ##
  gcv_k_vec <- colSums(rss_mat)/((1-1:npc/n_compl)^2)
  K.gcv     <- which.min(gcv_k_vec)
  ##
  return(K.gcv)
}