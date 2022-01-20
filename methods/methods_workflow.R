methods_workflow <- function(b,
                             T.period,
                             t.points,
                             T_hp,
                             log.Thp,
                             curves,
                             events,
                             B,
                             xlist,
                             blist,
                             method = c('interpolation', 'interpolation-noweight',
                                        'Kraus1', 'Kraus2', 'KLNoAl' ,'KLAl'),
                             loc=NULL)
{

  ## Methods
  source('methods/find_obs_inc.R')
  source('methods/interpolation.R')
  source('methods/create_weights.R')
  source('methods/wt_bsplinesmoothing.R')
  source('methods/Reconstruction/my_reconstructKraus1.R')
  source('methods/Reconstruction/finegrid_evaluation.R')
  source('methods/Reconstruction/my_reconstructKneipLiebl.R')
  
  source('methods/Regression/lambda_select.R')
  source('methods/Regression/weighted_fRegress.R')
  source('methods/Regression/my_predict_fRegress.R')
  source('methods/Regression/eval_MSE.R')
  
  ## Utilities for the identification of the b-th batch
  n      <- dim(curves)[2]
  
  #### Separate training and test set ------------------------------------------
  unique.events <- unique(events)
  n.events <- length(unique.events)
  n.test.events <- floor(n.events/B)
  
  set.seed(140996)
  shuffled.unique.events <- sample(unique.events, n.events)
  
  idxs <- ((b-1)*n.test.events + 1):min((b*n.test.events),n.events)
  test.events <- shuffled.unique.events[idxs]
  
  test <- (events %in% test.events)
  train <- (!test)
  
  curves.test    <- curves[,test]
  curves.train   <- curves[,train]
  
  # Create xlist.test and xlist.train
  xlist.test  <- list()
  xlist.train <- list()
  for(i in 1:length(xlist))
  {
    xlist.test[[i]]  <- xlist[[i]]
    xlist.test[[i]]$coefs  <- xlist[[i]]$coefs[,test]
    
    xlist.train[[i]] <- xlist[[i]]
    xlist.train[[i]]$coefs <- xlist[[i]]$coefs[,train]
  }
  
  reconst_fcts   <- find_obs_inc(Y = curves.train)
  
  T_hp.test      <- T_hp[test]
  T_hp.train     <- T_hp[train]
  
  log.Thp.test   <- log.Thp[test]
  log.Thp.train  <- log.Thp[train]
  
  #### Interpolation -----------------------------------------------------------
  if(any(method == 'interpolation'))
  {
    interpolate   <- interpolation(curves       = curves.train,
                                   t.points     = T.period,
                                   T_hp         = T_hp.train,
                                   reconst_fcts = reconst_fcts)
    curves.interp <- interpolate$curves.rec
    
    ## Construction of the weights
    wgt           <- create_weights(curves        = curves.interp,
                                    t.points      = t.points,
                                    loc.par       = loc,
                                    reconst_fcts  = reconst_fcts,
                                    T_hp          = log.Thp.train)
    
    ## Smoothing
    smth             <- wt_bsplinesmoothing(curves   = curves.interp,
                                            wgts.obs = wgt$wgts.obs,
                                            t.points = t.points,
                                            lambda   = 1e-5,
                                            set.cb   = FALSE)
    curves.interp.fd <- smth$curves.fd
    
    ### B-list
    # fPar  <- smth$fPar
    # blist.fixed <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)
    # 
    # for(q in 1:dim(XX)[2])
    # {
    #   blist.fixed[[q]]$lambda <- lambda.fixed[[q]]
    # }
    # blist <- blist.fixed
    
    ### Selection of the optimal lambda for the regressors that are
    ### affected the most by its value, in terms of resulting MSE.
    ### These regressors are:
    ### 1. intercept
    ### 2. low magnitudes
    ### 7. geometrical spreading
    
    ## Selection of the optimal lambda
    #blist <- lambda_select(curves.fd     = curves.interp.fd,
    #                       xlist         = xlist.train,
    #                       t.points      = T.period,
    #                       reg.idxs      = 7,
    #                       blist.fixed   = blist.fixed,
    #                       XX            = X.train,
    #                       wgts.fd       = wgt$wgts.fd,
    #                       interpolation = TRUE)
    
    ## Regression and beta estimation
    mod <- weighted_fRegress(y            = curves.interp.fd,
                             xfdlist      = xlist.train,
                             betalist     = blist,
                             wgts         = wgt$wgts.fd)
    
    ## y_hat and MSE
    curves.interp.hat   <- my_predict_fRegress(mod          = mod,
                                               xlist        = xlist.test,
                                               t.points     = t.points)
    curves.interp.hat.v <- eval.fd(t.points, curves.interp.hat)
    
    evalMSE    <- eval_MSE(curves     = curves.test,
                           curves.hat = curves.interp.hat.v,
                           t.points   = t.points)
    MSE <- evalMSE$MSE
    MSE_reconstruction <- evalMSE$MSE_reconstruction
    
  }
  
  #### Interpolation and NO WEIGHT ---------------------------------------------
  if(any(method == 'interpolation-noweight'))
  {
    interpolate   <- interpolation(curves       = curves.train,
                                   t.points     = T.period,
                                   T_hp         = T_hp.train,
                                   reconst_fcts = reconst_fcts)
    curves.interp <- interpolate$curves.rec

    ## Smoothing
    basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
    fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-5)
    curves.interp.fd <- smooth.basis(t.points, curves.interp, fPar)$fd
    
    ### B-list
    # fPar  <- smth$fPar
    # blist.fixed <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)
    # 
    # for(q in 1:dim(XX)[2])
    # {
    #   blist.fixed[[q]]$lambda <- lambda.fixed[[q]]
    # }
    # blist <- blist.fixed
    
    ### Selection of the optimal lambda for the regressors that are
    ### affected the most by its value, in terms of resulting MSE.
    ### These regressors are:
    ### 1. intercept
    ### 2. low magnitudes
    ### 7. geometrical spreading
    
    ## Selection of the optimal lambda
    #blist <- lambda_select(curves.fd     = curves.interp.fd,
    #                       xlist         = xlist.train,
    #                       t.points      = T.period,
    #                       reg.idxs      = 7,
    #                       blist.fixed   = blist.fixed,
    #                       XX            = X.train,
    #                       wgts.fd       = wgt$wgts.fd,
    #                       interpolation = TRUE)
    
    ## Regression and beta estimation
    mod <- fRegress(y            = curves.interp.fd,
                    xfdlist      = xlist.train,
                    betalist     = blist,
                    returnMatrix = FALSE,
                    method       = 'fRegress',
                    sep          = '.')
    
    ## y_hat and MSE
    curves.interp.hat   <- my_predict_fRegress(mod          = mod,
                                               xlist        = xlist.test,
                                               t.points     = t.points)
    curves.interp.hat.v <- eval.fd(t.points, curves.interp.hat)
    
    evalMSE    <- eval_MSE(curves     = curves.test,
                           curves.hat = curves.interp.hat.v,
                           t.points   = t.points)
    MSE <- evalMSE$MSE
    MSE_reconstruction <- evalMSE$MSE_reconstruction
    
  }
  
  #### Kraus 1 -----------------------------------------------------------------
  if(any(method=='Kraus1'))
  {
    log.tfine         <- seq(range(t.points)[1], range(t.points)[2], by=0.01)
    tfine             <- 10^log.tfine
    tfine[1]          <- 0
    curves.train.fine <- finegrid_evaluation(curves   = curves.train,
                                             t.points = t.points,
                                             lambda   = 1e-5,
                                             tfine    = log.tfine)
    
    Kraus1            <- my_reconstructKraus1(X_mat        = curves.train.fine,
                                              tfine        = log.tfine,
                                              alpha        = 0.002,
                                              reconst_fcts = reconst_fcts)
    curves.Kraus1     <- Kraus1[['X_compl_mat']]
    
    ## Smoothing
    basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
    fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-5)
    curves.Kraus1.fd <- smooth.basis(t.points, curves.Kraus1, fPar)$fd
    
    # ## B-list
    # fPar  <- smth$fPar
    # blist.fixed <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)
    # 
    # for(q in 1:dim(XX)[2])
    # {
    #   blist.fixed[[q]]$lambda <- lambda.fixed[[q]]
    # }
    # blist <- blist.fixed
    # 
    ## Selection of the optimal lambda
    #blist <- lambda_select(curves.fd   = curves.Kraus1.fd,
    #                       xlist       = xlist.train,
    #                       t.points    = T.period,
    #                       reg.idxs    = c(2,7,9),
    #                       blist.fixed = blist.fixed,
    #                       XX          = X.train)
    
    ## Regression and beta estimation
    mod <- fRegress(y            = curves.Kraus1.fd,
                    xfdlist      = xlist.train,
                    betalist     = blist,
                    returnMatrix = FALSE,
                    method       = 'fRegress',
                    sep          = '.')
    
    ## y_hat and MSE
    curves.Kraus1.hat   <- my_predict_fRegress(mod          = mod,
                                               xlist        = xlist.test,
                                               t.points     = t.points)
    curves.Kraus1.hat.v <- eval.fd(t.points, curves.Kraus1.hat)
    
    evalMSE    <- eval_MSE(curves     = curves.test,
                           curves.hat = curves.Kraus1.hat.v,
                           t.points   = t.points)
    MSE <- evalMSE$MSE
    MSE_reconstruction <- evalMSE$MSE_reconstruction

  }
  
  
  #### Kraus 2 --------------------------------------------------------------------
  if(any(method == 'Kraus2'))
  {
    Kraus2 <- ReconstPoFD::reconstructKraus(X_mat        = curves.train,
                                            #alpha        = 0.001894,
                                            reconst_fcts = reconst_fcts)
    
    curves.Kraus2.recon          <- Kraus2[['X_reconst_mat']]
    curves.Kraus2                <- curves.train
    curves.Kraus2[,reconst_fcts] <- curves.Kraus2.recon
    
    ## Smoothing
    basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
    fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-5)
    curves.Kraus2.fd <- smooth.basis(t.points, curves.Kraus2, fPar)$fd
    
    # ## B-list
    # fPar  <- smth$fPar
    # blist.fixed <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)
    # 
    # for(q in 1:dim(XX)[2])
    # {
    #   blist.fixed[[q]]$lambda <- lambda.fixed[[q]]
    # }
    # #blist <- blist.fixed
    # 
    # ## Selection of the optimal lambda
    # blist <- lambda_select(curves.fd     = curves.Kraus2.fd,
    #                        xlist         = xlist.train,
    #                        t.points      = T.period,
    #                        reg.idxs      = 7,
    #                        blist.fixed   = blist.fixed,
    #                        XX            = X.train,
    #                        wgts.fd       = wgt$wgts.fd,
    #                        interpolation = TRUE)
    
    ## Regression and beta estimation
    mod <- fRegress(y            = curves.Kraus2.fd,
                    xfdlist      = xlist.train,
                    betalist     = blist,
                    returnMatrix = FALSE,
                    method       = 'fRegress',
                    sep          = '.')
    ## y_hat and MSE
    curves.Kraus2.hat   <- my_predict_fRegress(mod          = mod,
                                               xlist        = xlist.test,
                                               t.points     = t.points)
    curves.Kraus2.hat.v <- eval.fd(t.points, curves.Kraus2.hat)
    
    evalMSE    <- eval_MSE(curves     = curves.test,
                           curves.hat = curves.Kraus2.hat.v,
                           t.points   = t.points)
    MSE <- evalMSE$MSE
    MSE_reconstruction <- evalMSE$MSE_reconstruction
    
  }
  
  
  
  
  #### Kneip-Liebl No Alignment -----------------------------------------------------------
  if(any(method == 'KLNoAl'))
  {
    Y_list <- lapply(seq_len(ncol(curves.train)), function(i) curves.train[!is.na(curves.train[,i]),i])
    U_list <- lapply(seq_len(ncol(curves.train)), function(i) t.points[!is.na(curves.train[,i])])
    
    KL_NoAl_CES <- my_reconstructKneipLiebl(Ly           = Y_list,
                                            Lu           = U_list,
                                            method       = 'Error>0_AlignNO_CEscores',
                                            K            = 2,
                                            reconst_fcts = reconst_fcts,
                                            nRegGrid     = NULL,
                                            maxbins      = NULL,
                                            progrbar     = FALSE)
    curves.KL.recon <- matrix(unlist(KL_NoAl_CES[['Y_reconst_list']]),
                              nrow = length(t.points), ncol = length(reconst_fcts))
    curves.KL       <- curves.train
    curves.KL[,reconst_fcts] <- curves.KL.recon
    
    ## Smoothing
    basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
    fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-5)
    curves.KL.fd <- smooth.basis(t.points, curves.KL, fPar)$fd
    
    # ## B-list
    # fPar  <- smth$fPar
    # blist.fixed <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)
    # 
    # for(q in 1:dim(XX)[2])
    # {
    #   blist.fixed[[q]]$lambda <- lambda.fixed[[q]]
    # }
    # #blist <- blist.fixed
    # 
    # ## Selection of the optimal lambda
    # blist <- lambda_select(curves.fd     = curves.KL.fd,
    #                        xlist         = xlist.train,
    #                        t.points      = T.period,
    #                        reg.idxs      = 7,
    #                        blist.fixed   = blist.fixed,
    #                        XX            = X.train,
    #                        wgts.fd       = wgt$wgts.fd,
    #                        interpolation = TRUE)
    
    ## Regression and beta estimation
    mod <- fRegress(y            = curves.KL.fd,
                    xfdlist      = xlist.train,
                    betalist     = blist,
                    returnMatrix = FALSE,
                    method       = 'fRegress',
                    sep          = '.')
    
    ## y_hat and MSE
    curves.KL.hat   <- my_predict_fRegress(mod          = mod,
                                           xlist        = xlist.test,
                                           t.points     = t.points)
    curves.KL.hat.v <- eval.fd(t.points, curves.KL.hat)
    
    evalMSE    <- eval_MSE(curves     = curves.test,
                           curves.hat = curves.KL.hat.v,
                           t.points   = t.points)
    MSE <- evalMSE$MSE
    MSE_reconstruction <- evalMSE$MSE_reconstruction
    
  }
  
  
  #### Kneip-Liebl Yes Alignment -----------------------------------------------------------
  if(any(method == 'KLAl'))
  {
    Y_list <- lapply(seq_len(ncol(curves.train)), function(i) curves.train[!is.na(curves.train[,i]),i])
    U_list <- lapply(seq_len(ncol(curves.train)), function(i) t.points[!is.na(curves.train[,i])])
    KL_Al_CES <- my_reconstructKneipLiebl(Ly             = Y_list,
                                          Lu             = U_list,
                                          method         = 'Error>0_AlignYES_CEscores',
                                          K              = 2,
                                          reconst_fcts   = reconst_fcts,
                                          nRegGrid       = NULL,
                                          maxbins        = NULL,
                                          progrbar       = FALSE)
    curves.KLAL.recon <- matrix(unlist(KL_Al_CES[['Y_reconst_list']]),
                                nrow = length(t.points), ncol = length(reconst_fcts))
    curves.KLAL       <- curves.train
    curves.KLAL[,reconst_fcts] <- curves.KLAL.recon
    
    ## Smoothing
    basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
    fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-5)
    curves.KLAL.fd <- smooth.basis(t.points, curves.KLAL, fPar)$fd
    
    # ## B-list
    # fPar  <- smth$fPar
    # blist.fixed <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)
    # 
    # for(q in 1:dim(XX)[2])
    # {
    #   blist.fixed[[q]]$lambda <- lambda.fixed[[q]]
    # }
    # #blist <- blist.fixed
    # 
    # ## Selection of the optimal lambda
    # blist <- lambda_select(curves.fd     = curves.KLAL.fd,
    #                        xlist         = xlist.train,
    #                        t.points      = T.period,
    #                        reg.idxs      = 7,
    #                        blist.fixed   = blist.fixed,
    #                        XX            = X.train,
    #                        wgts.fd       = wgt$wgts.fd,
    #                        interpolation = TRUE)
    
    
    ## Regression and beta estimation
    mod <- fRegress(y            = curves.KLAL.fd,
                    xlist        = xlist.train,
                    betalist     = blist,
                    returnMatrix = FALSE,
                    method       = 'fRegress',
                    sep          = '.')
    
    ## y_hat and MSE
    curves.KLAL.hat   <- my_predict_fRegress(mod          = mod,
                                             xlist        = xlist.test,
                                             t.points     = t.points)
    curves.KLAL.hat.v <- eval.fd(t.points, curves.KLAL.hat)
    
    evalMSE    <- eval_MSE(curves     = curves.test,
                           curves.hat = curves.KLAL.hat.v,
                           t.points   = t.points)
    MSE <- evalMSE$MSE
    MSE_reconstruction <- evalMSE$MSE_reconstruction
    
  }
  
  ## Output --------------------------------------------------------------
  out.list <- list(MSE                = MSE,
                   MSE_reconstruction = MSE_reconstruction)
  
  return(out.list)
  
}