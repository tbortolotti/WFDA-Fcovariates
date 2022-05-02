methods_workflow <- function(b,
                             T.period,
                             t.points,
                             breaks,
                             T_hp,
                             log.Thp,
                             curves,
                             events,
                             B,
                             xlist,
                             blist,
                             method = c('extrapolation', 'KLAl'),
                             fix.par=NULL,
                             wgts.flag = TRUE)
                                 
{
  
  ## Methods
  source('methods/find_obs_inc.R')
  source('methods/extrapolation.R')
  source('methods/create_zero_weights.R')
  source('methods/create_weights.R')
  source('methods/wt_bsplinesmoothing.R')
  source('methods/Reconstruction/my_reconstructKneipLiebl.R')
  
  source('methods/Regression/weighted_fRegress.R')
  source('methods/Regression/unwgt_fRegress.R')
  source('methods/Regression/my_predict_fRegress.R')
  source('methods/Regression/eval_MSE.R')
  
  ## Utilities for the identification of the b-th batch
  n      <- dim(curves)[2]
  p      <- length(xlist)
  
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
  
  n.test <- sum(test)
  
  curves.test    <- curves[,test]
  curves.train   <- curves[,train]
  
  # Create xlist.test and xlist.train
  xlist.test  <- list()
  xlist.train <- list()
  
  for(i in 1:length(xlist))
  {
    if (inherits(xlist[[i]], "fd"))
    {
      xlist.test[[i]]  <- xlist[[i]]
      xlist.test[[i]]$coefs  <- xlist[[i]]$coefs[,test]
      
      xlist.train[[i]] <- xlist[[i]]
      xlist.train[[i]]$coefs <- xlist[[i]]$coefs[,train]
      
    } else if (inherits(xlist[[i]], "numeric")) {
      
      xlist.test[[i]]  <- xlist[[i]][test]
      xlist.train[[i]] <- xlist[[i]][train]
      
    } else if (inherits(xlist[[i]], "matrix" )) {
      
      xlist.test[[i]]  <- xlist[[i]][test,1]
      xlist.train[[i]] <- xlist[[i]][train,1]
    }
  }
  
  reconst_fcts   <- find_obs_inc(Y = curves.train)
  
  T_hp.test      <- T_hp[test]
  T_hp.train     <- T_hp[train]
  
  log.Thp.test   <- log.Thp[test]
  log.Thp.train  <- log.Thp[train]
  
  #### Extrapolation -----------------------------------------------------------
  if(any(method == 'extrapolation'))
  {
    extrapolate   <- extrapolation(curves       = curves.train,
                                   t.points     = T.period,
                                   T_hp         = T_hp.train,
                                   reconst_fcts = reconst_fcts)
    curves.train.rec <- extrapolate$curves.rec
  }
  
  #### Kneip-Liebl Alignment -----------------------------------------------------------
  
  if(any(method == 'KLAl'))
  {
    load('Results/Reconstruction-methods-comparison/K_Al.RData')
    K.train <- K_Al[train]
    
    Y_list <- lapply(seq_len(ncol(curves.train)), function(i) curves.train[!is.na(curves.train[,i]),i])
    U_list <- lapply(seq_len(ncol(curves.train)), function(i) t.points[!is.na(curves.train[,i])])
    
    reconstruction <- my_reconstructKneipLiebl(Ly             = Y_list,
                                               Lu             = U_list,
                                               method         = 'Error>0_AlignYES_CEscores',
                                               K              = K.train[reconst_fcts], #a preliminar step was run to fix this term over CV iterations
                                               reconst_fcts   = reconst_fcts,
                                               nRegGrid       = NULL,
                                               maxbins        = NULL,
                                               progrbar       = FALSE)
    curves.recon <- matrix(unlist(reconstruction[['Y_reconst_list']]),
                                nrow = length(t.points), ncol = length(reconst_fcts))
    
    if(wgts.flag)
    {
      for(ii in 1:dim(curves.recon)[2])
      {
        curves.recon[!is.na(curves.train[,reconst_fcts[ii]]),ii] <- curves.train[!is.na(curves.train[,reconst_fcts[ii]]),reconst_fcts[ii]]
      }
    }
    curves.train.rec <- curves.train
    curves.train.rec[,reconst_fcts] <- curves.recon
  }
  
  ## Construction of the weights -----------------------------------------------
  if(is.numeric(fix.par))
  {
    wgt       <- create_weights(curves.rec    = curves.train.rec,
                                t.points      = t.points,
                                breaks        = breaks,
                                fix.par       = fix.par,
                                reconst_fcts  = reconst_fcts,
                                Thp           = log.Thp.train)
    
  } else {
    wgt       <- create_zero_weights(curves.rec    = curves.train.rec,
                                     t.points      = t.points,
                                     breaks        = breaks,
                                     reconst_fcts  = reconst_fcts,
                                     Thp           = log.Thp.train)
  }
  
  ## Smoothing -----------------------------------------------------------------
  if(wgts.flag)
  {
    smth      <- wt_bsplinesmoothing(curves   = curves.train.rec,
                                     wgts.obs = wgt$wgts.obs,
                                     t.points = t.points,
                                     breaks   = breaks,
                                     lambda   = 1e-5)
    curves.fd <- smth$curves.fd
    fPar  <- smth$fPar
  } else {
    basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
    esp        <- seq(-7,1, by=1)
    lambda.vec <- sort(10^esp)
    gcv.vec    <- numeric(length(lambda.vec))
    for(j in 1:length(lambda.vec))
    {
      fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
      gcv.vec[j] <- sum(smooth.basis(t.points, curves.train.rec, fPar)$gcv)/n
    }
    lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
    fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
    curves.fd <- smooth.basis(t.points, curves.train.rec, fPar)$fd
  }
  
  ## Regression and beta estimation --------------------------------------------
  if(wgts.flag)
  {
    mod <- weighted_fRegress(y            = curves.fd,
                             xfdlist      = xlist.train,
                             betalist     = blist,
                             wgts         = wgt$wgts.fd)
  } else {
    mod <- unwgt_fRegress(y            = curves.fd,
                          xfdlist      = xlist.train,
                          betalist     = blist,
                          returnMatrix = FALSE,
                          method       = 'fRegress',
                          sep          = '.')
  }
  
  beta_estimates <- mod$betaestlist
  
  ## Prepare the list of functional covariates which we use for prediction on
  ## test set
  onebasis <- create.constant.basis(range(t.points))
  onesfd   <- fd(1,onebasis)
  xfdlist.test  <- xlist.test
  xerror <- FALSE
  for (j in 1:p) {
    xfdj <- xfdlist.test[[j]]
    if (inherits(xfdj, "numeric")) {
      if (!is.matrix(xfdj)) xfdj = as.matrix(xfdj)
      Zdimj <- dim(xfdj)
      if (Zdimj[1] != n.test && Zdimj != 1) {
        print(paste("Vector in XFDLIST[[",j,"]] has wrong length."))
        xerror = TRUE 
      } 
      if (Zdimj[2] != 1) {
        print(paste("Matrix in XFDLIST[[",j,"]] has more than one column."))
        xerror = TRUE 
      } 
      xfdlist.test[[j]] <- fd(matrix(xfdj,1,n.test), onebasis)
    }
  }
  
  ## y_hat and MSE
  curves.hat   <- my_predict_fRegress(mod          = mod,
                                      xlist        = xfdlist.test,
                                      t.points     = t.points)
  curves.hat.vals <- eval.fd(t.points, curves.hat)
  
  ## MSE pointwise
  diff2 <- (curves.hat.vals - curves.test)^2
  
  # Per ogni riga della matrice, voglio sommare i quadrati e dividere per il numero di
  # osservazioni che ho su quella riga
  
  MSE_pw <- numeric(N)
  for(t in 1:N)
  {
    MSE_pw[t] <- sum(diff2[t,!is.na(diff2[t,])])/sum(!is.na(diff2[t,]))
  }
  
  MSE_glob <- numeric(n.test)
  MSE_glob_bis <- numeric(n.test)
  for(n in 1:n.test)
  {
    MSE_glob[n] <- sum(diff2[!is.na(diff2[,n]),n])/sum(!is.na(diff2[,n]))
    MSE_glob_bis[n] <- sum(diff2[!is.na(diff2[,n]),n])/(sum(!is.na(diff2[,n]))^2)
  }
  
  MSE_part <- numeric(0)
  MSE_part_bis <- numeric(0)
  t.idxs <- seq(32,37)
  for(n in 1:n.test)
  {
    if(sum(!is.na(diff2[t.idxs,n]))!=0)
    {
      temp <- sum(diff2[!is.na(diff2[t.idxs,n]),n])/sum(!is.na(diff2[t.idxs,n]))
      temp2 <- sum(diff2[!is.na(diff2[t.idxs,n]),n])/(sum(!is.na(diff2[t.idxs,n]))^2)
      MSE_part <- c(MSE_part, temp)
      MSE_part_bis <- c(MSE_part_bis, temp2)
    }
  }
  ## Output --------------------------------------------------------------
  out.list <- list(MSE_pw         = MSE_pw,
                   MSE_glob       = MSE_glob,
                   MSE_glob_bis   = MSE_glob_bis,
                   MSE_part       = MSE_part,
                   MSE_part_bis   = MSE_part_bis,
                   beta_estimates = beta_estimates)
  
  return(out.list)
}