workflow_weighted_analysis <- function(b,
                                       B,
                                       t.points,
                                       breaks,
                                       T_hp,
                                       curves,
                                       curves.true.fd,
                                       xlist,
                                       method = c('Kraus', 'KLNoAl' ,'KLAl', NULL),
                                       fix.par = 10,
                                       wgts.flag = TRUE)
{
  
  ## Methods
  source('methods/find_obs_inc.R')
  source('Simulation/methods/create_weights.R')
  source('Simulation/methods/create_zero_weights.R')
  source('methods/wt_bsplinesmoothing.R')
  source('methods/Reconstruction/my_reconstructKneipLiebl.R')
  
  source('methods/Regression/weighted_fRegress.R')
  source('methods/Regression/my_predict_fRegress.R')
  source('methods/Regression/eval_MSE.R')
  source('Simulation/methods/eval_MSE_functional.R')
  
  ## Utilities for the identification of the b-th batch
  n <- dim(curves)[2]
  p <- length(xlist)
  #### Separate training and test set ------------------------------------------
  n.test <- floor(n/B)
  idxs <- ((b-1)*n.test + 1):min((b*n.test),n)
  
  test <- rep(FALSE, n)
  test[idxs] <- TRUE
  train <- (!test)
  
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

  #### Reconstruction ----------------------------------------------------------
  ## Kraus method
  if(any(method == 'Kraus'))
  {
    reconstruction <- ReconstPoFD::reconstructKraus(X_mat        = curves.train,
                                                    alpha        = NULL,
                                                    reconst_fcts = reconst_fcts)
    curves.recon <- reconstruction[['X_reconst_mat']]
  }
  
  ## Kneip-Liebl Yes Alignment
  if(any(method == 'KLAl'))
  {
    Y_list <- lapply(seq_len(ncol(curves.train)), function(i) curves.train[!is.na(curves.train[,i]),i])
    U_list <- lapply(seq_len(ncol(curves.train)), function(i) t.points[!is.na(curves.train[,i])])
    
    reconstruction <- my_reconstructKneipLiebl(Ly             = Y_list,
                                               Lu             = U_list,
                                               method         = 'Error>0_AlignYES_CEscores',
                                               K              = NULL,
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

  }
  
  ## Kneip-Liebl No Alignment
  if(any(method == 'KLNoAl'))
  {
    Y_list <- lapply(seq_len(ncol(curves.train)), function(i) curves.train[!is.na(curves.train[,i]),i])
    U_list <- lapply(seq_len(ncol(curves.train)), function(i) t.points[!is.na(curves.train[,i])])
    
    reconstruction <- my_reconstructKneipLiebl(Ly           = Y_list,
                                               Lu           = U_list,
                                               method       = 'Error>0_AlignNO_CEscores',
                                               K            = NULL,
                                               reconst_fcts = reconst_fcts,
                                               nRegGrid     = NULL,
                                               maxbins      = NULL,
                                               progrbar     = FALSE)
    
    curves.recon <- matrix(unlist(reconstruction[['Y_reconst_list']]),
                              nrow = length(t.points), ncol = length(reconst_fcts))
    if(wgts.flag)
    {
      for(ii in 1:dim(curves.recon)[2])
      {
        curves.recon[!is.na(curves.train[,reconst_fcts[ii]]),ii] <- curves.train[!is.na(curves.train[,reconst_fcts[ii]]),reconst_fcts[ii]]
      }
    }
    
  }
  
  curves.train.rec <- curves.train
  curves.train.rec[,reconst_fcts] <- curves.recon
  
  # pdf(file = paste0("useful-pics/po-curve-and-reconstruction.pdf"), width = 8, height = 5)
  # par(mar=c(4.5, 4.5, 2.5, 1)+.1)
  # plot(t.points, curves.train[, reconst_fcts[1]], type='l', lwd=3, col='darkblue',
  #      main="(a) Partially observed curve", ylim=c(-3.5,-1), xlab="t", ylab="y(t)",
  #      cex.main=1.8, cex.axis=1.8, cex.lab=1.8)
  # lines(t.points[t.points>=1.75], curves.train.rec[(t.points>=1.75), reconst_fcts[1]],
  #       col='lightblue', lwd=3)
  # points(t.points, curves.train[, reconst_fcts[1]],pch=16, col='darkblue')
  # grid()
  # dev.off()
  
  ## Construction of the weights -----------------------------------------------
  if(is.numeric(fix.par))
  {
    wgt       <- create_weights(curves.rec    = curves.train.rec,
                                t.points      = t.points,
                                breaks        = breaks,
                                fix.par       = fix.par,
                                reconst_fcts  = reconst_fcts,
                                Thp           = T_hp.train)
    
    # plot(wgt$wgts.fd, ylim=c(0,1))
    # title(main=paste0('Parameter ', fix.par))
  } else {
    wgt       <- create_zero_weights(curves.rec    = curves.train.rec,
                                     t.points      = t.points,
                                     breaks        = breaks,
                                     reconst_fcts  = reconst_fcts,
                                     Thp           = T_hp.train)
  }

  ## Smoothing -----------------------------------------------------------------
  if(wgts.flag) # && is.numeric(fix.par)
  {
    smth      <- wt_bsplinesmoothing(curves   = curves.train.rec,
                                     wgts.obs = wgt$wgts.obs,
                                     t.points = t.points,
                                     breaks   = breaks,
                                     lambda   = NULL)
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

  ## B-list
  blist <- list(fPar,fPar,fPar)
  
  ## Regression and beta estimation --------------------------------------------
  if(wgts.flag)
  {
    mod <- weighted_fRegress(y            = curves.fd,
                             xfdlist      = xlist.train,
                             betalist     = blist,
                             wgts         = wgt$wgts.fd)
  } else {
    mod <- fRegress(y            = curves.fd,
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
  
  
  ## Evaluation of the MSE as integral mean of the squared functional errors
  curves.true.fd.test <- curves.true.fd
  curves.true.fd.test$coefs <- as.matrix(curves.true.fd$coefs[,test], ncol=n.test)

  MSE <- eval_MSE_functional(curves     = curves.true.fd.test,
                             curves.hat = curves.hat,
                             t.points   = t.points)

  
  ## Output --------------------------------------------------------------
  out.list <- list(MSE                = MSE,
                   beta_estimates     = beta_estimates)
  
  return(out.list)
  
}