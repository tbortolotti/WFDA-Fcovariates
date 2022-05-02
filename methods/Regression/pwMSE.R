pwMSE <- function(curves, curves.fd, xlist, t.points, events, blist, B,
                        wgts.fd, set.ITA18 = FALSE, data.f = NULL, wgts.flag=TRUE)
{
  # Methods ------------------------------------------------------------------
  source('methods/Regression/weighted_fRegress.R')
  source('methods/Regression/unwgt_fRegress.R')
  source('methods/Regression/my_predict_fRegress.R')
  source('methods/fit_ITA18.R')
  source('methods/predict_ITA18.R')
  
  # Utilities ------------------------------------------------------------------
  n <- dim(curves)[2]
  N <- dim(curves)[1]
  p <- length(xlist)
  
  #### Separate training and test set ------------------------------------------
  unique.events <- unique(events)
  n.events <- length(unique.events)
  n.test.events <- floor(n.events/B)
  
  set.seed(140996)
  shuffled.unique.events <- sample(unique.events, n.events)
  
  MSE_eval <- matrix(data=0, nrow=N, ncol=B)
  MSE_ita18 <- matrix(data=0, nrow=N, ncol=B)
  pb <- progress_bar$new(total=B, format = "  computing [:bar] :percent eta: :eta")
  for(b in 1:B)#b=1
  {
    idxs <- ((b-1)*n.test.events + 1):min((b*n.test.events),n.events)
    test.events <- shuffled.unique.events[idxs]
    
    test <- (events %in% test.events)
    train <- (!test)
    
    n.test <- sum(test)
    n.train <- sum(train)
    
    curves.test    <- curves[,test]
    curves.train   <- curves[,train]
    
    curves.fd.train <- curves.fd
    curves.fd.train$coefs <- curves.fd$coefs[,train]
    
    wgts.fd.train <- wgts.fd
    wgts.fd.train$coefs <- wgts.fd$coefs[,train]
    
    data.test      <- data.f[test,]
    data.train     <- data.f[train,]
    
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
    
    ## Regression and beta estimation
    if(wgts.flag)
    {
      mod <- weighted_fRegress(y            = curves.fd.train,
                               xfdlist      = xlist.train,
                               betalist     = blist,
                               wgts         = wgts.fd.train)
    } else {
      mod <- unwgt_fRegress(y            = curves.fd.train,
                            xfdlist      = xlist.train,
                            betalist     = blist,
                            returnMatrix = FALSE,
                            method       = 'fRegress',
                            sep          = '.')
    }
    
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
    curves.hat.v <- eval.fd(t.points, curves.hat)
    
    n.T <- numeric(N)
    for(t in 1:N)
    {
      n.T[t] <- sum(!is.na(curves.test[t,]))
    }
    
    diff2.mat <- (curves.test - curves.hat.v)^2
    diff2.mat[is.na(diff2.mat)] <- 0
    MSE_eval[,b] <- rowSums(diff2.mat)/n.T
    
    
    ##ITA18
    if(set.ITA18)
    {
      coefs.ITA18      <- fit_ITA18(data.train, curves.train)$coefs.ITA18
      curves.hat.ita18 <- predict_ITA18(data.test, coefs.ITA18)
      diff2.mat.18     <- (curves.test - curves.hat.ita18)^2
      diff2.mat.18[is.na(diff2.mat.18)] <- 0
      MSE_ita18[,b]    <- rowSums(diff2.mat.18)/n.T
    }
    
    pb$tick()
    
  }

  MSE_t <- rowSums(MSE_eval)/B
  MSE_t18 <- rowSums(MSE_ita18)/B
  

  out <- list(MSE_t   = MSE_t,
              MSE_t18 = MSE_t18)

  
  
  return(out)
}