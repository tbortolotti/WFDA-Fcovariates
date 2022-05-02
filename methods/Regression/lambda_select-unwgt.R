lambda_select_unwgt <- function(b, B, curves, curves.fd, xlist, t.points, fPar, exp.vec, events)
{
  # Methods ------------------------------------------------------------------
  source('methods/Regression/unwgt_fRegress.R')
  source('methods/Regression/my_predict_fRegress.R')
  source('methods/Regression/eval_MSE_functional.R')
  
  # Utilities ----------------------------------------------------------------
  n <- dim(curves.fd$coefs)[2]
  p <- length(xlist)
  onebasis <- create.constant.basis(range(t.points))
  onesfd   <- fd(1,onebasis)
  
  ## Create betalist ---------------------------------------------------------
  # First round
  blist      <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)
  
  # Second round
  #load('blist_options/log/MSE_unwgt_fPar2.RData')
  
  ## Create a vector of lambdas to be tested
  lambda.vec <- 10^exp.vec
  
  ## Separate training and test set -----------------------------------------
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
  
  curves.fd.test <- curves.fd
  curves.fd.test$coefs <- curves.fd$coefs[,test]
  curves.fd.train <- curves.fd
  curves.fd.train$coefs <- curves.fd$coefs[,train]
  
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
  
  ## Test
  MSE_list <- list()
  
  pb <- progress_bar$new(total=p)
  for(reg in 1:p) # reg=1
  {
    MSE_vec <- numeric(length(lambda.vec))
    for(l in 1:length(lambda.vec)) # l=1 
    {
      lambda <- lambda.vec[l]
      
      ## Change the lambda value of the corresponding functionalPar in blist
      blist[[reg]]$lambda <- lambda
      
      ## Regression and beta estimation
      mod <- unwgt_fRegress(y            = curves.fd.train,
                            xfdlist      = xlist.train,
                            betalist     = blist,
                            returnMatrix = FALSE,
                            method       = 'fRegress',
                            sep          = '.')
      
      ## Prepare the list of functional covariates which we use for prediction on
      ## test set
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
      MSE <- eval_MSE_functional(curves     = curves.fd.test,
                                 curves.hat = curves.hat,
                                 t.points   = t.points)
      
      MSE_vec[l] <- MSE
      
    }
    
    # set the lambda correspondent to the current regressor as the one minimizing the MSE
    blist[[reg]]$lambda <- lambda.vec[which(MSE_vec==min(MSE_vec))]
    
    MSE_list[[reg]] <- MSE_vec
    pb$tick()
  }
  
  return(MSE_list)
}