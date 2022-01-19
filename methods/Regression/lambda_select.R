lambda_select <- function(b, B, curves, curves.fd, xlist, t.points, fPar, events, wgts.fd)
{
  # Methods ------------------------------------------------------------------
  source('methods/Regression/weighted_fRegress.R')
  source('methods/Regression/my_predict_fRegress.R')
  source('methods/Regression/eval_MSE.R')
  
  # Utilities ----------------------------------------------------------------
  n <- dim(curves.fd$coefs)[2]
  q <- length(xlist)
  onebasis <- create.constant.basis(range(t.points))
  onesfd   <- fd(1,onebasis)
  
  ## Create betalist ---------------------------------------------------------
  blist      <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar)
  
  ## Create a vector of lambdas to be tested
  exp.vec    <- seq(log10(fPar$lambda)-1, max(log10(fPar$lambda)+6,2), by=1)
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

  curves.test    <- curves[,test]
  curves.train   <- curves[,train]
  
  curves.fd.train <- curves.fd
  curves.fd.train$coefs <- curves.fd$coefs[,train]
  
  wgts.fd.train <- wgts.fd
  wgts.fd.train$coefs <- wgts.fd$coefs[,train]
  
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
  
  ## Test
  MSE_list <- list()
  
  pb <- progress_bar$new(total=q)
  for(reg in 1:q)
  {
    MSE_vec <- numeric(length(lambda.vec))
    for(l in 1:length(lambda.vec))
    {
      lambda <- lambda.vec[l]
      
      ## Change the lambda value of the corresponding functionalPar in blist
      blist[[reg]]$lambda <- lambda
      
      ## Regression and beta estimation
      mod <- weighted_fRegress(y            = curves.fd.train,
                               xfdlist      = xlist.train,
                               betalist     = blist,
                               wgts         = wgts.fd.train)
      
      ## y_hat and MSE
      curves.hat <- my_predict_fRegress(mod          = mod,
                                        xlist        = xlist.test,
                                        t.points     = t.points)
      curves.hat.v <- eval.fd(t.points, curves.hat)
      
      err2.mat <- (curves.test - curves.hat.v)^2
      err.fd   <- smooth.basis(t.points, err2.mat, fPar)$fd
      
      MSE <- inprod(err.fd, onesfd, 0, 0, rng=range(t.points))/(range(t.points)[2]-range(t.points)[1])
      MSE <- sum(MSE)/n
      
      MSE_vec[l] <- MSE
    
    }
    
    # set the lambda correspondent to the current regressor as the one minimizing the MSE
    blist[[reg]]$lambda <- lambda.vec[which(MSE_vec==min(MSE_vec))]
    
    MSE_list[[reg]] <- MSE_vec
    pb$tick()
  }
  
  return(MSE_list)
}