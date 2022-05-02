eval_MSE_functional <- function(curves, curves.hat, t.points)
{
  n <- dim(curves$coefs)[2]
  
  delta.T <- range(t.points)[2] - range(t.points)[1]
  MSE.vec <- numeric(n)
  for(i in 1:n) #i=1
  {
    yi <- curves
    yi$coefs <- as.matrix(curves$coefs[,i], ncol=1)
    
    yi.hat <-  curves.hat
    yi.hat$coefs <- as.matrix(curves.hat$coefs[,i], ncol=1)
    
    err.i <- yi - yi.hat
    MSE.vec[i] <- diag(inprod(err.i, err.i, 0, 0, rng=range(t.points)))
    MSE.vec[i] <- MSE.vec[i]/delta.T
  }
  
  MSE     <- sum(MSE.vec)/n

  return(MSE)
}
