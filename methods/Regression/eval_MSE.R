eval_MSE <- function(curves, curves.hat, t.points)
{
  ## Utilities
  n <- dim(curves)[2]
  
  ## MSE_b
  ## per ogni funzione nel test
  ## valuto la distanza di quella predetta da quella vera, poi integro
  ## ma soltanto nella parte di dominio dove questa funzione è osservata
  MSE <- 0
  MSE_reconstruction <- 0
  
  range_reconstruction <- c(log10(5),log10(10))
  
  n.complete <- 0
  
  for(j in 1:n)
  {
    curve       <- curves[,j]
    M_bool_vec  <- is.na(curve)
    O_bool_vec  <- !M_bool_vec
    
    t.O         <- t.points[O_bool_vec]
    curve.O     <- curve[O_bool_vec]
    
    curve.hat.O <- curves.hat[O_bool_vec,j]
    
    err.O       <- curve.O - curve.hat.O
    
    basis    <- create.bspline.basis(rangeval=range(t.O), breaks=t.O, norder=4)
    fPar     <- fdPar(fdobj=basis, Lfdobj=2, lambda=1e-5)
    err.O.fd <- smooth.basis(t.O, err.O, fPar)$fd

    MSE_j    <- inprod(err.O.fd, err.O.fd, 0, 0, rng=range(t.O))
    
    Tf <- tail(t.O, n=1)
    Ti <- head(t.O, n=1)
    
    delta.T <- Tf - Ti
    MSE     <- MSE + MSE_j/delta.T
    
    if(sum(O_bool_vec)==length(t.points))
    {
      MSE_reconstruction_j <- inprod(err.O.fd, err.O.fd, 0, 0, rng=range_reconstruction)
      MSE_reconstruction <- MSE_reconstruction + MSE_reconstruction_j/(range_reconstruction[2]-range_reconstruction[1])
      n.complete <- n.complete + 1
    }
  }
  
  MSE <- MSE/n
  MSE_reconstruction <- MSE_reconstruction/n.complete
  
  out <- list(MSE                = MSE,
              MSE_reconstruction = MSE_reconstruction)
  
  return(out)
}