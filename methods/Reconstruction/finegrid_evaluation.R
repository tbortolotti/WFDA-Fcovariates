finegrid_evaluation <- function(curves,t.points,lambda=NULL, tfine)
{
  
  n      <- dim(curves)[2]
  breaks <- t.points
  basis  <- create.bspline.basis(rangeval=range(t.points),
                                 breaks=breaks, norder=4)
  
  if(is.null(lambda)) # perform a GCV procedure to select lambda.opt
  {
    esp        <- seq(0,7, by=1)
    lambda.vec <- 10^-esp
    gcv.mat    <- matrix(data=0, nrow=dim(curves)[2], ncol=length(lambda.vec))
    for(j in 1:length(lambda.vec))
    {
      functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
      for(i in 1:n)
      {
        curve          <- curves[,i]
        M_bool_vec     <- is.na(curve)
        O_bool_vec     <- !M_bool_vec
        gcv.mat[i,j]   <- smooth.basis(t.points[O_bool_vec], curve[O_bool_vec], functionalPar)$gcv
      }
    }
    gcv <- colSums(gcv.mat)
    
    lambda.opt <- lambda.vec[which(gcv == min(gcv))]
    
  } else {
    lambda.opt = lambda
  }
  functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
  
  ## Fine grid evaluation of the curves
  curves.fine <- matrix(nrow=length(tfine), ncol=n)
  for(i in 1:n)
  {
    curve            <- curves[,i]
    M_bool_vec       <- is.na(curve)
    O_bool_vec       <- !M_bool_vec
    curve.s          <- smooth.basis(t.points[O_bool_vec], curve[O_bool_vec], functionalPar)
    curves.fine[,i]  <- eval.fd(tfine, curve.s$fd, Lfdobj=0, returnMatrix=FALSE)
    t.last           <- tail(t.points[O_bool_vec], n=1)
    curves.fine[which(tfine > t.last),i] <- NA
  }
  
  ## Output

  return(curves.fine)
  
}