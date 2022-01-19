#' Reconstruct functional data through interpolation
#'
#' This function allows to reconstruct partially observed functional data through an 
#' interpolation. Each incomplete functional datum, observed up to a sample point T*,
#' is continuously extended with a straight line with slope the mean slope that all other
#' complete observations have from T* up to the right extrema of the complete observation
#' interval.
#' @param  curves   pxn matrix, p: number of discretization points, n=number of functions
#' @export out      a list containing the original and the reconstructed curves


interpolation <- function(curves, t.points, T_hp, reconst_fcts)
{

  ## Evaluation of mean slope -------------------------------------------------
  util       <- curves[,-reconst_fcts]
  T.minor    <- sort(unique(T_hp[reconst_fcts]))
  mean.slope <- numeric(length=length(T.minor))
  Tf.idx     <- length(t.points)
  Tf         <- t.points[Tf.idx]

  for(j in 1:length(T.minor))
  {
    Ti.idx        <- tail(which(t.points<=T.minor[j]), n=1)
    Ti            <- t.points[Ti.idx]
    mean.slope[j] <- mean(util[Tf.idx,]-util[Ti.idx,])/(Tf-Ti)
  }
  
  ## All observations
  curves.rec <- curves
  
  for(i in 1:length(reconst_fcts))
  {
    Y         <- curves[,reconst_fcts[i]]
    t.lastreg <- T_hp[reconst_fcts[i]]
    slope     <- mean.slope[which(T.minor==t.lastreg)]
    i.lastobs <- tail(which(t.points<=t.lastreg),n=1)
    t.lastobs <- t.points[i.lastobs]
    Y.lastobs <- Y[i.lastobs]
    Y.rec     <- Y
    for(j in 1:(length(t.points)))
    {
      t <- t.points[j]
      if(t > t.lastobs)
      {
        Y.rec[j] <- Y.lastobs + slope*(t-t.lastobs)
      }
    }
    curves.rec[,reconst_fcts[i]] <- Y.rec 
  }

  ## Output ------------------------------------------------------------------
  
  out <- list(curves, curves.rec)
  names(out) <- c('curves', 'curves.rec')
  
  return(out)
  
}

