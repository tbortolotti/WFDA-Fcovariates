#' Reconstruct functional data through extrapolation
#'
#' This function allows to reconstruct partially observed functional data through an 
#' extrapolation. Each incomplete functional datum, observed up to a sample point T*,
#' is continuously extended with a straight line with slope the mean slope that all other
#' complete observations have from T* up to the right extrema of the complete observation
#' interval.
#' @param  curves   pxn matrix, p: number of discretization points, n=number of functions
#' @export out      a list containing the original and the reconstructed curves


extrapolation <- function(curves, t.points, T_hp, reconst_fcts)
{
  
  N <- length(t.points)
  ## Evaluation of mean slope --------------------------------------------------
  mean.slope <- numeric(N-1)
  for(j in 1:(N-1))
  {
    mean.slope[j] <- mean(curves[N,-reconst_fcts]-curves[j,-reconst_fcts])/(t.points[N]-t.points[j])
  }

  ## All observations
  curves.rec <- curves
  
  for(i in 1:length(reconst_fcts))
  {
    Y         <- curves[,reconst_fcts[i]]
    id.na     <- is.na(Y)
    t.lastobs <- tail(t.points[!id.na], n=1)
    slope     <- tail(mean.slope[!id.na], n=1)
    Y.lastobs <- tail(Y[!id.na], n=1)
    curves.rec[id.na,reconst_fcts[i]] <- Y.lastobs + slope*(t.points[id.na]-t.lastobs)
  }

  ## Output --------------------------------------------------------------------
  
  out <- list(curves, curves.rec)
  names(out) <- c('curves', 'curves.rec')
  
  return(out)
  
}

