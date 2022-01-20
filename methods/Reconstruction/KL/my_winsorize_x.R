my_winsorize_x <- function(x, cut = 0.005){
  x <- x[!is.na(x)]
  ##
  cut_point_top    <- stats::quantile(x, 1 - cut)
  cut_point_bottom <- stats::quantile(x,     cut)
  i <-  which(x >= cut_point_top)
  j <-  which(x <= cut_point_bottom)
  x[i] <- cut_point_top
  x[j] <- cut_point_bottom
  ##
  return(x)
}