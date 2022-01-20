my_PACE_fun <- function(mu, argvals, scoresP, efunctionsP, K=NULL){
  ##
  K        <- min(ncol(efunctionsP), K)
  ##
  reconstr <- unname(c( t(as.matrix(mu)) + t(scoresP[1:K]) %*% t(efunctionsP[,1:K]) ))
  ##
  ## ######################
  return(list("y_reconst"  = c(reconstr), "x_reconst"  = c(argvals)))
  ## ######################
}