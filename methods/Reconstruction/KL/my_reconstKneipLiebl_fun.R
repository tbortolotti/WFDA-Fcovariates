my_reconstKneipLiebl_fun <- function(mu, argvals, argvalsO, scoresO, efun_reconst, fragmO=NULL, K=NULL){
  ##
  K        <- min(ncol(efun_reconst), K)
  locO     <- match(argvalsO, argvals)
  ##
  reconstr <- unname(c( t(as.matrix(mu)) + t(scoresO[1:K, drop=FALSE]) %*% t(efun_reconst[,1:K, drop=FALSE]) ))
  ##
  if( !is.null(fragmO) ){
    ## Aligning the reconstructed parts to fragm0
    reconstr[1:min(locO)]               <- reconstr[1:min(locO)]               + fragmO[1]              - reconstr[min(locO)]
    reconstr[max(locO):length(argvals)] <- reconstr[max(locO):length(argvals)] + fragmO[length(fragmO)] - reconstr[max(locO)]
    reconstr[locO]                      <- fragmO
  }
  ##
  ## ######################
  return(list("y_reconst"  = c(reconstr), "x_reconst"  = c(argvals)))
  ## ######################
}
