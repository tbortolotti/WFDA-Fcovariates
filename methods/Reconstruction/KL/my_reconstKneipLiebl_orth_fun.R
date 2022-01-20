my_reconstKneipLiebl_orth_fun <- function(mu, argvals, CE_scores_orth , efun_reconst_orth, K=NULL){
  ##
  K        <- min(ncol(efun_reconst_orth), K)
  ##
  reconstr <- unname(c( t(as.matrix(mu)) + t(CE_scores_orth[1:K]) %*% t(efun_reconst_orth[,1:K]) ))
  ##
  ## ######################
  return(list("y_reconst"  = c(reconstr), "x_reconst"  = c(argvals)))
  ## ######################
}