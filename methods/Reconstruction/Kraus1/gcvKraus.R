gcvKraus <- function(cov_mat, mean_vec, X_Compl_mat, M_bool_vec, tfine, alpha){
  
  ## Useful functions ------------------------------------------------
  n_Compl  <- ncol(X_Compl_mat)
  rss_vec  <- rep(NA,n_Compl)
  ##
  for(j in 1:n_Compl){
    X_gcv             <-  X_Compl_mat[,j]
    X_gcv[M_bool_vec] <- NA
    ##
    result_tmp <- reconstKraus1_fun(cov_mat    = cov_mat, 
                                    X_cent_vec = c(X_gcv - mean_vec), 
                                    alpha      = alpha,
                                    tfine      = tfine)
    ##
    X_fit      <- c(result_tmp[['X_cent_reconst_vec']] + mean_vec)
    ##
    rss_vec[j] <- sum((X_fit[M_bool_vec] - X_Compl_mat[M_bool_vec,j])^2)
  }
  gcv <- sum(rss_vec)/((1-result_tmp[['df']]/n_Compl)^2)
  ##
  return(gcv)
}
