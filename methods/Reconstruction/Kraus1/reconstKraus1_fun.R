reconstKraus1_fun <- function(cov_mat, X_cent_vec, alpha=1e-4, tfine){
  
  delta           <- tfine[2]-tfine[1]

  M_bool_vec      <- is.na(X_cent_vec)
  O_bool_vec      <- !M_bool_vec
  p               <- nrow(cov_mat)

  Phi             <- base::eigen(cov_mat)$vectors * delta^(-1/2)
  Phi_O           <- Phi[O_bool_vec,]
  Phi_M           <- Phi[M_bool_vec,]
  covMO_mat       <- cov_mat[M_bool_vec, O_bool_vec]
  covOM_mat       <- cov_mat[O_bool_vec, M_bool_vec]
  covOO_mat       <- cov_mat[O_bool_vec, O_bool_vec]
  covOO_a_mat     <- covOO_mat + alpha * diag(p)[O_bool_vec, O_bool_vec]
  covOO_a_mat_inv <- solve(covOO_a_mat)

  Aa_mat          <- covOO_a_mat_inv %*% covOM_mat %*% Phi_M
  
  beta.O     <- list()
  beta.M     <- list()
  beta.list  <- list()
  X_cent_rec <- matrix(data=0, nrow=dim(Phi)[1], ncol=1)
  for(j in 1:dim(Phi)[2])
  {
    #beta.O[[j]]    <- integrate(my.fun, lower=tOrange[1], upper=tOrange[2],
    #                         X_cent_vec.fd, Phi.fd)$value
    #beta.M[[j]]    <- integrate(my.fun, lower=tOrange[1], upper=tOrange[2],
    #                         X_cent_vec.fd, Aa.fd)$value
    beta.O[[j]]    <- X_cent_vec[O_bool_vec] %*% Phi_O[,j] * delta
    beta.M[[j]]    <- X_cent_vec[O_bool_vec] %*% Aa_mat[,j] * delta
    
    beta.list[[j]] <- beta.O[[j]] + beta.M[[j]]
    X_cent_rec     <- X_cent_rec + as.numeric(beta.list[[j]])*Phi[,j]
  }

  ## compute df
  #df <- sum(diag(covOO_a_mat_inv %*% covOO_mat))
  lam00 <- base::eigen(covOO_mat)$values * delta
  lam00 <- lam00[lam00>0]
  df    <- sum(lam00/(lam00+alpha))

  return(list("X_cent_reconst_vec" = X_cent_rec,
              "beta.list"          = beta.list,
              "df"                 = df,
              "alpha"              = alpha))
}
