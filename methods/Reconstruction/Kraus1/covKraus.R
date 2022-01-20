covKraus <- function(X_mat){
  p <- nrow(X_mat)
  n <- ncol(X_mat)
  ##
  X_cent_mat  <- X_mat - rowMeans(X_mat, na.rm = TRUE)
  ##
  covKraus_mat <- matrix(NA, ncol = p, nrow = p)	
  for(s in seq(1, p)){
    for(t in seq(s, p)){
      X_cent_s  <- X_cent_mat[s, ]
      X_cent_t  <- X_cent_mat[t, ]
      n_na      <- sum(is.na(c(X_cent_s * X_cent_t)))
      if(n-n_na == 0){
        covKraus_mat[s,t] <- NA
      }else{
        covKraus_mat[s,t] <- mean(X_cent_s * X_cent_t, na.rm = TRUE)
      }
      covKraus_mat[t,s] <- covKraus_mat[s,t]
    }
  }
  return(covKraus_mat)
}