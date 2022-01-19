find_obs_inc <- function(Y)
{
  n            <- dim(Y)[2]
  incomplete   <- c()
  for(i in 1:n)
  {
    incomplete[i] <- any(is.na(Y[,i]))
  }
  reconst_fcts <- which(incomplete==TRUE)
  return(reconst_fcts)
}