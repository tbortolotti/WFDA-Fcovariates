fit_ITA18 <- function(data.f, curves)
{
  # Data -----------------------------------------------------------------------
  dJB      <- data.f$dJB
  MAG      <- data.f$MAG
  SoF      <- data.f$SoF
  VS30     <- data.f$VS30
  
  q <- 9
  
  # Mh, Mref, h ----------------------------------------------------------------
  load('DATA/ITA18_regressors.RData')
  Mh.vec   <- ITA18.regressors$Mh.vec
  Mref.vec <- ITA18.regressors$Mref.vec
  h.vec    <- ITA18.regressors$h.vec
  
  n <- dim(curves)[2]
  N <- dim(curves)[1]
  
  coefs.ITA18 <- as.data.frame(matrix(data=0, nrow=N, ncol=q))
  y.hat.ITA18 <- as.data.frame(matrix(data=NA, nrow=N, ncol=n))
  sigma.ITA18 <- numeric(N)
  
  # for every period, repeat
  for(t in 1:N) #t=1
  {
    y        <- curves[t,]
    y.in     <- (!is.na(y))
    n.in     <- sum(y.in)
    
    M.h.T    <- Mh.vec[t]
    M.ref.T  <- Mref.vec[t]
    h.T      <- h.vec[t]
    
    # generate dataset for the ITA18 model at time t
    R           <- sqrt(dJB[y.in]^2 + h.T^2)
    reg.M_l     <- ifelse(MAG[y.in]<=M.h.T, MAG[y.in] - M.h.T, 0)
    reg.M_h     <- ifelse(MAG[y.in]>M.h.T, MAG[y.in] - M.h.T, 0)
    reg.SS      <- ifelse(SoF[y.in]=="SS", 1, 0)
    reg.TF      <- ifelse(SoF[y.in]=="TF", 1, 0)
    reg.D_1     <- (MAG[y.in] - M.ref.T)*log10(R)
    reg.D_2     <- log10(R)
    reg.D_3     <- R
    reg.S       <- ifelse(VS30[y.in]<=1500, log10(VS30[y.in]/800), log10(1500/800))
    
    mod.ITA18.T <- lm(y[y.in] ~ reg.M_l + reg.M_h + reg.SS + reg.TF + reg.D_1 + reg.D_2 +
                      reg.D_3 + reg.S)
    
    coefs.ITA18[t,] <- mod.ITA18.T$coefficients
    y.hat.ITA18[t,y.in] <- mod.ITA18.T$fitted.values
    
    res <- y[y.in] - y.hat.ITA18[t,y.in]
    rss <- sum(res^2)
    delta1 <- n.in-q-1
    sigma.ITA18[t] <- sqrt(rss/delta1)
    
  }
  names(coefs.ITA18) <- c('a0', 'b1', 'b2', 'f1', 'f2', 'c1', 'c2', 'c3', 'k0')
  
  outlist <- list(coefs.ITA18 = coefs.ITA18,
                  sigma.ITA18 = sigma.ITA18,
                  y.hat.ITA18 = y.hat.ITA18)
  
  return(outlist)
}