predict_ITA18 <- function(data, coefs.ITA18)
{
  
  # Data ----------------------------------------------------------------------
  dJB      <- data$dJB
  MAG      <- data$MAG
  SoF      <- data$SoF
  VS30     <- data$VS30
  
  # ITA18 regressors -------------------------------------------------
  a0       <- coefs.ITA18$a0
  b1       <- coefs.ITA18$b1
  b2       <- coefs.ITA18$b2
  c1       <- coefs.ITA18$c1
  c2       <- coefs.ITA18$c2
  c3       <- coefs.ITA18$c3
  k0       <- coefs.ITA18$k0
  f1       <- coefs.ITA18$f1
  f2       <- coefs.ITA18$f2
  
  load('DATA/ITA18_regressors.RData')
  
  Mh.vec   <- ITA18.regressors$Mh.vec
  Mref.vec <- ITA18.regressors$Mref.vec
  h.vec    <- ITA18.regressors$h.vec
  
  n <- length(dJB)
  N <- length(a0)
  
  SA.ITA18 <- matrix(data=0, nrow=N, ncol=n)
  
  # for every period, repeat
  for(t in 1:N) #t=1
  {

    a0.ITA18 <- a0[t]
    b1.ITA18 <- b1[t]
    b2.ITA18 <- b2[t]
    c1.ITA18 <- c1[t]
    c2.ITA18 <- c2[t]
    c3.ITA18 <- c3[t]
    k0.ITA18 <- k0[t]
    f1.ITA18 <- f1[t]
    f2.ITA18 <- f2[t]
    
    beta.T.ITA18 <- as.matrix(c(a0.ITA18, b1.ITA18, b2.ITA18, f1.ITA18,
                                f2.ITA18, c1.ITA18, c2.ITA18, c3.ITA18, k0.ITA18))
    
    M.h.T    <- Mh.vec[t]
    M.ref.T  <- Mref.vec[t]
    h.T      <- h.vec[t]

    # generate dataset for the ITA18 model at time t
    R           <- sqrt(dJB^2 + h.T^2)
    intercept   <- rep(1, n)
    reg.M_l     <- ifelse(MAG<=M.h.T, MAG - M.h.T, 0)
    reg.M_h     <- ifelse(MAG>M.h.T, MAG - M.h.T, 0)
    reg.SS      <- ifelse(SoF=="SS", 1, 0)
    reg.TF      <- ifelse(SoF=="TF", 1, 0)
    reg.D_1     <- (MAG - M.ref.T)*log10(R)
    reg.D_2     <- log10(R)
    reg.D_3     <- R
    reg.S       <- ifelse(VS30<=1500, log10(VS30/800), log10(1500/800))
    
    XX.T        <- cbind(intercept, reg.M_l, reg.M_h, reg.SS, reg.TF, reg.D_1, reg.D_2, reg.D_3, reg.S)
    
    SA.ITA18[t,] <- XX.T%*%beta.T.ITA18
    
  }
  
  return(SA.ITA18)
}