#### Construction of Sphi, dim(Sphi) = nL * nT

build_Sphi <- function(y2cmaps, curves)
{
  
  library(Matrix)
  
  ## Utilities --------------------------------------
  L  <- dim(y2cmaps)[1]
  TT <- dim(y2cmaps)[2]
  n  <- dim(y2cmaps)[3]
 
  # -------------------------------------------------
  
  Sphi <- Matrix(nrow = (n*L), ncol = (n*TT), data = 0, sparse = TRUE)
  #Sphi <- as(Sphi, "dgTMatrix") # by default, Matrix() returns dgCMatrix
  
  pb <- progress_bar$new(total=L, format = "  computing [:bar] :percent eta: :eta")
  for(l in 1:L)
  {
    r.idxs    <- ((l-1)*n+1):(l*n) 
    for(t in 1:TT)
    {
      c.idxs   <- ((t-1)*n+1):(t*n)
      diag(Sphi[r.idxs,c.idxs]) <- y2cmaps[l,t,]
    }
    pb$tick()
  }
  
  save(Sphi, file="DATA/Sphi.RData")
  
}
