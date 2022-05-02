workflow_weighted_smoothing <- function(t.points,
                                        breaks,
                                        T_hp,
                                        curves,
                                        curves.true.fd,
                                        method = c('Kraus', 'KLNoAl' ,'KLAl', NULL),
                                        fix.par = 10,
                                        wgts.flag = TRUE)
{
  
  ## Methods
  source('methods/find_obs_inc.R')
  source('Simulation/methods/create_weights_simulation.R')
  source('Simulation/methods/create_old_weights.R')
  source('Simulation/methods/create_zero_weights.R')
  source('methods/wt_bsplinesmoothing.R')
  source('methods/Reconstruction/my_reconstructKneipLiebl.R')
  
  source('Simulation/methods/eval_MSE_functional.R')
  
  ## Utilities for the identification of the b-th batch
  n <- dim(curves)[2]

  #### Separate training and test set ------------------------------------------
  reconst_fcts   <- find_obs_inc(Y = curves)
  
  #### Reconstruction ----------------------------------------------------------
  ## Kraus method
  if(any(method == 'Kraus'))
  {
    reconstruction <- ReconstPoFD::reconstructKraus(X_mat        = curves,
                                                    alpha        = NULL,
                                                    reconst_fcts = reconst_fcts)
    curves.recon <- reconstruction[['X_reconst_mat']]
  }
  
  ## Kneip-Liebl Yes Alignment
  if(any(method == 'KLAl'))
  {
    Y_list <- lapply(seq_len(ncol(curves)), function(i) curves[!is.na(curves[,i]),i])
    U_list <- lapply(seq_len(ncol(curves)), function(i) t.points[!is.na(curves[,i])])
    
    reconstruction <- my_reconstructKneipLiebl(Ly             = Y_list,
                                               Lu             = U_list,
                                               method         = 'Error>0_AlignYES_CEscores',
                                               K              = NULL,
                                               reconst_fcts   = reconst_fcts,
                                               nRegGrid       = NULL,
                                               maxbins        = NULL,
                                               progrbar       = FALSE)
    
    curves.recon <- matrix(unlist(reconstruction[['Y_reconst_list']]),
                           nrow = length(t.points), ncol = length(reconst_fcts))
    if(wgts.flag)
    {
      for(ii in 1:dim(curves.recon)[2])
      {
        curves.recon[!is.na(curves[,reconst_fcts[ii]]),ii] <- curves[!is.na(curves[,reconst_fcts[ii]]),reconst_fcts[ii]]
      }
    }
    
  }
  
  ## Kneip-Liebl No Alignment
  if(any(method == 'KLNoAl'))
  {
    Y_list <- lapply(seq_len(ncol(curves)), function(i) curves[!is.na(curves[,i]),i])
    U_list <- lapply(seq_len(ncol(curves)), function(i) t.points[!is.na(curves[,i])])
    
    reconstruction <- my_reconstructKneipLiebl(Ly           = Y_list,
                                               Lu           = U_list,
                                               method       = 'Error>0_AlignNO_CEscores',
                                               K            = NULL,
                                               reconst_fcts = reconst_fcts,
                                               nRegGrid     = NULL,
                                               maxbins      = NULL,
                                               progrbar     = FALSE)
    
    curves.recon <- matrix(unlist(reconstruction[['Y_reconst_list']]),
                           nrow = length(t.points), ncol = length(reconst_fcts))
    if(wgts.flag)
    {
      for(ii in 1:dim(curves.recon)[2])
      {
        curves.recon[!is.na(curves[,reconst_fcts[ii]]),ii] <- curves[!is.na(curves[,reconst_fcts[ii]]),reconst_fcts[ii]]
      }
    }
    
  }
  
  curves.rec <- curves
  curves.rec[,reconst_fcts] <- curves.recon
  
  ## Construction of the weights -----------------------------------------------
  if(is.numeric(fix.par))
  {
    wgt       <- create_old_weights(curves.rec    = curves.rec,
                                    t.points      = t.points,
                                    breaks        = breaks,
                                    fix.par       = fix.par,
                                    reconst_fcts  = reconst_fcts,
                                    Thp           = T_hp)

  } else if (fix.par=="0-wgts") {
    wgt       <- create_weights_simulation(curves.rec    = curves.rec,
                                           t.points      = t.points,
                                           breaks        = breaks,
                                           fix.par       = 100,
                                           reconst_fcts  = reconst_fcts,
                                           Thp           = T_hp)
  } else {
    wgt       <- create_zero_weights(curves.rec    = curves.rec,
                                     t.points      = t.points,
                                     breaks        = breaks,
                                     reconst_fcts  = reconst_fcts,
                                     Thp           = T_hp)
  }
  
  ## Smoothing -----------------------------------------------------------------
  if(wgts.flag) # && is.numeric(fix.par)
  {
    smth      <- wt_bsplinesmoothing(curves   = curves.rec,
                                     wgts.obs = wgt$wgts.obs,
                                     t.points = t.points,
                                     breaks   = breaks,
                                     lambda   = NULL)
    curves.fd <- smth$curves.fd
    fPar  <- smth$fPar
  } else {
    basis <- create.bspline.basis(rangeval=range(t.points), breaks=breaks, norder=4)
    esp        <- seq(-7,1, by=1)
    lambda.vec <- sort(10^esp)
    gcv.vec    <- numeric(length(lambda.vec))
    for(j in 1:length(lambda.vec))
    {
      fPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.vec[j])
      gcv.vec[j] <- sum(smooth.basis(t.points, curves.rec, fPar)$gcv)/n
    }
    lambda.opt <- lambda.vec[which(gcv.vec == min(gcv.vec))]
    fPar  <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda.opt)
    curves.fd <- smooth.basis(t.points, curves.rec, fPar)$fd
  }
  
  MSE <- eval_MSE_functional(curves     = curves.true.fd,
                             curves.hat = curves.fd,
                             t.points   = t.points)
  
  
  ## Output --------------------------------------------------------------
  out.list <- list(MSE = MSE)
  
  return(out.list)
  
}