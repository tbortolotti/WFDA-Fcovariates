#' Reconstruct partially observed functions
#'
#' This function allows you to reconstruct the missing parts of a function given the observed parts.
#' @param Ly           List of Y-values. The ith (i=1,...,n) list-element contains \eqn{Y_{i1},\dots,Y_{im}}{Y_{i1},...,Y_{im}}
#' @param Lu           List of U-values. The ith (i=1,...,n) list-element contains \eqn{U_{i1},\dots,U_{im}}{U_{i1},...,U_{im}}
#' @param reconst_fcts A vector specifying the list elements in Ly which need to be reconstructed. Default (reconst_fcts=NULL) will reconstruct all functions.
#' @param method       One of the following options: 'Error=0_AlignYES_CommonGrid', 'Error>0_AlignYES', 'Error>=0_AlignNO', 'Error>0_AlignYES_CEscores', 'Error>0_AlignNO_CEscores', 'PACE', 'Error=0_PACE'.
#' @param K            Truncation parameter. If K=NULL (default), K is determined using GCV.
#' @param nRegGrid     Number of gridpoints within [a,b] used for the reconstruction result.
#' @param maxbins      If maxbins=NULL (default), maxbins is set to 1000. For speeding up simulations, use, for instance maxbins=200.
#' @param progrbar     Show progress bar (TRUE) or not (FALSE, default)
#' @export reconstructKneipLiebl
#' @examples  
#' 
#' a <- 0; b <- 1
#' set.seed(223109)
#' 
#' ## Generate partially observed functional data with error
#' SimDat        <- simuldata(n = 50, m=15, a = a, b = b, DGP="DGP1")
#' ## 
#' Y_list   <- SimDat[['Y_list']]; Y_mat <- SimDat[['Y_mat']]
#' U_list   <- SimDat[['U_list']]; U_mat <- SimDat[['U_mat']]
#' ##
#' ## Reconstruction with alignments of reconstructed parts and presmoothed fragments
#' ## using integral scores
#' reconst_result_1 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>0_AlignYES', reconst_fcts = 1:3)
#' Y_reconst_mat_1  <- matrix(unlist(reconst_result_1[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_1  <- matrix(unlist(reconst_result_1[['U_reconst_list']]), ncol=3) 
#' ##
#' ## Reconstruction with out alignments
#' ## using integral scores
#' reconst_result_2 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>=0_AlignNO', reconst_fcts = 1:3)
#' Y_reconst_mat_2  <- matrix(unlist(reconst_result_2[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_2  <- matrix(unlist(reconst_result_2[['U_reconst_list']]), ncol=3) 
#' ##
#' par(mfrow=c(1,3))
#' matplot(x=U_mat[,1:3], y=Y_mat[,1:3], ylab="", col=gray(.5), type="l", 
#' main="Orig. Data", xlim=c(a,b))
#' matplot(x=U_reconst_mat_1, y=Y_reconst_mat_1, col=gray(.5), 
#' type="l", main="With Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2) 
#' matplot(x=U_reconst_mat_2, y=Y_reconst_mat_2, col=gray(.5), 
#' type="l", main="Without Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2)
#' par(mfrow=c(1,1))
#' 
#' 
#' ## Reconstruction with alignments of reconstructed parts and presmoothed fragments
#' ## using conditional expectation scores
#' reconst_result_1 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>0_AlignYES_CEscores', reconst_fcts = 1:3)
#' Y_reconst_mat_1  <- matrix(unlist(reconst_result_1[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_1  <- matrix(unlist(reconst_result_1[['U_reconst_list']]), ncol=3) 
#' ##
#' ## Reconstruction with out alignments
#' ## using conditional expectation scores
#' reconst_result_2 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>0_AlignNO_CEscores', reconst_fcts = 1:3)
#' Y_reconst_mat_2  <- matrix(unlist(reconst_result_2[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_2  <- matrix(unlist(reconst_result_2[['U_reconst_list']]), ncol=3) 
#' ##
#' par(mfrow=c(1,3))
#' matplot(x=U_mat[,1:3], y=Y_mat[,1:3], ylab="", col=gray(.5), type="l", 
#' main="Orig. Data", xlim=c(a,b))
#' matplot(x=U_reconst_mat_1, y=Y_reconst_mat_1, col=gray(.5), 
#' type="l", main="With Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2) 
#' matplot(x=U_reconst_mat_2, y=Y_reconst_mat_2, col=gray(.5), 
#' type="l", main="Without Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2)
#' par(mfrow=c(1,1))
#' 
#' 
#' ## Generate partially observed functional data without error 
#' SimDat        <- simuldata(n = 50, a = a, b = b, DGP="DGP2")
#' ## 
#' Y_list   <- SimDat[['Y_list']]; Y_mat <- SimDat[['Y_mat']]
#' U_list   <- SimDat[['U_list']]; U_mat <- SimDat[['U_mat']]
#' ##
#' ## Reconstruction with alignments of reconstructed parts and observed fragments
#' reconst_result_1 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error=0_AlignYES_CommonGrid', reconst_fcts = 1:3)
#' Y_reconst_mat_1  <- matrix(unlist(reconst_result_1[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_1  <- matrix(unlist(reconst_result_1[['U_reconst_list']]), ncol=3) 
#' ##
#' ## Reconstruction without alignments
#' reconst_result_2 <- reconstructKneipLiebl(Ly = Y_list, Lu = U_list, 
#' method = 'Error>=0_AlignNO', reconst_fcts = 1:3)
#' Y_reconst_mat_2  <- matrix(unlist(reconst_result_2[['Y_reconst_list']]), ncol=3) 
#' U_reconst_mat_2  <- matrix(unlist(reconst_result_2[['U_reconst_list']]), ncol=3) 
#' ##
#' par(mfrow=c(1,3))
#' matplot(x=U_mat[,1:3], y=Y_mat[,1:3], ylab="", col=gray(.5), type="l", 
#' main="Orig. Data", xlim=c(a,b))
#' matplot(x=U_reconst_mat_1, y=Y_reconst_mat_1, col=gray(.5), 
#' type="l", main="With Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2) 
#' matplot(x=U_reconst_mat_2, y=Y_reconst_mat_2, col=gray(.5), 
#' type="l", main="Without Alignment", ylab="", xlab="", xlim=c(a,b))
#' matlines(x=U_mat[,1:3], y=Y_mat[,1:3], col=gray(.2), lwd=2)
#' par(mfrow=c(1,1))

my_reconstructKneipLiebl <- function(Ly,
                                     Lu, 
                                     reconst_fcts = NULL,
                                     method       = c('Error=0_AlignYES_CommonGrid',
                                                      'Error>0_AlignYES',
                                                      'Error>=0_AlignNO',
                                                      'Error>0_AlignYES_CEscores',
                                                      'Error>0_AlignNO_CEscores',
                                                      'PACE',
                                                      'Error=0_PACE'),
                                     K            = NULL,
                                     nRegGrid     = NULL,
                                     maxbins      = NULL, 
                                     progrbar     = FALSE){
  
  ## Useful methods ------------------------------------------------
  source('KL/my_fpca.R')
  source('KL/my_gcvKneipLiebl.R')
  source('KL/my_reconstKneipLiebl_fun.R')
  source('KL/my_PACE_fun.R')
  
  ##
  method <- switch(method, 
                   'Error=0_AlignYES_CommonGrid' = 1, 
                   'Error>0_AlignYES'            = 2,
                   'Error>=0_AlignNO'            = 3,
                   'Error>0_AlignYES_CEscores'   = 4,
                   'Error>0_AlignNO_CEscores'    = 5,
                   'PACE'                        = 6,
                   'Error=0_PACE'                = 7)
  ##
  center <- TRUE

  ##
  Y_reconst_list <- vector("list", length(reconst_fcts))
  U_reconst_list <- vector("list", length(reconst_fcts))
  K_vec          <- rep(NA,        length(reconst_fcts))
  ##
  if(method == 1){# 'Error=0_AlignYES_CommonGrid'
    ## method 1: 
    ## Requires: error=0, common argvalues (i.e., data structure as in Kraus JRSSB)
    ## Uses:     Classical (integral) scores. Alignment of the reconstructed parts and the fully observed fragments.
    fpca_obj <- my_fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if( any(fpca_obj$obs_argvalsO[[i]] != fpca_obj$argvalsO[[i]]) ){stop("The fragment must be fully observed (obs_argvalsO == argvalsO).") }
      ##
      if(is.null(K)){
        K_vec[i]   <- my_gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 1,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      fragmO     <- c(stats::na.omit(c(fpca_obj$Y[i,])))
      ##
      result_tmp <- my_reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$scoresO[[i]], 
                                          efun_reconst= fpca_obj$efun_reconst[[i]],
                                          fragmO      = fragmO, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
    }
  }
  if(method == 2){# 'Error>0_AlignYES'
    ## method 2: 
    ## Requires: error>0, random or common argvalues (needs a large number of observed argvalues per function)
    ## Does:     Classical (integral) scores. Pre-smoothing of observed fragments. Alignment of reconstructed parts and pre-smoothed fragments.
    fpca_obj <- my_fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if( !all(range(fpca_obj$obs_argvalsO[[i]]) == range(fpca_obj$argvalsO[[i]])) ){
        stop("The range of obs_argvalsO of the fragment must equal the range of argvalsO.") }
      ##
      if(is.null(K)){
        K_vec[i]   <- my_gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 2,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      smooth.fit        <- suppressMessages(stats::smooth.spline(y=c(stats::na.omit(c(fpca_obj$Y[i,]))), x=fpca_obj$obs_argvalsO[[i]]))
      fragmO_presmooth  <- stats::predict(smooth.fit, fpca_obj$argvalsO[[i]])$y
      ##
      result_tmp <- my_reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$scoresO[[i]], 
                                          efun_reconst= fpca_obj$efun_reconst[[i]],
                                          fragmO      = fragmO_presmooth, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
    }
  }
  if(method == 3){# 'Error>=0_AlignNO'
    ## method 3: 
    ## Requires: error>=0, random or common argvalues
    ## Uses:     Classical (integral) scores. No alignment 
    fpca_obj <- my_fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = FALSE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      ##
      if(is.null(K)){
        K_vec[i]   <- my_gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 3,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- my_reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                             argvals     = fpca_obj$argvals, 
                                             argvalsO    = fpca_obj$argvalsO[[i]], 
                                             scoresO     = fpca_obj$scoresO[[i]], 
                                             efun_reconst= fpca_obj$efun_reconst[[i]],
                                             fragmO      = NULL, 
                                             K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
    }
  }
  if(method == 4){# 'Error>0_AlignYES_CEscores'
    ## method 4: 
    ## Requires: error>0, random or common argvalues (needs a large number of observed argvalues per function for the pre-smoothing)
    ## Uses:     CEscores (PACE). Pre-smoothing of observed fragments. Alignment of reconstructed parts with the pre-smoothed fragments.
    fpca_obj <- my_fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = TRUE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if( !all(range(fpca_obj$obs_argvalsO[[i]]) == range(fpca_obj$argvalsO[[i]])) ){
        stop("The range of obs_argvalsO of the fragment must equal the range of argvalsO.") }
      ##
      if(is.null(K)){
        K_vec[i]   <- my_gcvKneipLiebl(fpca_obj = fpca_obj, 
                                       argvalsO = fpca_obj$argvalsO[[i]], 
                                       method   = 4,
                                       progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      smooth.fit        <- suppressMessages(stats::smooth.spline(y=c(stats::na.omit(c(fpca_obj$Y[reconst_fcts[i],]))), x=fpca_obj$obs_argvalsO[[i]]))
      fragmO_presmooth  <- stats::predict(smooth.fit, fpca_obj$argvalsO[[i]])$y
      ##
      result_tmp <- my_reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                             argvals     = fpca_obj$argvals, 
                                             argvalsO    = fpca_obj$argvalsO[[i]], 
                                             scoresO     = fpca_obj$CE_scoresO[[i]], 
                                             efun_reconst= fpca_obj$efun_reconst[[i]],
                                             fragmO      = fragmO_presmooth, 
                                             K           = K_vec[i])
      ##
      Y_reconst_list[[i]]   <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]   <- result_tmp[['x_reconst']]
      
      #x11(width=8000, height=5000)
      #par(mfrow=c(1,2))
      #matplot(fpca_obj$argvalsO[[i]], Y_list[[i]], type='l', xlim=c(0,10))
      #matplot(fpca_obj$argvals, Y_reconst_list[[i]], type='l', xlim=c(0,10))
    }
  }
  if(method == 5){# 'Error>0_AlignNO_CEscores'
    ## method 5: 
    ## Requires: error>0, random or common argvalues (works also for a few observed argvalues per function, since no pre-smoothing)
    ## Uses:     CEscores (PACE). No alignment.
    fpca_obj <- my_fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        CEscores     = TRUE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      ##
      if(is.null(K)){
        K_vec[i]   <- my_gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 5,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- my_reconstKneipLiebl_fun(mu          = fpca_obj$mu, 
                                          argvals     = fpca_obj$argvals, 
                                          argvalsO    = fpca_obj$argvalsO[[i]], 
                                          scoresO     = fpca_obj$CE_scoresO[[i]], 
                                          efun_reconst= fpca_obj$efun_reconst[[i]],
                                          fragmO      = NULL, 
                                          K           = K_vec[i])
      ##
      Y_reconst_list[[i]]  <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]  <- result_tmp[['x_reconst']]
    }
  }
  if(method == 6){# 'PACE'
    fpca_obj <- my_fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        PACE         = TRUE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      ##
      if(is.null(K)){
        K_vec[i]   <- my_gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 6,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- my_PACE_fun(mu          = fpca_obj$mu, 
                             argvals     = fpca_obj$argvals, 
                             scoresP     = fpca_obj$scoresP[[i]], 
                             efunctionsP = fpca_obj$efunctionsP, 
                             K           = K_vec[i])
      ##
      Y_reconst_list[[i]]  <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]  <- result_tmp[['x_reconst']]
    }
  } 
  if(method == 7){# 'PACE_Error==0'
    fpca_obj <- my_fpca(Ly           = Ly, 
                        Lu           = Lu, 
                        reconst_fcts = reconst_fcts, 
                        PACE_E0      = TRUE, 
                        center       = center, 
                        maxbins      = maxbins)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      ##
      if(is.null(K)){
        K_vec[i]   <- my_gcvKneipLiebl(fpca_obj = fpca_obj, 
                                    argvalsO = fpca_obj$argvalsO[[i]], 
                                    method   = 7,
                                    progrbar = progrbar)
      }else{K_vec[i] <- K}
      ##
      result_tmp <- my_PACE_fun(mu          = fpca_obj$mu, 
                             argvals     = fpca_obj$argvals, 
                             scoresP     = fpca_obj$scoresP[[i]], 
                             efunctionsP = fpca_obj$efunctionsP, 
                             K           = K_vec[i])
      ##
      Y_reconst_list[[i]]  <- result_tmp[['y_reconst']]
      U_reconst_list[[i]]  <- result_tmp[['x_reconst']]
    }
  }
  ##
  if(!is.null(nRegGrid)){
    ## Evaluate the reconstruced functions at a regular gird of length nRegGrid
    xout <- seq(from = min(fpca_obj$argvals), to = max(fpca_obj$argvals), len=nRegGrid)
    ##
    for(i in 1:length(reconst_fcts)){ # i <- 1
      if(all(is.na(Y_reconst_list[[i]]))){
        Y_reconst_list[[i]]  <- rep(NA, length(xout))
        U_reconst_list[[i]]  <- xout
      } else {
        Reconstr_on_RegGrid  <- stats::spline(y = Y_reconst_list[[i]], x = U_reconst_list[[i]], xout = xout)
        Y_reconst_list[[i]]  <- Reconstr_on_RegGrid$y
        U_reconst_list[[i]]  <- Reconstr_on_RegGrid$x
      }
    }
  }
  ##
  return(list(
    "Y_reconst_list"  = Y_reconst_list,
    "U_reconst_list"  = U_reconst_list,
    "K"               = K_vec
  ))
}



