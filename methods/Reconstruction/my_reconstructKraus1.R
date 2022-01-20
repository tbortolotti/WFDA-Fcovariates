#' Reconstruct principal component scores as proposed by Kraus (JRSSB, 2015).
#'
#' This function allows to reconstruct the principal component scores of partially
#' observed functional data as proposed in: 
#' Kraus, D. (2015). Components and completion of partially observed functional data. 
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology), 77(4), 777-801. 
#' @param X_mat        pxn matrix (p: number of discretization points, n=number of functions
#' @param alpha        Ridge parameter. If alpha = NULL (default), an optimal alpha is determined by GCV
#' @param reconst_fcts A vector specifying the list elements in Ly which need to be reconstructed.
#'                     Default (reconst_fcts=NULL) will reconstruct all functions.
#' @export reconstructKraus
#' @examples  
#' SimDat       <- simuldata(n = 50, a = 0, b = 1, DGP="DGP3")
#' Y_mat        <- SimDat[['Y_mat']]
#' U_mat        <- SimDat[['U_mat']]
#' U_true_mat   <- SimDat[['U_true_mat']]
#' ##
#' result        <- reconstructKraus(X_mat = Y_mat)
#' Y_reconst_mat <- result[['X_reconst_mat']]
#' ##
#' par(mfrow=c(2,1))
#' matplot(x=U_mat[,1:5], y=Y_mat[,1:5], col=gray(.5), type="l", 
#' main="Original Data", ylab="", xlab="")
#' matplot(x=U_true_mat[,1:5], y=Y_reconst_mat[,1:5], col=gray(.5), 
#' type="l", main="Kraus (2015)", ylab="", xlab="")
#' par(mfrow=c(1,1))
my_reconstructKraus1 <- function(X_mat, tfine, alpha = NULL, reconst_fcts = NULL)
{
  
  ## Immagino che X_mat sia una matrice di funzioni valutate su una griglia
  ## denza di punti ben equispaziati
  
  ## Load useful functions -----------------------------------------
  source('methods/Reconstruction/Kraus1/meanKraus.R')
  source('methods/Reconstruction/Kraus1/covKraus.R')
  source('methods/Reconstruction/Kraus1/reconstKraus1_fun.R')
  source('methods/Reconstruction/Kraus1/gcvKraus.R')
  gcvKraus <- Vectorize(FUN = gcvKraus, vectorize.args = "alpha")
  
  ##
  mean_vec      <- meanKraus(X_mat)
  cov_mat       <- covKraus(X_mat)
  n             <- ncol(X_mat)

  if(is.null(reconst_fcts)){
    reconst_fcts <- 1:n
  }
  X_reconst_mat <- X_mat[, reconst_fcts, drop=FALSE]

  NonNA_fcts    <- apply(X_mat,2,function(x)!any(is.na(x)))
  X_Compl_mat   <- X_mat[,NonNA_fcts]
  alpha_vec     <- rep(NA, length(reconst_fcts)) 
  df_vec        <- rep(NA, length(reconst_fcts)) 

  for(i in 1:length(reconst_fcts)){
    X_tmp      <- X_mat[,reconst_fcts[i]]
    ##
    M_bool_vec <- is.na(X_tmp)
    O_bool_vec <- !M_bool_vec
    ##
    if(is.null(alpha)){
      alpha_vec[i] <- optimize(f = function(alpha){gcvKraus(cov_mat     = cov_mat, 
                                                            mean_vec    = mean_vec, 
                                                            X_Compl_mat = X_Compl_mat, 
                                                            M_bool_vec  = M_bool_vec,
                                                            tfine       = tfine,
                                                            alpha       = alpha)},
                               interval = c(.Machine$double.eps, sum(diag(cov_mat))*n), maximum = FALSE)$minimum
    }else{
      alpha_vec[i] <- alpha
    }
    ##
    result_tmp <- reconstKraus1_fun(cov_mat   = cov_mat, 
                                   X_cent_vec = c(X_tmp - mean_vec), 
                                   alpha      = alpha_vec[i],
                                   tfine      = tfine)
    ##
    X_reconst_mat[,i]  <- c(result_tmp[['X_cent_reconst_vec']] + mean_vec)
    df_vec[i]          <- result_tmp[['df']]
  }
  
  X_compl_mat                <- X_mat
  X_compl_mat[,reconst_fcts] <- X_reconst_mat
  
  return(list("X_compl_mat"   = X_compl_mat,
              "alpha"         = alpha_vec, 
              "df"            = df_vec
  ))
}
