my_irreg2mat <- function(ydata, binning = FALSE, maxbins = 1000){
  ##
  ydata <- ydata[stats::complete.cases(ydata), ]
  nobs  <- length(unique(ydata$.id))
  newid <- as.numeric(as.factor(ydata$.id))
  bins  <- sort(unique(ydata$.index))
  if(binning && (length(bins) > maxbins)){
    binvalues <- seq((1 - 0.001 * sign(bins[1])) * bins[1], 
                     (1 + 0.001 * sign(bins[length(bins)])) * bins[length(bins)], 
                     l = maxbins + 1)
    bins      <- binvalues
    binvalues <- utils::head(stats::filter(binvalues, c(0.5, 0.5)), -1)
  }else{
    binvalues <- bins
    bins      <- c((1 - 0.001 * sign(bins[1])) * bins[1], bins[-length(bins)], 
                   (1 + 0.001 * sign(bins[length(bins)])) * bins[length(bins)])
    if(bins[1] == 0           ){bins[1]            <- -0.001}
    if(bins[length(bins)] == 0){bins[length(bins)] <-  0.001}
  }
  newindex         <- cut(ydata$.index, breaks = bins, include.lowest = TRUE)
  Y                <- matrix(NA, nrow = nobs, ncol = nlevels(newindex))
  colnames(Y)      <- binvalues
  attr(Y, "index") <- binvalues
  ## If there are more than one data-point within a bin, 
  ## then only one of these is used (the last one).
  Y[cbind(newid, as.numeric(newindex))] <- ydata$.value
  ##
  return(Y)
}
