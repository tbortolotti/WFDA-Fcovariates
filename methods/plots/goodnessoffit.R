goodnessoffit <- function(mod.fit, input, data, name_dir)
{
  library(ggplot2)
  library(ggpubr)
  
  source('methods/Regression/my_predict_fRegress.R')
  
  # create the directory to save plots
  dir.current <- getwd()
  new.dir <- paste0(dir.current,"/Results/",name_dir)
  dir.create(new.dir)
  
  library(wesanderson)
  pal <- wes_palette('Cavalcanti1')
  # Data ----------------------------------------------------------------------
  
  X <- input$xlist
  names(X) <- c('Intercept', 'Low Magnitude', 'High Magnitude', 'SS', 'TF', 'Geometrical spreading (M)', 'Geometrical spreading', 'Anelastic spreading', 'VS30')
  blist <- input$blist
  t.points <- input$t.points
  basis <- mod.fit$yhatfdobj$basis
  theta.t <- eval.basis(t.points, basis)
  theta <- t(theta.t)
  n <- dim(X[[1]]$coefs)[2]
  q <- length(X)
  
  # Fitted values  --------------------------------------------
  y.hat     <- mod.fit$yhatfdobj
  y         <- mod.fit$yfdobj
  rangeval  <- y$basis$rangeval
  res       <- y - y.hat
  y.vals    <- theta.t %*% y$coefs
  yhat.vals <- theta.t %*% y.hat$coefs
  
  a0.est <- mod.fit$betaestlist[[1]]$fd
  b1.est <- mod.fit$betaestlist[[2]]$fd
  b2.est <- mod.fit$betaestlist[[3]]$fd
  f1.est <- mod.fit$betaestlist[[4]]$fd
  f2.est <- mod.fit$betaestlist[[5]]$fd
  c1.est <- mod.fit$betaestlist[[6]]$fd
  c2.est <- mod.fit$betaestlist[[7]]$fd
  c3.est <- mod.fit$betaestlist[[8]]$fd
  k0.est <- mod.fit$betaestlist[[9]]$fd
  
  # Evaluation of important quantities --------------------------------
  ## for each residual, evaluate the integral. Generate a vector to
  ## save the values of these integrals
  onesbasis <- create.constant.basis(rangeval)
  onesfd    <- fd(1,onesbasis)
  
  res2 <- res*res
  
  res.integrals  <- inprod(res, onesfd, 0, 0, rangeval)
  yfit.integrals <- inprod(y.hat, onesfd, 0, 0, rangeval)
  res.norms <- inprod(res2, onesfd, 0, 0, rangeval)

  # Boxplots residuals vs regressors ----------------------------------
  
  box.dir <- paste0(dir.current,"/Results/",name_dir,"/boxplot_residuals")
  dir.create(box.dir)
  
  #reg.idx <- c(2,3,6,7,8,9)
  names(data) <- c('Joyner-Boore distance', 'Magnitude', 'Style-of-faulting', 'Shear-wave velocity')
  
  for(j in 1:length(data))
  {
    if(!(j==3))
    {
      
      #jdx <- reg.idx[j]
      reg <- data[[j]]
      xx <- as.numeric(quantile(reg, c(0,0.20,0.40,0.60,0.80,1)))
      my.list <- list()
      for(i in 1:(length(xx)-1))
      {
        x.break_1 <- xx[i]
        x.break_2 <- xx[i+1]
        idxs <- which(reg>x.break_1 & reg<=x.break_2)
        integr <- numeric(length(idxs))
        for(kk in 1:length(idxs))
        {
          k <- idxs[kk]
          integr[kk] <- res.integrals[k]/10
        }
        my.list[[i]] <- integr
      }
      
      x.val <- round((xx[2] + xx[1])/2,digits=2)
      df <- cbind(my.list[[1]], rep(x.val, time=length(my.list[[1]])))
      for(i in 2:length(my.list))
      {
        x.val <- round((xx[i+1] + xx[i])/2, digits=2)
        df.add <- cbind(my.list[[i]], rep(x.val, time=length(my.list[[i]])))
        df <- rbind(df,df.add)
      }
      df <- as.data.frame(df)
      names(df) <- c('y','x')
      
      ggplot(df, aes(x, y, group=x, fill=as.factor(x))) + 
        geom_boxplot(alpha=0.4) +
        scale_fill_brewer(palette="YlGn") +
        geom_hline(yintercept=0, color=pal[1], size=1, linetype=5) +
        theme_bw() +
        theme(legend.position = "none") + 
        #ggtitle(names(data)[j]) +
        #theme(plot.title = element_text(hjust = 0.5)) +
        xlab(names(data)[j]) + ylab("Residuals integrals")
      
      ggsave(filename = paste0(j,".png"),
             plot = last_plot(),
             width = 8,
             height = 5,
             units = "in",
             device = NULL,
             path = box.dir,
             scale = 1,
             limitsize = TRUE,
             dpi = 300)
      
      #dev.off()
      
    }
    
  }
  
  ## Points integral means vs regressors ------------------------------------
  box.dir <- paste0(dir.current,"/Results/",name_dir,"/dots_residuals")
  dir.create(box.dir)
  
  #reg.idx <- c(2,3,6,7,8,9)
  names(data) <- c('Joyner-Boore distance', 'Magnitude', 'Style-of-faulting', 'Shear-wave velocity')
  
  for(j in 1:length(data))
  {
    if(!(j==3))
    {
      
      #jdx <- reg.idx[j]
      reg <- data[[j]]
      res.int <- res.integrals/10
      
      png(file = paste0(box.dir,"/",j,".png"), width = 8000, height = 5000, units = "px", res=800)
      par(mar=c(5,6,4,1)+.1)
      plot(reg, res.int, pch=1, col='grey70', xlab=names(data)[j], ylab="Residuals integrals", cex.lab=2.2, cex.axis=2.2)
      abline(h=0, type='l', col=pal[2], lty=2, lwd=3)
      dev.off()
    }
    
  }

  
  ## Norm of residuals vs regressors ---------------------------------
  
  box.dir <- paste0(dir.current,"/Results/",name_dir,"/boxplot_normresiduals")
  dir.create(box.dir)
  
  names(data) <- c('Joyner-Boore distance', 'Magnitude', 'Style-of-faulting', 'Shear-wave velocity')
  
  for(j in 1:length(data))
  {
    if(!(j==3))
    {
      #jdx <- reg.idx[j]
      reg <- data[[j]]
      xx <- as.numeric(quantile(reg, c(0,0.20,0.40,0.60,0.80,1)))
      my.list <- list()
      for(i in 1:(length(xx)-1))
      {
        x.break_1 <- xx[i]
        x.break_2 <- xx[i+1]
        idxs <- which(reg>x.break_1 & reg<=x.break_2)
        integr <- numeric(length(idxs))
        for(kk in 1:length(idxs))
        {
          k <- idxs[kk]
          integr[kk] <- sqrt(res.norms[k])
        }
        my.list[[i]] <- integr
      }
      
      x.val <- round((xx[2] + xx[1])/2,digits=2)
      df <- cbind(my.list[[1]], rep(x.val, time=length(my.list[[1]])))
      for(i in 2:length(my.list))
      {
        x.val <- round((xx[i+1] + xx[i])/2, digits=2)
        df.add <- cbind(my.list[[i]], rep(x.val, time=length(my.list[[i]])))
        df <- rbind(df,df.add)
      }
      df <- as.data.frame(df)
      names(df) <- c('y','x')
      
      ggplot(df, aes(x, y, group=x, fill=as.factor(x))) + 
        geom_boxplot(alpha=0.4) +
        scale_fill_brewer(palette="YlGn") +
        geom_hline(yintercept=0, color=pal[1], size=1, linetype=5) +
        theme_bw() +
        theme(legend.position = "none") + 
        #ggtitle(names(data)[j]) +
        #theme(plot.title = element_text(hjust = 0.5)) +
        xlab(names(data)[j]) + ylab("")
      
      ggsave(filename = paste0(j,".png"),
             plot = last_plot(),
             device = NULL,
             width = 8,
             height = 5,
             units = "in",
             path = box.dir,
             scale = 1,
             limitsize = TRUE,
             dpi = 300)
      
    }
    
  }
  
  
  # ## Linearity check ------------------------------------------------------
  # # Per ogni regressore, voglio togliere alle y l'effetto del regressore e plottare
  # # il residuo vs la x, per vedere di non aver lasciato gi? della dipendenza non
  # # lineare
  # 
  # lin.dir <- paste0(dir.current,"/Results/",name_dir,"/linearity_check")
  # dir.create(lin.dir)
  # 
  # reg.idx <- c(2,3,6,7,8,9)
  # 
  # for(j in 1:length(reg.idx))
  # {
  #   
  #   X.reg <- matrix(data=0, nrow=n, ncol=q)
  #   
  #   jdx <- reg.idx[j]
  #   reg <- X[,jdx]
  #   
  #   X.reg[,jdx] <- X[,jdx]
  #   xlist.reg <- lapply(seq_len(ncol(X.reg)), function(i) X.reg[,i])
  #   
  #   y.reg     <- my_predict_fRegress(mod          = mod.fit,
  #                                    xlist        = xlist.reg,
  #                                    t.points     = t.points)
  #   
  #   res.reg <- y - y.reg
  #   res.reg.int  <- inprod(res.reg, onesfd, 0, 0, rangeval)
  #   
  #   
  #   xx <- as.numeric(quantile(reg, c(0,0.20,0.40,0.60,0.80,1)))
  #   my.list <- list()
  #   for(i in 1:(length(xx)-1))
  #   {
  #     x.break_1 <- xx[i]
  #     x.break_2 <- xx[i+1]
  #     idxs <- which(reg>x.break_1 & reg<=x.break_2)
  #     integr <- numeric(length(idxs))
  #     for(kk in 1:length(idxs))
  #     {
  #       k <- idxs[kk]
  #       integr[kk] <- res.reg.int[k]/10
  #     }
  #     my.list[[i]] <- integr # lista che salva gli integrali relativi al batch i
  #   }
  #   
  #   x.val <- round((xx[2] + xx[1])/2,digits=2)
  #   df <- cbind(my.list[[1]], rep(x.val, time=length(my.list[[1]])))
  #   for(i in 2:length(my.list))
  #   {
  #     x.val <- round((xx[i+1] + xx[i])/2, digits=2)
  #     df.add <- cbind(my.list[[i]], rep(x.val, time=length(my.list[[i]])))
  #     df <- rbind(df,df.add)
  #   }
  #   df <- as.data.frame(df)
  #   names(df) <- c('y','x')
  #   
  #   ggplot(df, aes(x, y, group=x, fill=as.factor(x))) + 
  #     geom_boxplot(alpha=0.4) +
  #     scale_fill_brewer(palette="YlGn") +
  #     geom_hline(yintercept=0, color=pal[1], size=1, linetype=5) +
  #     theme_bw() +
  #     theme(legend.position = "none") + 
  #     ggtitle(paste0(names(X)[jdx])) +
  #     theme(plot.title = element_text(hjust = 0.5))
  #   
  #   ggsave(filename = paste0(jdx,".png"),
  #          plot = last_plot(),
  #          width = 8,
  #          height = 5,
  #          units = "in",
  #          device = NULL,
  #          path = lin.dir,
  #          scale = 1,
  #          limitsize = TRUE,
  #          dpi = 300)
  # }

  # Residuals integrals vs fitted values integrals ------------------------
  
  reg <- yfit.integrals
  xx <- as.numeric(quantile(reg, c(0,0.20,0.40,0.60,0.80,1)))
  my.list <- list()
  for(i in 1:(length(xx)-1))
  {
    x.break_1 <- xx[i]
    x.break_2 <- xx[i+1]
    idxs <- which(reg>x.break_1 & reg<=x.break_2)
    integr <- numeric(length(idxs))
    for(kk in 1:length(idxs))
    {
      k <- idxs[kk]
      integr[kk] <- res.integrals[k]/10
    }
    my.list[[i]] <- integr
  }
  
  x.val <- round((xx[2] + xx[1])/2,digits=2)
  df <- cbind(my.list[[1]], rep(x.val, time=length(my.list[[1]])))
  for(i in 2:length(my.list))
  {
    x.val <- round((xx[i+1] + xx[i])/2, digits=2)
    df.add <- cbind(my.list[[i]], rep(x.val, time=length(my.list[[i]])))
    df <- rbind(df,df.add)
  }
  df <- as.data.frame(df)
  names(df) <- c('y','x')
  
  ggplot(df, aes(x, y, group=x, fill=as.factor(x))) + 
    geom_boxplot(alpha=0.4) +
    scale_fill_brewer(palette="YlGn") +
    geom_hline(yintercept=0, color=pal[1], size=1, linetype=5) +
    theme_bw() +
    theme(legend.position = "none") +
    #ggtitle("Residual integrals vs. Fitted values integrals") +
    #theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Residual integrals") + xlab("Fitted values integrals")
  ggsave(filename = paste0("resintegrals_vs_fittedintegrals.png"),
         plot = last_plot(),
         width = 8,
         height = 5,
         units = "in",
         device = NULL,
         path = new.dir,
         scale = 1,
         limitsize = TRUE,
         dpi = 300)
  
  
  ## Points residuals integrals vs fitted integrals --------------------------
  box.dir <- paste0(dir.current,"/Results/",name_dir,"/dots_residuals")
  
  yfit.int <- yfit.integrals/10
  res.int <- res.integrals/10
  
  png(file = paste0(box.dir,"/resvsfittedvals.png"), width = 8000, height = 5000, units = "px", res=800)
  par(mar=c(5,6,4,1)+.1)
  plot(yfit.int, res.int, pch=1, col='grey70', xlab="Fitted-values integrals", ylab="Residuals integrals", cex.lab=2.2, cex.axis=2.2)
  abline(h=0, type='l', col=pal[2], lty=2, lwd=3)
  dev.off()
  
}