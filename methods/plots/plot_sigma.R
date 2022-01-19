plot_sigma <- function(SigmaE, sigma.ITA18, t.points, name_dir)
{
  library(wesanderson)
  pal <- wes_palette('Cavalcanti1')
  
  cur.dir <- getwd()
  new.dir <- paste0(cur.dir,'/Results/',name_dir)
  dir.create(new.dir)
  
  sigma.t <- sqrt(diag(SigmaE))
  
  df <- data.frame(t.points, sigma.ITA18, sigma.t)
  x.ticks <- seq(-3,1,by=1)
  
  t.plot <- t.points
  
  ggplot(df) +
    geom_line(aes(x = t.plot, y = sigma.t,
                  color = pal[2]), size=1) +
    geom_line(aes(x = t.plot, y = sigma.ITA18,
                  color = pal[1]), size=1) +
    scale_color_manual(name="Model:",values=c(pal[2],pal[1]), labels = c("Functional", "Scalar")) +
    scale_x_continuous(breaks=x.ticks, labels=10^x.ticks) +
    scale_y_continuous(breaks=seq(0.25,0.35, by=0.05)) +
    #ggtitle(paste0("MSE -") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(text = element_text(size = 20)) +
    #theme(axis.title.x = element_text(size = 16)) +
    theme(axis.text.x = element_text(size = 20)) +
    #theme(axis.title.y = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 20)) +
    theme(legend.text = element_text(size = 20)) +
    xlab("Period [s]") + ylab(TeX('$\\hat{\\sigma}$'))
  
  ggsave(filename = paste0("sigma_comparison.png"),
         plot = last_plot(),
         width = 8,
         height = 3,
         units = "in",
         device = NULL,
         path = new.dir,
         scale = 1,
         limitsize = TRUE,
         dpi = 300)
  
}



