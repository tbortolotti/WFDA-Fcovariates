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
                  color = "darkorange"), size=1) +
    scale_color_manual(name="Model:",values=c(pal[2],"darkorange"), labels = c("F-ITA18", "ITA18")) +
    scale_x_continuous(breaks=x.ticks, labels=10^x.ticks) +
    scale_y_continuous(breaks=seq(0.25,0.35, by=0.05)) +
    theme_bw() +
    labs(x="Period [s]", y=TeX(r'($\hat{\sigma}$)'), title='(b) Comparison of residuals standard deviation') +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
          text = element_text(size = 22),
          axis.text.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 19))
  
  ggsave(filename = paste0("sigma_comparison.pdf"),
         plot = last_plot(),
         width = 9,
         height = 5,
         units = "in",
         device = NULL,
         path = new.dir,
         scale = 1,
         limitsize = TRUE,
         dpi = 300)
  
}



