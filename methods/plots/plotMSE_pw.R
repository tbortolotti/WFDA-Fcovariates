plotMSE_pw <- function(MSE.vec, MSE.ita18 = rep(0,37), t.points, name_dir)
{
  
  library(wesanderson)
  library(latex2exp)
  pal  <- wes_palette("Cavalcanti1")
  
  # create the directory to save plots
  dir.current <- getwd()
  new.dir <- paste0(dir.current,"/Results/",name_dir)
  dir.create(new.dir)
  
  df <- data.frame(t.points, MSE.vec, MSE.ita18)
  x.ticks <- seq(-3,1,by=1)

  t.plot <- t.points
  
  ggplot(df) +
    geom_line(aes(x = t.plot, y = MSE.vec,
                  color = pal[2]), size=1) +
    geom_line(aes(x = t.plot, y = MSE.ita18,
                  color = "darkorange"), size=1) +
    scale_color_manual(name="Model:",values=c(pal[2],"darkorange"), labels = c("F-ITA18", "ITA18")) +
    scale_x_continuous(breaks=x.ticks, labels=10^x.ticks) +
    scale_y_continuous(breaks=seq(0.08,0.18, by=0.02)) +
    theme_bw() +
    labs(x="Period [s]", y="MSE", title='(a) Comparison of point-wise MSE') +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size=22),
          text = element_text(size = 22),
          axis.text.x = element_text(size = 22),
          axis.text.y = element_text(size = 22),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 19))
  
  ggsave(filename = paste0("MSE_pointwise.pdf"),
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