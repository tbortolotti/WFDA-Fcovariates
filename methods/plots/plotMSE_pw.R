plotMSE_pw <- function(MSE.vec, MSE.ita18 = rep(0,37), t.points, name_dir)
{
  
  library(wesanderson)
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
                  color = pal[1]), size=1) +
    scale_color_manual(name="Model:",values=c(pal[2],pal[1]), labels = c("Functional", "Scalar")) +
    scale_x_continuous(breaks=x.ticks, labels=10^x.ticks) +
    scale_y_continuous(breaks=seq(0.08,0.18, by=0.02)) +
    #ggtitle(paste0("MSE -") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(text = element_text(size = 12)) +
    #theme(axis.title.x = element_text(size = 16)) +
    theme(axis.text.x = element_text(size = 12)) +
    #theme(axis.title.y = element_text(size = 16)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    xlab("Period [s]") + ylab("MSE")
  
  ggsave(filename = paste0("MSE_pointwise.png"),
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