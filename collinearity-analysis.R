setwd('/Users/teresabortolotti/Documents/R/WFDA-Fcovariates')

library(fda)
library(fdakma)
library(roahd)
library(coda)
library(devtools)
library(fastmatrix)
library(R.matlab)
library(latex2exp)
library(calculus)
library(ReconstPoFD)
library(tidyverse)
library("xtable")
library(snowfall)
library(psych)
library(progress)
library(wesanderson)
library(ggplot2)
library(ggcorrplot)
library(reshape2)

rm(list=ls())
graphics.off()
cat("\014")

## Load Data -------------------------------------------------------------------
load('DATA/curves.RData')
load('DATA/t_period.RData')
load('DATA/obs.RData')
load('DATA/T_hp.RData')
load('DATA/xlist-logg.RData')
load('DATA/data.RData')
load('DATA/events.RData')

t.points <- log10(T.period)
t.points[1] <- -2.5

## Analysis of multicollinearity in the covariates -----------------------------
x7 <- xlist[[7]]
x8 <- xlist[[8]]
cor.mat.78 <- cor.fd(t.points, x7, t.points, x8)

x6 <- xlist[[6]]
x3 <- xlist[[3]]
cor.mat.36 <- cor.fd(t.points, x3, t.points, x6)

x2 <- xlist[[2]]
cor.mat.23 <- cor.fd(t.points, x2, t.points, x3)

cor.mat.26 <- cor.fd(t.points, x2, t.points, x6)

cor.mat.37 <- cor.fd(t.points, x3, t.points, x7)

xticks <- seq(-2,1,1)
ext.ticks <- c(-2.5,xticks)

levs <- 5

pdf(file = "/Users/teresabortolotti/Desktop/flash/correlations/geom-and-anel.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.78,
        xlab="Geometric Attenuation",
        ylab="Anelastic Attenuation",
        #color.palette = function(n) hcl.colors(n, "RdYlBu", rev = TRUE),
        main=paste("(a) Correlation across periods for\n",
                   "Geometric and Anelastic Attenuation"),
        cex.main=0.8, axes=FALSE, nlevels = levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
dev.off()

pdf(file = "/Users/teresabortolotti/Desktop/flash/correlations/highmag-and-geom.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.36,
        xlab="High Magnitudes",
        ylab="Geometric Attenuation",
        main=paste("(d) Correlation across periods for\n",
                   "High Magnitudes and Geometric Attenuation"),
        cex.main=0.8, axes=FALSE, nlevels=levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks)
dev.off()

pdf(file = "/Users/teresabortolotti/Desktop/flash/correlations/lowmag-and-geom.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.26,
        xlab="Low Magnitudes",
        ylab="Geometric Attenuation",
        main=paste("(c) Correlation across periods for\n",
                   "Low Magnitudes and Geometric Attenuation"),
        cex.main=0.8, axes=FALSE, nlevels=levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks)
dev.off()

pdf(file = "/Users/teresabortolotti/Desktop/flash/correlations/lowmag-and-highmag.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.23,
        xlab="Low Magnitudes",
        ylab="High Magnitudes",
        main=paste("(b) Correlation across periods for\n",
                   "Low and High Magnitudes"),
        cex.main=0.8, axes=FALSE, nlevels=levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks)
dev.off()

pdf(file = "/Users/teresabortolotti/Desktop/flash/correlations/highmag-and-anel.pdf", width = 4, height = 4)
contour(t.points, t.points, cor.mat.37,
        xlab="High Magnitudes",
        ylab="Anelastic Attenuation",
        #color.palette = function(n) hcl.colors(n, "RdYlBu", rev = TRUE),
        main=paste("(e) Correlation across periods for\n",
                   "High Magnitudes and Anelastic Attenuation"),
        cex.main=0.8, axes=FALSE, nlevels = levs,
        lwd=2, labcex=0.8, vfont=c("sans serif", "bold italic"))
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
dev.off()

## FILLED CONTOUR 
levs <- 5
palette <- wes_palette("Zissou1", 15, type = "continuous")
filled.contour(t.points, t.points, cor.mat,
               xlab="Period",
               ylab="Period",
               zlim=c(-1,1),
               color.palette = function(n) hcl.colors(n, "RdYlBu", rev = TRUE),
               #color.palette = wes_palette("Zissou1", 10, type = "continuous"),
               main=paste("Correlation across periods for\n",
                          "geometric and anelastic attenuation"),
               cex.main=1, axes=FALSE, nlevels = levs,
               lwd=2)
axisIntervals(1, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )
axisIntervals(2, atTick1=ext.ticks, atTick2=NA,
              atLabels=xticks,
              labels=10^xticks )



## CORRPLOT

pal <- wes_palette("Zissou1")
ggcorrplot(cor.mat, lab=TRUE,
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c(pal[1], "white", pal[5]))

ggsave(filename = paste0("corrmatrix.png"),
       plot = last_plot(),
       device = NULL,
       path = new.dir,
       scale = 1,
       limitsize = FALSE,
       dpi = 320)

dev.off()

