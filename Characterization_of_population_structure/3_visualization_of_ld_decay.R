# Script for visualizing LD dacay for LSP individuals to choose an arbitrary value for LD prunning

rm(list = ls())
setwd("/home/mullerta/work/SP_WP/Analyses/LD_prunning/ld_decay/")
library(tidyverse)

# set path
my_bins_chr1 <- "ld_decay_thin0.1_chrA.ld_decay_bins"
my_bins_chrZ <- "ld_decay_thin0.1_chrZ.ld_decay_bins"


# read in data
ld_bins_chr1 <- read_tsv(my_bins_chr1)
ld_bins_chrZ <- read_tsv(my_bins_chrZ)

ld<-rbind(ld_bins_chr1,ld_bins_chrZ)
ld$chr <- c(rep("chr1",nrow(ld_bins_chr1)),rep("chrZ",nrow(ld_bins_chrZ)))

# plot LD decay (Figure S4)
ggplot(data = ld, aes(distance, avg_R2, color=chr))+
  geom_line()+
  scale_color_manual(values = c("turquoise","red"))+
  xlab("Distance (in bp)") + ylab(expression(italic(r)^2))+xlim(0,100000)+ylim(0,0.80)+theme_minimal(
  )+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(size = 18),  
    axis.title.y = element_text(size = 18, margin = margin(r = 40)),  
    axis.text.x = element_text(size = 18),  
    axis.title.x = element_text(size = 18),  
    axis.ticks = element_line(size = 0.5),  
    axis.line = element_line(size = 0.5, color = "black"), 
    plot.title = element_text(hjust = 0.5, size = 17, face = "bold"),  
    legend.text = element_text(size = 18),  
    legend.title = element_text(size = 18)) 

