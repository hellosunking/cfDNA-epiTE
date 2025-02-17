rm(list=ls())
#load packages
library(pROC)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggprism)
library(ggbeeswarm)

args <- commandArgs(trailing=T)
aa<- args[1]
bb<-args[2]

size <- read.table(aa,header=T)

pdf(paste0( bb ,".Chromatin.State.pdf"), width=12, height=3 )
  df2 <- melt(size)
  p <- ggplot(df2, aes(x = variable, y = value,colour=variable)) +
    geom_boxplot(color = "gray")+theme_bw()+
    geom_beeswarm(method = "swarm", corral = "random", corral.width = 0.9) +
    labs(title = bb, x=NULL, y="Density")+
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "#000000",linewidth = 0.3),
          axis.text = element_text(colour = "#000000"),
          axis.text.x = element_text(angle = 45,vjust = 1.0,hjust = 1.0),
          #legend.direction = "top"
    )+stat_compare_means(method="anova")+
    scale_color_manual(values = c("#f89588","#f89588","#f89588","#9987ce","#9987ce","#9987ce","#9987ce","#9987ce","#efa666","#efa666","#63b2ee",
                                  "#63b2ee","#63b2ee","#63b2ee","#63b2ee"))
print(p)

dev.off()
