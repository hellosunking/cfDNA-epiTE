rm(list=ls())
#setwd("C:/Users/Sunlab/Desktop")

#load packages
library(pROC)
library(ggplot2)
library(reshape2)
library(ggsci)
library(ggpubr)
library(ggprism)
##load data
args <- commandArgs(trailing=T)
aa<- args[1]
bb<-args[2]
size <- read.table(aa,header=T)

pdf( paste0(bb,".boxplot.pdf"), width=15, height=4.5 ) 
library(forcats)
library(ggbeeswarm)

g<-as.vector(row.names(table(size$Group)))

n=length(g)

mycolor<-pal_npg("nrc", alpha = 0.8)(n-1)

##size
head(size)
df2 <- melt(size)
df2 %>%
  mutate(name = fct_relevel(variable, "SINE.Alu",
                            "SINE.MIR",
                            "LTR.Gypsy",
                            "LTR.ERVL_MaLR",
                            "LTR.ERVL",
                            "LTR.ERVK",
                            "LTR.ERV1",
                            "LINE.RTE_X",
                            "LINE.L2",
                            "LINE.L1",
                            "LINE.CR1",
                            "DNA.TcMar_Tigger",
                            "DNA.TcMar_Mariner",
                            "DNA.hAT_Tip100",
                            "DNA.hAT_Charlie",
                            "DNA.hAT_Blackjack",
                            "Genomewide"
  )) %>%
  ggplot( aes(x = name, y = value,fill=Group,color=Group)) +
  geom_boxplot(outlier.shape = NA,color = "black")+
  theme_bw()+
  #geom_beeswarm(dodge.width=1.0) +
  labs(x=NULL, y=NULL, title = bb )+
  theme(panel.grid = element_blank(),
        axis.line = element_line(colour = "#000000",linewidth = 0.3),
        axis.text = element_text(colour = "#000000"),
        axis.text.x = element_text(angle = 45, vjust = 1.0,hjust = 1.0),
        #legend.direction = "top"
  )+
  scale_fill_manual(values = c("gray",mycolor))+
  labs(fill = "class")
#  +stat_compare_means(size=2.5,
#                     method = "wilcox.test",
#                     label = "p.signif")
dev.off()
