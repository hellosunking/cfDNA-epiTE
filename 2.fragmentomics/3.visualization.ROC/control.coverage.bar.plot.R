rm(list=ls())
#setwd("C:/Users/Sunlab/Desktop")

#load packages
library(pROC)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggprism)
library(ggbeeswarm)
library(dplyr)
library(forcats)
#library(reshape)

args <- commandArgs(trailing=T)
aa<- args[1]
bb<-args[2]

size <- read.table(aa,header=T)
#motif <- read.xlsx(bb,header=T)
#Diversity <- read.xlsx(cc,header=T)
#Cover <- read.xlsx(dd,header=T)

pdf( paste0(bb,".pdf"), width=10, height=4 )

head(size)

g<-as.factor(table(size$sid))

n=length(g)


df <- melt(size)
fct <- reorder(df$variable, -df$value)
colour=rainbow(16)                   
df$G<-rep(c("SINE","LINE","SINE","LINE","LTR","DNA","LTR","LTR","DNA","LINE","DNA","DNA","LTR","DNA","LINE","LTR"), each = n)
ggplot(df,aes(x=fct,y=value,fill=G))+ 
  geom_bar(stat="summary",fun=mean,position="dodge",width = 0.6)+
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", width = 0.3,linewidth=0.3)+ 
  labs(x=NULL,y=NULL,title = "your title")+ 
  theme_bw()+theme(panel.grid = element_blank(),
                   legend.position = "none",
                   legend.title = element_text(size=10),
                   legend.key.size = unit(8, "pt"),
                   axis.line = element_line(colour = "#000000",linewidth = 0.3),
                   axis.text = element_text(colour = "#000000"),
                   axis.text.x = element_text(angle = 45,vjust = 1.0,hjust = 1.0),
                   #legend.direction = "top"
  )+
  scale_y_continuous(expand = c(0,0))+ 
  #coord_cartesian(ylim=c(0,41))+ 
  labs(fill = "subtype")+ 
  scale_fill_manual(values = c("#f8cb7f","#63b2ee","#f89588","#7cd6cf"))
  #stat_compare_means(label = "p.signif", method = "t.test",paired = T,
                    # ref.group = "Genomewide")
#r=compare_means(value~variable, df,method = "t.test",paired = T) 
dev.off()
