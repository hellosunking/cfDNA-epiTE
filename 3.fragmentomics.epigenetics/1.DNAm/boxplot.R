rm(list=ls())
#load packages
library(pROC)
library(cowplot)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggprism)

args <- commandArgs(trailing=T)
aa<-args[1]
bb<-args[2]
cc<-args[3]

pic<-list()
pdf(paste0(bb,".",cc,".boxplot.pdf") , width=2, height=2.5 )  
 
size.Low<-list()
size.High<-list()
size <- read.table(aa,header=T)

size.Alu<-size[size$Group==cc,]

size.min<-min(min(size.Alu$Low),min(size.Alu$High),min(size.Alu$Partial))
size.max<-max(max(size.Alu$Low),max(size.Alu$High),max(size.Alu$Partial))
    
    size.Alu<-size[size$Group==cc,]
    #pvalue<-t.test(x=size.Alu$Low,y=size.Alu$High,paired = T)$p.value
    pvalue1<-t.test(x=size.Alu$Low,y=size.Alu$Partial,paired = T)$p.value
	pvalue2<-t.test(x=size.Alu$Low,y=size.Alu$High,paired = T)$p.value
	pvalue3<-t.test(x=size.Alu$Partial,y=size.Alu$High,paired = T)$p.value
	df2 <- melt(size.Alu)
    p<- ggboxplot(df2, x = "variable", y = "value",
                   color  = "gray", add = "point",width = 0.4)+
      labs(title=paste0(bb, ", p1=",sprintf("%.3e",pvalue1), ", p2=",sprintf("%.3e",pvalue2), ", p3=",sprintf("%.3e",pvalue3)), x=NULL, y="DNA methylation")+
      geom_line(aes(group = sid), color = 'gray', lwd = 0.5) +  
      geom_point(aes(color = variable), size = 0.5, shape = 1,stroke = 1)+
      theme_bw()+
      theme(panel.grid = element_blank(),
            legend.position = "none",
            plot.title=element_text(size=10), 
            axis.line = element_line(colour = "#000000",linewidth = 0.3),
            axis.text = element_text(colour = "#000000"),
            axis.text.y  = element_text(angle = 90,vjust = 1.0,hjust = 0.5),
      )+
      scale_color_manual(values = c("#00468BFF","#5F559BFF","#ED0000FF")) +ylim(size.min,size.max)
  print(p)

dev.off()
