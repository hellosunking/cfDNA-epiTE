rm(list=ls())
#load packages
library(pROC)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggprism)

args <- commandArgs(trailing=T)
aa<- args[1]
bb<-args[2]

pROC.p = function( rr ) {
	v  = var( rr );
	b  = rr$auc - 0.5;
	se = sqrt(v);
	z  = b / se;
	p  = 2 * pt(-abs(z), df=Inf);
	p;
}

data <- read.table(aa,header=T)

pdf( paste0(bb,".ROC.pdf"), width=6, height=6 )

n=length(data[1,])

g <- c("SINE.Alu", "LINE.L1" ,"SINE.MIR","LINE.L2","LTR.ERVL_MaLR","DNA.hAT_Charlie",
	 "LTR.ERV1","LTR.ERVL","DNA.TcMar_Tigger" ,"LINE.CR1", "DNA.hAT_Blackjack","DNA.hAT_Tip100","LTR.Gypsy",
	 "DNA.TcMar_Mariner","LINE.RTE_X","LTR.ERVK")

type <- as.vector(row.names(table(data$Group)))
print(type)
if (n==18){
#type <- as.vector(row.names(table(Cover$Group)))
	lenn<-length(type)
	for (ca in c(2:lenn)){
	data.cancer <- data[data$Group == type[ca],]
	data.c <- data[data$Group=="0.Control",]
	head(data.cancer)
	b=vector()
	for (i in c(1:16)){
		pvalue<-wilcox.test(data.cancer[,g[i]],data.c[,g[i]])$p.value
		if (pvalue<0.05){
			b<-append(b,g[i])
		}
	}
	col=rainbow(length(b),alpha = 0.6)
	auc=vector()
	p=vector()
	num=length(b)
	for (i in c(1:num)){
		data.type=data[data$Group == type[ca] | data$Group=="0.Control",]
		data.roc <- roc(data.type$Group , data.type[,b[i]])
		auc[i]=round(data.roc$auc,3)
		p[i] = sprintf("%.2e", pROC.p(data.roc))
		if(i==1){
			plot(data.roc,col=col[i],legacy.axes=F,lwd= 2,main = paste0(type[ca],",",bb,"-P<0.05"),cex.main=0.8)
		}else{
			plot(data.roc,add=TRUE,col=col[i],lwd= 2)
		}
	}
	legend("bottomright",legend=c(paste0(b," AUC=",auc,", p=",p)),col=col,lty=1,bty="n",cex=0.7)
	}
}else{
	lenn<-length(type)
	for (ca in c(2:lenn)){
		data.cancer <- data[data$Group == type[ca],]
		data.c <- data[data$Group=="0.Control",]
		head(data.cancer)
		b=vector()
		for (i in c(1:16)){
			pvalue<-wilcox.test(data.cancer[,g[i]],data.c[,g[i]])$p.value
			if (pvalue<0.05){
				b<-append(b,g[i])
			}
		}
		col=rainbow(length(b),alpha = 0.6)
		auc=vector()
		p=vector()
		num=length(b)
		delong.p=vector()
		data.type=data[data$Group == type[ca] | data$Group=="0.Control",]
		roc.overall<-roc(data.type$Group , data.type[,"Genomewide"])
		for (i in c(1:num)){
			#data.type=data[data$Group == type[ca] | data$Group=="0.Control",]
			data.roc <- roc(data.type$Group , data.type[,b[i]])
			auc[i]=round(data.roc$auc,3)
			p[i] = sprintf("%.2e", pROC.p(data.roc))
			delong=roc.test( data.roc,  roc.overall )
			delong.p[i] = sprintf("%.3e", delong$p.value)
			if(i==1){
				plot(data.roc,col=col[i],legacy.axes=F,lwd= 2,main = paste0(type[ca],",",bb,"-P<0.05"),cex.main=0.8)
			}else{
				plot(data.roc,add=TRUE,col=col[i],lwd= 2)
			}
		}
		plot(roc.overall,add=TRUE,col="#808080",lwd= 2)
		auc[num+1]=round(roc.overall$auc,3)
		p[num+1] = sprintf("%.3e", pROC.p(roc.overall));
		delong.p[num+1]=NA
		legend("bottomright",legend=c(paste0(append(b,"Genomewide")," AUC=",auc,", p=",p,",Delong.p=", delong.p)),col=c(col[1:16],"#808080"),lty=1,bty="n",cex=0.7)
	}
}
dev.off()

pdf(paste0(bb,".te.ROC.pdf"),h=6,w=6)

for (n in c(3:18)) {
	#subset(size,size$Group == "control" | size$Group == "HCC")
	num=length(type)
	b=type[2:num]
	col=rainbow(num,alpha = 0.6)
	# col=c("#df9e9b","#eaaa60","#d8e7ca","#BC3C28","#99badf","#999acd","#fcd1ec")
	auc=vector()
	p=vector()
	lenn=num-1
	for (i in c(1:lenn)){
		data.te <- data[data$Group == "0.Control" | data$Group == b[i],][,c(1,2,n)]
		head(data.te)
		te.name<-colnames(data.te)[3]

		data.roc <- roc(data.te$Group , data.te[,3])
		auc[i]=round(data.roc$auc,3)
		p[i] = sprintf("%.3e", pROC.p(data.roc));
		if(i==1){
			  plot(data.roc,col=col[i],legacy.axes=F,lwd= 2,main = paste0("bb.",te.name),cex.main=0.8)
		}else{
			  plot(data.roc,add=TRUE,col=col[i],lwd= 2)
		}
	}
	legend("bottomright",legend=c(paste0(b," AUC=",auc,", p=",p)),col=col,lty=1,bty="n",cex=0.7)
}
dev.off()
