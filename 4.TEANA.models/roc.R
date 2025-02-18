#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

library(pROC)
library(caret)
library(tidyverse)

args <- commandArgs(trailing=T)
        aa<- args[1]
        bb<-args[2]
		cc<-args[3]

setwd(bb)
test_pred<-read.table(file = aa,header =TRUE,row.names=1)

pdf(paste0(cc,".test.final.pdf"),h=5,w=5)
test_r <- roc( test_pred$Type, test_pred$pred)
test_r
#hist(test_r)
plot(test_r,legacy.axes=T,lwd= 2)
abline(v=0.95,col="gray60",lty=2)
auc=round(test_r$auc,3)
legend("bottomright",legend=c(paste0(cc," AUC=",auc)),lty=1,bty="n",cex=0.8)
print(ci.se(test_r))
print(ci.se(test_r,specificities=c(0.95,0.98,1.0)))
dev.off()

#data<-ifelse(data$Type!="Control",1,0)

#set.seed(123)
#head(data)

#legend( 'topleft', c('A', 'B'), col=c('red','blue'), lty=c(1,1), bty='n', cex=1.1 );

#dev.off();



