#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

library(pROC)
#library(caret)
library(tidyverse)

args <- commandArgs(trailing=T)
test_pred<-read.table(args[1], header=T, row.names=1)
pdf(args[3],h=5,w=5)
setwd(args[2])

test_r <- roc(test_pred$Type, test_pred$pred)
test_r
plot(test_r,legacy.axes=T,lwd= 2)
abline(v=0.95,col="gray60",lty=2)
auc=round(test_r$auc,3)
legend("bottomright",legend=c(paste0(args[3]," AUC=",auc)),lty=1,bty="n",cex=0.8)
#print(ci.se(test_r))
print(ci.se(test_r,specificities=0.95))
dev.off()

