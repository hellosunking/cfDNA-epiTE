#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

library(pROC)
library(caret)
library(tidyverse)

argv <- commandArgs(trailing=T)
pred.tbl <- read.table(argv[1], head=T, row.names=1)

p2=data.frame(Bile=pred.tbl$Bile, Breast=pred.tbl$Breast, Colorectal=pred.tbl$Colorectal,
			  Gastric=pred.tbl$Gastric, Lung=pred.tbl$Lung, Ovarian=pred.tbl$Ovarian, Pancreatic=pred.tbl$Pancreatic)

b = c("Bile","Breast","Colorectal","Gastric","Lung","Ovarian","Pancreatic")
predlab = b[apply(p2,1,which.max)]
cm = confusionMatrix(data=as.factor(predlab), reference=as.factor(pred.tbl$Type))

cat("Prediction\\Reference\t", file=argv[2])
write.table( cm$table, file=argv[2], append=T, quote=F, sep="\t" )

