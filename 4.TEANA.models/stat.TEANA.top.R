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
        aa <- args[1]
        bb <- args[2]

pred.tbl<-read.table(file = aa,header =TRUE,row.names=1)

p2=data.frame(Bile=pred.tbl$Bile,Breast=pred.tbl$Breast,Colorectal=pred.tbl$Colorectal,Gastric=pred.tbl$Gastric,Lung=pred.tbl$Lung,Ovarian=pred.tbl$Ovarian,Pancreatic=pred.tbl$Pancreatic)

multiclass.roc(response=pred.tbl$Type,predictor=p2)

b=c("Bile","Breast","Colorectal","Gastric","Lung","Ovarian","Pancreatic")
predlab=b[apply(p2,1,which.max)]
confusionMatrix(data=as.factor(predlab),reference=as.factor(pred.tbl$Type),mode="everything")
