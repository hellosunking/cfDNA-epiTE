#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

rm(list=ls())
library(pROC);library(dplyr);library(verification)
library(binom)
args <- commandArgs(trailing=T)

aa <- args[1]
bb <- args[2]
cc <- args[3]

test <- read.table(file =  aa,header=T)
size <- read.table(bb,header=T,sep = "\t")

a <- merge(test,size,by = intersect(names(size)[1],names(test)[1]))[,c(1,2,3,5)]
a <- na.omit(a)
head(a)

pROC.p = function( rr ) {
	v  = var( rr );
	b  = rr$auc - 0.5;
	se = sqrt(v);
	z  = b / se;
	p  = 2 * pt(-abs(z), df=Inf);
	p;
}

name=paste0(cc,".test",".stage.roc.pdf")
pdf(name,h=5,w=5)

#b=c("BRCA","COREAD","ESCA","LIHC","NSCLC","PACA","STAD")
b=c("Stage I","Stage II","Stage III","Stage IV","Stage X")

col=c("#df9e9b","#eaaa60","#BC3C28","#99badf","#999acd")
#col<-rainbow(le)
results <- list()
ci.auc.low=vector()
ci.auc.up=vector()

auc=vector()
pvalue=vector()
p=vector()
for (i in c(1:5) ){
	c <- a[a$stage==b[i] | a$stage== "0.Healthy",]
	test_r<-roc(c$Type,c$pred)
    
	# Get the threshold corresponding to 95% specificity
    coords <- coords(test_r, "all", ret = c("threshold", "sensitivity", "specificity", "precision"))
    target_coord <- coords[coords[, "specificity"] >= 0.95, ][which.max(coords[coords[, "specificity"] >= 0.95, "sensitivity"]), ]

	threshold_raw <- target_coord["threshold"]
	sensitivity_raw <- target_coord["sensitivity"]
	specificity_raw <- target_coord["specificity"]

	threshold <- threshold_raw$threshold
	sensitivity <- sensitivity_raw$sensitivity
	specificity <- specificity_raw$specificity

	
	# Calculate the number of positive and negative cases detected at the determined threshold
	num_positive_detected <- sum(c$pred >= threshold & c$Type == 1)
	num_negative_detected <- sum(c$pred >= threshold & c$Type == 0)

	
	# Calculate 95% confidence intervals for sensitivity
	ci <- binom.confint(num_positive_detected, length(test_r$case), conf.level = 0.95, methods = "exact")
	#ci <- binom.test(num_positive_detected, length(test_r$case), conf.level = 0.95)
	results[[i]] <- list(sensitivity = sensitivity, specificity = specificity,
						 ci_lower = ci$lower, ci_upper = ci$upper,
						 num_positive_detected = num_positive_detected)
	cat(paste0("Group ", i, ": Sensitivity at 95% specificity = ", sensitivity,", Threshold = ", threshold,
			   ", 95% CI = [", ci$lower, ", ", ci$upper,
			   "], Positive cases detected = ", num_positive_detected,
			   ", Negative cases detected = ",num_negative_detected, "\n"))

	auc[i]=round(test_r$auc,3)
    ci.auc(test_r)
	ci.auc.low[i]=round(ci.auc(test_r)[1],3)
	ci.auc.up[i]=round(ci.auc(test_r)[3],3)
	p[i] = sprintf("%.3e", pROC.p(test_r))
	
	print(ci.se(test_r,specificities=c(0.95,0.98)))

	if(i==1){
		plot(test_r,col=col[i],legacy.axes=T,lwd=2)
	}else{
		plot(test_r,add=TRUE,col=col[i],lwd= 2)
	}
}

abline(v=0.95,col="gray60",lty=2)
legend("bottomright",legend=c(paste0(b," AUC=",auc,",p=",p)),col=col,lty=1,bty="n",cex=0.7)
dev.off()

