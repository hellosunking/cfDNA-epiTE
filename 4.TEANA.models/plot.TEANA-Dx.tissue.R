#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

rm(list=ls())
library(pROC);library(dplyr);library(verification);library(binom)
args <- commandArgs(trailing=T)
predict <- read.table(args[1], header=T, row.names=1)
sample  <- read.table(args[2], header=T, row.names=1, sep="\t")
sample  <- sample[ rownames(predict), ]
pdf( args[3], h=5, w=5 )

#a <- merge(predict,size,by = intersect(names(size)[1],names(predict)[1]))[,c(1,2,3,5)]
#a <- na.omit(a)
a = data.frame( predict$Type, predict$pred, sample$Group )
colnames(a) = c("Type", "pred", "Group")

pROC.p = function( rr ) {
	v  = var( rr );
	b  = rr$auc - 0.5;
	se = sqrt(v);
	z  = b / se;
	p  = 2 * pt(-abs(z), df=Inf);
	p;
}

b <- as.vector(row.names(table(sample$Group)))
num <- length(b)

col=c("gray","#df9e9b","#eaaa60","#d8e7ca","#BC3C28","#99badf","#999acd","#fcd1ec")
results <- list()
ci.auc.low=vector()
ci.auc.up=vector()

auc=vector()
pvalue=vector()
p=vector()
for (i in c(2:num) ){
	c <- subset(a, Group==b[i] | Group== "0.Healthy" )
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

	cat(paste0(b[i], ": sens at 95% spec = ", sensitivity, ", 95% CI = [", ci$lower, ", ", ci$upper, "]\n"))

	auc[i]=round(test_r$auc,3)
	ci.auc(test_r)
	ci.auc.low[i]=round(ci.auc(test_r)[1],3)
	ci.auc.up[i] =round(ci.auc(test_r)[3],3)
	p[i] = sprintf("%.1e", pROC.p(test_r));
	
#	print(ci.se(test_r, specificities=0.95))
	if(i==2){
		plot(test_r, col=col[i], legacy.axes=T, lwd=2)
	}else{
		plot(test_r, add=T, col=col[i], lwd=2)
	}
}

abline(v=0.95, col="gray60", lty=2)
legend("bottomright",legend=c(paste0(b[2:num]," AUC=",auc[2:num],",p=",p[2:num])),col=col,lty=1,bty="n",cex=0.7)
dev.off()

