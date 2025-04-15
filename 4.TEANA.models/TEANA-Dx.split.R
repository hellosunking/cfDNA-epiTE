library(pROC)
library(caret)
library(tidyverse)
library(foreach)
library(doParallel)

args <- commandArgs(trailing=T)
all <- read.table( args[1], header=T, row.names=1)
test.list <- read.table( args[2])
setwd( args[3] )
set.seed(123)

validation <- all[test.list[,1],]
## get the training subset
data <- all[!rownames(all) %in% test.list[,1],]

test_auc_value <- {}
train_auc_value <- {}
prime_numbers <- c(31,43,47,59,67,23,19,73,79,83)
registerDoParallel(cores=as.integer(detectCores()/10))

for(rep in 1:10) {
	set.seed(prime_numbers[rep])
	folds <- createFolds(data$Type, k = 10)
	foreach(i = 1:10, .packages = c("caret", "pROC"), .verbose = F ) %dopar% {
		train_indices <- unlist(folds[-i])
		train <- data[train_indices, ]
		test_indices <- unlist(folds[i])
		test <- data[ test_indices, ]

		ctrl <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions=T,
							 classProbs=T, search="random", summaryFunction=twoClassSummary, verboseIter=F)
		#train the model for each parameter combination in the grid, using CV to evaluate
		param_grid <- expand.grid( n.trees= c(50,100,200,300), interaction.depth = c(1,3,5,7),
								shrinkage=c(0.01,0.05,0.1,0.2,0.3), n.minobsinnode = c(5,10,15,20,30) )

		set.seed(123)
		gbm_grid <- train(Type ~ ., data=train, method="gbm", trControl=ctrl, tuneGrid=param_grid, metric="ROC")
		param_grid <- expand.grid( n.trees=gbm_grid$bestTune$n.trees,
						 interaction.depth=gbm_grid$bestTune$interaction.depth,
						 shrinkage=gbm_grid$bestTune$shrinkage,
						 n.minobsinnode=gbm_grid$bestTune$n.minobsinnode )

		set.seed(123)
		model<-train( Type ~ ., data=train, method="gbm",
						trControl=trainControl(method="none", savePredictions=T, classProbs=T,
									summaryFunction=twoClassSummary, verboseIter=F),
						tuneGrid=param_grid, metric="ROC" )
		saveRDS(model, paste0("GBM.model",i,"_rep",rep,".rds") )

		## predict the testing part. Not necessary here
		test_pred <- predict(model, newdata=test, type="prob")
		test_Type = ifelse(test$Type=="Control",0,1)
		test_pred.t <- data.frame(sid=rownames(test), Type=test_Type, pred=test_pred$cancer)
		write.table(test_pred.t, file=paste0("ncv_pred",i,"_rep",rep,".txt"), sep="\t", col.names=T, row.names=F, quote=F)

		## predict the non-overlapping validation subset
		val_pred <- predict(model, newdata=validation, type="prob")
		val_Type = ifelse(validation$Type=="Control",0,1)
		val_pred.t <- data.frame(sid=rownames(validation), Type=val_Type, pred=val_pred$cancer)
		write.table(val_pred.t, file=paste0("testing_pred",i,"_rep",rep,".txt"), sep="\t", col.names=T, row.names=F, quote=F)
	}
}

stopImplicitCluster()

