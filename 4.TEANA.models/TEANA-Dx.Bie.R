library(pROC)
library(caret)
library(tidyverse)
library(foreach)
library(doParallel)

args <- commandArgs(trailing=T)

data<-read.table(args[1], header=T, row.names=1)
setwd(args[2])

set.seed(123)
test_auc_value <- {}
train_auc_value <- {}
prime_numbers <- c(71,89,37,11,53,41,7,29,61,13)
registerDoParallel(cores = as.integer(detectCores() / 10) )

for (rep in 1:10) {
	set.seed(prime_numbers[rep])
	folds <- createFolds(data$Type, k = 10)
	foreach(i = 1:10, .packages = c("caret", "pROC"), .verbose=F ) %dopar% {
		train_indices <- unlist(folds[-i])
		test_indices <- unlist(folds[i])
		train <- data[train_indices, ]
		test <- data[test_indices, ]

		ctrl <- trainControl(method="repeatedcv", number=10, repeats=10, savePredictions=T, classProbs=T,
						 search = "random", summaryFunction=twoClassSummary, verboseIter=F)

		# train the model for each parameter combination in the grid, using CV to evaluate
		param_grid <- expand.grid( n.trees= c(50,100,200,300), interaction.depth = c(1,3,5,7),
									shrinkage=c(0.01,0.05,0.1,0.2,0.3), n.minobsinnode = c(5,10,15,20,30) )

		set.seed(123)
		gbm_grid <- train(Type ~ ., data=train, method="gbm", trControl=ctrl, tuneGrid=param_grid, metric="ROC")

		param_grid <- expand.grid( n.trees=gbm_grid$bestTune$n.trees, interaction.depth=gbm_grid$bestTune$interaction.depth,
									shrinkage=gbm_grid$bestTune$shrinkage, n.minobsinnode=gbm_grid$bestTune$n.minobsinnode )

		set.seed(123)
		model <- train( Type ~ ., data=train, method="gbm",
						trControl=trainControl(method = "none", savePredictions = TRUE, classProbs = TRUE,
									summaryFunction = twoClassSummary, verboseIter = F),
						#preProcess = c("corr", "nzv"),
						tuneGrid=param_grid, metric="ROC" )
		saveRDS(model, paste0("GBM.model",i,"_rep",rep,".rds"))

		## predict on testing set
		test_pred <- predict(model, newdata = test,	type = "prob")
		test_Type=ifelse(test$Type=="Control",0,1)
		test_pred.t <- data.frame(sid = rownames(test), Type = test_Type, pred = test_pred$cancer)
#		test_r <- roc(test_Type, test_pred$cancer)
		write.table(test_pred.t,file = paste("test_pred",i,"_rep",rep,".txt",sep=""),
					sep="\t", col.names=T, row.names=F, quote=F )
	}
}

stopImplicitCluster()

