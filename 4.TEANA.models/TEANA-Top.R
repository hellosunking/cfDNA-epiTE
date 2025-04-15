library(pROC)
library(caret)
library(tidyverse)
library(foreach)
library(doParallel)

args <- commandArgs(trailing=T)
data <- read.table( args[1], header=T, row.names=1)
setwd( args[2] )
registerDoParallel(cores=as.integer(detectCores()/10))

data <- data[data$Type!="Control",]
b=c("Bile","Breast","Colorectal","Gastric","Lung","Ovarian","Pancreatic")

for (i in c(1:7)){
	set<-data
	set$Type <- ifelse(set$Type!=b[i],"aControl",b[i])
	test_auc_value <- {}
	train_auc_value <- {}
	prime_numbers <- c(71,89,37,11,53,41,7,29,61,13)

	for(rep in c(1:10)) {
		set.seed(prime_numbers[rep])
		folds <- createFolds(set$Type, k = 10)
		foreach(n=1:10, .packages = c("caret", "pROC"), .combine = 'c', .verbose = F ) %dopar% {
			train_indices <- unlist(folds[-n])
			test_indices <- unlist(folds[n])
			train <- set[train_indices, ]
			test <- set[test_indices, ]

			ctrl <- trainControl(method="repeatedcv",number=10,repeats=10,savePredictions=T,
								 classProbs=T,search="random",summaryFunction=twoClassSummary,verboseIter=F)

			# train the model for each parameter combination in the grid, using CV to evaluate
			param_grid <-expand.grid(n.trees= c(50,100,150,200,300),
							interaction.depth = c(1,3,5), shrinkage=c(0.01,0.05,0.1,0.2,0.3), n.minobsinnode = c(10,15,20))
			set.seed(7)
			gbm_grid<-train(Type ~ ., data=train, method="gbm",	trControl=ctrl,
				#preProcess = c("corr", "nzv"),
				tuneGrid=param_grid, metric="ROC")

			param_grid <- expand.grid( n.trees=gbm_grid$bestTune$n.trees,
						 interaction.depth=gbm_grid$bestTune$interaction.depth,
						 shrinkage=gbm_grid$bestTune$shrinkage,
						 n.minobsinnode=gbm_grid$bestTune$n.minobsinnode )

			dir.create(b[i])
			set.seed(7)
			model <- train( Type ~ ., data=train, method="gbm",
						   trControl=trainControl(method = "none", savePredictions=T,
									classProbs=T, summaryFunction=twoClassSummary, verboseIter=F),
							tuneGrid=param_grid, metric="ROC" )
			saveRDS(model, paste0(b[i],"/GBM.model",n,"_rep",rep,".rds"))

			##test
			test_pred <- predict(model, newdata=test, type="prob")
			test_Type=ifelse(test$Type=="aControl",0,1)
			test_pred.t <- data.frame(sid=rownames(test), Type=test_Type, pred=test_pred[,2])
			#test_r <- roc(test_Type, test_pred[,2])
			write.table(test_pred.t, file=paste0(b[i],"/test_pred",n,"_rep",rep,".txt"),
							sep="\t", col.names=T, row.names=F, quote=F)
		}
	}
}

stopImplicitCluster()
