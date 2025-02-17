library(pROC)
library(caret)
library(tidyverse)
library(foreach)
library(doParallel)


args <- commandArgs(trailing=T)

aa<- args[1]
bb<-args[2]

setwd(bb)	
data<-read.table(file = aa,header =TRUE,row.names=1)


set.seed(123)
test_auc_value <- {}
train_auc_value <- {}
prime_numbers <- c(71,89,37,11,53,41,7,29,61,13)
registerDoParallel(cores = detectCores() - 1)

for (rep in 1:10) {
	set.seed(prime_numbers[rep])
	folds <- createFolds(data$Type, k = 10)
	foreach(i = 1:10, .packages = c("caret", "pROC"), .verbose = TRUE ) %dopar% {
	train_indices <- unlist(folds[-i])
	test_indices <- unlist(folds[i])

	train <- data[train_indices, ]
	test <- data[test_indices, ]

		head(test)
	ctrl <- trainControl(method="repeatedcv",number=10,repeats=10,savePredictions=TRUE,classProbs=TRUE,search = "random",summaryFunction = twoClassSummary,verboseIter=TRUE)
# train the model for each parameter combination in the grid, using CV to evaluate
param_grid <-expand.grid(
        n.trees= c(50,100,200,300),
        interaction.depth = c(1,3,5,7),
        shrinkage=c(0.01,0.05,0.1,0.2,0.3),
        n.minobsinnode = c(5,10,15,20,30)
)

set.seed(123)
gbm_grid<-train(Type ~ .,data=train,
				method="gbm",
				trControl=ctrl,
				#preProcess = c("corr", "nzv"),
				tuneGrid=param_grid,
				metric="ROC"
							)
print(gbm_grid)


param_grid <-expand.grid(
						 n.trees=gbm_grid$bestTune$n.trees,
						 interaction.depth=gbm_grid$bestTune$interaction.depth,
						 shrinkage=gbm_grid$bestTune$shrinkage,
						 n.minobsinnode=gbm_grid$bestTune$n.minobsinnode
						 )

set.seed(123)
model<-train(
			 Type ~ .,data=train,
			 method="gbm",
			 trControl=trainControl(method = "none",
									savePredictions = TRUE,
									classProbs = TRUE,
									summaryFunction = twoClassSummary,
									verboseIter = TRUE),
			 #preProcess = c("corr", "nzv"),
			 tuneGrid=param_grid,
			 metric="ROC"
			 )
print(model)

name=paste0("GBM.model",i,"_rep",rep,".rds")
saveRDS(model, name)


##train
train_pred <- predict(model,
						  newdata = train,
						  #n.trees = 29,
						  type = "prob",
						  )
train_Type=ifelse(train$Type=="Control",0,1)
train_pred.t <- data.frame(sid = rownames(train),
						   Type = train_Type,
						   pred = train_pred$cancer)
train_r <- roc(train_Type, train_pred$cancer)
#plot(test_r)

##test
test_pred <- predict(model, 
					newdata = test, 
					#n.trees = 29,
					type = "prob")
test_Type=ifelse(test$Type=="Control",0,1)
test_pred.t <- data.frame(sid = rownames(test),
						   Type = test_Type,
						   pred = test_pred$cancer)
test_r <- roc(test_Type, test_pred$cancer)
#plot(test_r)

#write.table(train_pred.t,file = paste("train_pred.final",".txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = FALSE)
#write.table(test_pred.t,file = paste("test_pred.final",".txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = FALSE)

write.table(train_pred.t,file = paste("train_pred",i,"_rep",rep,".txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = FALSE)
write.table(test_pred.t,file = paste("test_pred",i,"_rep",rep,".txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = FALSE)

train_r
test_r

#test_auc_value<- append(test_auc_value,as.numeric(auc(test_r)))
#train_auc_value<- append(train_auc_value,as.numeric(auc(train_r)))
}
}

#mean(test_auc_value)
#mean(train_auc_value)
#stopCluster(cl)
stopImplicitCluster()
