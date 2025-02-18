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

data <- data[data$Type!="Control",]
#data<-ifelse(data$Type!="Control",1,0)

b=c("Bile","Breast","Colorectal","Gastric","Lung","Ovarian","Pancreatic")


#set.seed(123)
head(data)

registerDoParallel(cores = detectCores() - 1)
#registerDoParallel(cores = detectCores() - 1)
#data$Type <- ifelse(data$Type == "Control" ,'0','1')

for (i in c(1:7)){
set<-data
set$Type <- ifelse(set$Type!=b[i],"aControl",b[i])
test_auc_value <- {}
train_auc_value <- {}
prime_numbers <- c(71,89,37,11,53,41,7,29,61,13)

for (rep in c(1:10)) {
	set.seed(prime_numbers[rep])
	folds <- createFolds(set$Type, k = 10)
	foreach(n = 1:10, .packages = c("caret", "pROC"), .combine = 'c',.verbose = TRUE ) %dopar% {
		train_indices <- unlist(folds[-n])
		test_indices <- unlist(folds[n])

		train <- set[train_indices, ]
		test <- set[test_indices, ]
		head(test)
		
		ctrl <- trainControl(method="repeatedcv",number=10,repeats=10,savePredictions=TRUE,classProbs=TRUE,search = "random",summaryFunction = twoClassSummary,verboseIter=TRUE)

# train the model for each parameter combination in the grid, using CV to evaluate
param_grid <-expand.grid(
        n.trees= c(50,100,150,200,300),
        interaction.depth = c(1,3,5),
        shrinkage=c(0.01,0.05,0.1,0.2,0.3),
        n.minobsinnode = c(10,15,20)
)

set.seed(7)
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

dir.create(b[i])
#write.table(param_grid,file = paste(b[i],"/best.param",".txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = FALSE)
set.seed(7)
model<-train(
			 Type ~ .,data=train,
			 #x=train[,-1],
			 #x=data_train,
			 #y=train$Type,
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

name=paste0(b[i],"/GBM.model",n,"_rep",rep,".rds")
saveRDS(model, name)


##train
train_pred <- predict(model,
						  newdata = train,
						  #n.trees = 29,
						  type = "prob",
						  )
head(train_pred)
train_Type=ifelse(train$Type=="aControl",0,1)
train_pred.t <- data.frame(sid = rownames(train),
						   Type = train_Type,
						   pred = train_pred[,2])
train_r <- roc(train_Type, train_pred[,2])
#plot(test_r)

##test
test_pred <- predict(model, 
					newdata = test, 
					#n.trees = 29,
					type = "prob")
print(test_pred)
test_Type=ifelse(test$Type=="aControl",0,1)
test_pred.t <- data.frame(sid = rownames(test),
						   Type = test_Type,
						   pred = test_pred[,2])
test_r <- roc(test_Type, test_pred[,2])
#plot(test_r)

write.table(train_pred.t,file = paste(b[i],"/train_pred",n,"_rep",rep,".txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = FALSE)
write.table(test_pred.t,file = paste(b[i],"/test_pred",n,"_rep",rep,".txt",sep=""),sep = "\t",col.names = T,row.names = F,quote = FALSE)

print(train_r)
print(test_r)

}

#mean(test_auc_value)
#mean(train_auc_value)
}

}

stopImplicitCluster()
