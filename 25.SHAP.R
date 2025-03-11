
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)


set.seed(12345)  
inputFile="merge.normalize.txt"      
geneFile="interGenes.txt"     
setwd("")     

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]
row.names(data)=gsub("-", "_", row.names(data))

data=t(data)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
data=as.data.frame(data)
data$Type=group

inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

yTestClass=test$Type
yTest=ifelse(yTestClass=="Control", 0, 1)
test<-test[,-ncol(test)]

control=trainControl(method="repeatedcv", number=5, savePredictions=TRUE)
methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)

modelList=list()
AUC=c()
ROCcolor=rainbow(nrow(methodRT))
for(i in 1:nrow(methodRT)){
	name=methodRT[i,"Name"]
	method=methodRT[i,"Method"]
	model=train(Type ~ ., data = train, method=method, trControl = control)
	if(name=="SVM"){
		model=train(Type ~ ., data = train, method=method, prob.model=TRUE, trControl = control)
	}
	pred=predict(model, newdata=test, type="prob")
	roc=roc(yTest, as.numeric(pred[,2]))
	AUC=c(AUC, paste0(name, ': ', sprintf("%.03f",roc$auc)))
	modelList[[method]]=as.numeric(roc$auc)
	if(i==1){
		pdf(file="ROC.pdf", width=5.5, height=5)
		plot(roc, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3)
	}else{
		plot(roc, print.auc=F, legacy.axes=T, main="", col=ROCcolor[i], lwd=3, add=T)
	}
}
legend('bottomright', AUC, col=ROCcolor, lwd=3, bty = 'n', cex=0.9)
dev.off()

aucValue=unlist(modelList)
bestMethod=names(which(aucValue==max(aucValue)))
bestMethod
train$Type=ifelse(train$Type=="Control", 0, 1)
model=train(Type ~., data = train, method = bestMethod, trControl=control)

fit=permshap(model, train[,-ncol(train)])
#fit=kernelshap(model, train[,-ncol(train)])      
shp <- shapviz(fit, X_pred = train[,-ncol(train)], X=train[,-ncol(train)], interactions=T)
important=sort(colMeans(abs(shp$S)), decreasing=T)
showVars=names(important)
pdf(file="barplot.pdf", width=6, height=6)
sv_importance(shp, kind="bar", show_numbers=TRUE)+theme_bw()
dev.off()
pdf(file="bee.pdf", width=7, height=6)
sv_importance(shp, kind = "bee", show_numbers=TRUE)+theme_bw()
dev.off()
pdf(file="dependence.pdf", width=9, height=6)
sv_dependence(shp, v = showVars)+theme_bw()
dev.off()
pdf(file="waterfall.pdf", width=7, height=5)
sv_waterfall(shp, row_id = 1)
dev.off()
pdf(file="force.pdf", width=9, height=5)
sv_force(shp, row_id = 1)
dev.off()


