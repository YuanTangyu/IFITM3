
library(glmnet)
library(pROC)

expFile="normalize.txt"  
geneFile="7Gene.txt"      
setwd("")    

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="Control", 0, 1)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)

aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
  k=k+1
  roc1=roc(y, as.numeric(rt[x,]))  
  if(k==1){
    pdf(file="ROC.genes.pdf", width=5, height=4.5)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3, add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

legend(x = 0.5, y = 0.5, aucText, lwd=3, bty="n", cex=0.8, col=bioCol[1:(ncol(rt)-1)])

dev.off()

