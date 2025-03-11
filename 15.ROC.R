
library(pROC)

rsFile="model.riskMatrix.txt"   
method="RF"          
setwd("")  

riskRT=read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)

CohortID=gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID=gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort=CohortID

for(Cohort in unique(riskRT$Cohort)){
  rt=riskRT[riskRT$Cohort==Cohort,]
  y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
  y=ifelse(y=="Control", 0, 1)
  
  roc1=roc(y, as.numeric(rt[,method]))    
  ci1=ci.auc(roc1, method="bootstrap")     
  ciVec=as.numeric(ci1)
  pdf(file=paste0("ROC.", Cohort, ".pdf"), width=5, height=4.75)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=Cohort)
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}


