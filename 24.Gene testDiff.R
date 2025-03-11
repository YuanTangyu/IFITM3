
library(limma)
library(ggpubr)

expFile="normalize.txt"    
geneFile="interGenes.txt"           
setwd("")   

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

Type=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=cbind(as.data.frame(data), Type)

group=levels(factor(data[,"Type"]))
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

for(i in colnames(data)[1:(ncol(data)-1)]){
  rt1=data[,c(i,"Type")]
  test=wilcox.test(rt1[,i] ~ Type)
  if(test$p.value<0.05){
    boxplot=ggviolin(rt1, x="Type", y=i, fill="Type",
                     xlab="",
                     ylab=paste0(i, " expression"),
                     legend.title="",
                     palette = c("#00FF92","#FF7F00"),
                     width=1, add = "boxplot", add.params = list(fill="white"))+ 
      stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")
    
    pdf(file=paste0("violin.",i,".pdf"), width=5, height=4.5)
    print(boxplot)
    dev.off()
  }
}


