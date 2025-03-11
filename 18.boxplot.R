
library(limma)
library(reshape2)
library(ggpubr)
library(PerformanceAnalytics)

expFile="normalize.txt"   
geneFile="modelGene.list.txt"    
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

Type=gsub("(.*)\\_(.*?)", "\\2", row.names(data))
treatData=data[Type=="CD",]
rt=cbind(as.data.frame(data), Type)

data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", fill = "Type",
            xlab="",
            ylab="Gene expression",
            legend.title="Type",
            palette = c("#0088FF", "#FF5555"), width=0.75)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()

pdf(file="cor.pdf", width=7, height=6.5)
chart.Correlation(treatData, histogram=TRUE, pch=19, method="pearson")
dev.off()

