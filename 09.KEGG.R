
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("")        
rt=read.table("id.txt",sep="\t",header=T,check.names=F) 
rt=rt[is.na(rt[,"entrezID"])==F,]                    
gene=rt$entrezID

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05) 
write.table(kk,file="KEGGId.txt",sep="\t",quote=F,row.names = F)                  

pdf(file="barplot.pdf",width = 9,height = 10)
barplot(kk, drop = TRUE, showCategory = 20)
dev.off()

pdf(file="bubble.pdf",width = 9,height = 10)
dotplot(kk, showCategory = 20)
dev.off()

