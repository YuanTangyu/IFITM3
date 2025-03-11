
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("")               
rt=read.table("id.txt",sep="\t",header=T,check.names=F) 
rt=rt[is.na(rt[,"entrezID"])==F,]                      
gene=rt$entrezID

kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)             
pdf(file="barplot.pdf",width = 8,height = 8)
barplot(kk, drop = TRUE, showCategory =5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
pdf(file="bubble.pdf",width = 8,height = 8)
dotplot(kk,showCategory = 5,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
