
library(limma)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(linkET)

expFile="normalize.txt"   
geneFile="modelGene.list.txt"      
immFile="CIBERSORT-Results.txt"    
setwd("")    

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

group=gsub("(.*)\\_(.*?)", "\\", colnames(data))
data=data[,group=="Treat",drop=F]
data=t(data)

immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]
immune=immune[,apply(immune,2,sd)>0]
geneLists=list()
for(i in 1:ncol(data)){geneLists[[colnames(data)[i]]]=i}

geneCor=data.frame()
for(cell in colnames(immune)){
  if(sd(immune[,cell])==0){next}
  for(gene in colnames(data)){
    x=as.numeric(immune[,cell])
    y=as.numeric(data[,gene])
    corT=cor.test(x, y, method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    geneCor=rbind(geneCor, cbind(spec=gene, env=cell, r=cor, p=pvalue))
  }
}
geneCor$r=as.numeric(geneCor$r)
geneCor$p=as.numeric(geneCor$p)
geneCor$pd=ifelse(geneCor$p<0.05, ifelse(geneCor$r>0, "Postive", "Negative"), "Not")
geneCor$r=abs(geneCor$r)
geneCor=geneCor %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                                    labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6",">= 0.6")))

qcorPlot=qcorrplot(correlate(immune, method="spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = geneCor, 
              curvature = nice_curvature()) +
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  scale_size_manual(values = c(0.5, 1.5, 2, 3)) +
  scale_colour_manual(values = c("#1B9E77", "#CCCCCC99", "#D95F02")) +
  guides(size = guide_legend(title = "abs(Cor)",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "pvalue", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Cell-cell cor", order = 3))

pdf(file="linkET.pdf", width=9, height=7)
print(qcorPlot)
dev.off()


