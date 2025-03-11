
library(VennDiagram)

setwd("D:\\CD+Z\\04.Veen")

data=read.table("1.txt", header=F, sep="\t", check.names=F)
Wilcoxon = data[,1]
data=read.table("diff.txt", header=F, sep="\t", check.names=F)
limma = data[,1]
data=read.table("Monocyte.txt", header=F, sep="\t", check.names=F)
edgeR = data[,1]

venn.diagram(
  x = list('Monocyte_DEGs' = edgeR, 'CD_DEGs' = limma, 'palmitoylation' = Wilcoxon),
  filename = 'VN-2.png',
  fill = c("dodgerblue", "goldenrod1", "green"),
  width = 8,  
  height = 8,  
  units = "in"  
)

intersectGenes1 = intersect(Wilcoxon,limma)
intersectGenes = intersect(intersectGenes1,edgeR)

write.table(file="intersectGenes.txt", intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)
