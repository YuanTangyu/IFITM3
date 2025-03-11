library(Seurat)
library(devtools)
library(BiocManager)
library(ComplexHeatmap)
library(ggunchull)
library(jjAnno)
library(scRNAtoolVis)


setwd("")
pbmc3k <- readRDS("6pbmc_SingleR.rds")
cc<-pbmc3k
cc$group = sample(c("HC","CD"),ncol(cc),replace = T)
table(cc$SingleR,cc$group )
cell<-c("T_cells","B_cell","NK_cell","Monocyte")
length(cell)
clusterdeg=data.frame()
for (i in 1:4) {
  Idents(cc)="SingleR"
  deg=FindMarkers(cc,ident.1 = "HC",ident.2 = "CD",
                  group.by = "group",subset.ident =cell[i])

  deg$cluster=cell[i]
  clusterdeg=rbind(deg,clusterdeg)
}
clusterdeg$gene=rownames(clusterdeg)
write.table(clusterdeg, file = "clusterdeg.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

col<-c(rgb(0.557,0.459,0.655),rgb(0.588,0.459,0.435),rgb(0.776,0.569,0.549),rgb(0.388,0.545,0.733))
p1<-jjVolcano(diffData = clusterdeg,
              log2FC.cutoff = 0.585, 
              size  = 3.5, #设置点的大小
              #aesCol = c('blue','orange'), 
              tile.col = col,
              #col.type = "adjustP", 
              topGeneN = 5 
)
p1

clusterdeg <- FindAllMarkers(cc, only.pos = FALSE,
                             min.pct = 0.1,
                             logfc.threshold = 0,group.by="SingleR")
p2<-jjVolcano(diffData = clusterdeg,
              log2FC.cutoff = 0.585,
              size  = 2.5, 
              #aesCol = c('blue','orange'), 
              #tile.col = col, 
              #col.type = "adjustP", 
              topGeneN = 5 
)

ggsave("DEGs.pdf", plot = p2, width = 15, height = 8)
write.table(clusterdeg, file = "clusterdeg-quanbu.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
