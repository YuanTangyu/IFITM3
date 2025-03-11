
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
#library(ROGUE)
library(clustree)
library(harmony)
library(SingleR)
library(celldex)

setwd(" ")
pbmc = readRDS("sce.all.filt.all.rds")
load("ref_Human_all.RData")
testdata = GetAssayData(object = pbmc@assays$RNA,layer = "counts")
clusters <- pbmc@meta.data$seurat_clusters
table(ref_Human_all@colData@listData[["label.main"]])
cellpred <-SingleR(test = testdata,ref = ref_Human_all,clusters = clusters,assay.type.test = "logcounts",
                   labels = ref_Human_all@colData@listData[["label.main"]],assay.type.ref = "logcounts")
celltype = data.frame(clusterID=rownames(cellpred),celltype= cellpred$labels,stringsAsFactors = F)
pbmc@meta.data$SingleR = "NA"
for (i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc$seurat_clusters==celltype$clusterID[i]),'SingleR']<-celltype$celltype[i]
}
p = plotScoreHeatmap(cellpred)
ggsave("1.pdf",p,width = 12,height = 6)

p1 <- DimPlot(pbmc,group.by = "SingleR",label = T,reduction = "tsne")
ggsave("tsne.pdf",p1,width = 12,height=6)

saveRDS(pbmc,"6pbmc_SingleR.rds")


