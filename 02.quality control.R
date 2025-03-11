getwd()
setwd(" ")
rm(list=ls())
options(stringsAsFactors = F) 
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(data.table)
library(dplyr)
library(harmony)
sce.all <- readRDS("sce.all.rds")

mito_genes=rownames(sce.all)[grep("^MT-", rownames(sce.all),ignore.case = T)] 
print(mito_genes) 
#sce.all=PercentageFeatureSet(sce.all, "^MT-", col.name = "percent_mito")
sce.all=PercentageFeatureSet(sce.all, features = mito_genes, col.name="percent_mito")
fivenum(sce.all@meta.data$percent_mito)


ribo_genes=rownames(sce.all)[grep("^Rp[sl]", rownames(sce.all),ignore.case = T)]
print(ribo_genes)
sce.all=PercentageFeatureSet(sce.all, features = ribo_genes, col.name = "percent_ribo")
fivenum(sce.all@meta.data$percent_ribo)

Hb_genes=rownames(sce.all)[grep("^Hb[^(p)]", rownames(sce.all),ignore.case = T)]
print(Hb_genes)
sce.all=PercentageFeatureSet(sce.all,  features = Hb_genes,col.name = "percent_hb")
fivenum(sce.all@meta.data$percent_hb)
head(sce.all@meta.data)

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
           "percent_ribo", "percent_hb")
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(sce.all, group.by="orig.ident", features=feats, pt.size=0, ncol=2) + 
  NoLegend()
w=length(unique(sce.all$orig.ident))/3+5;w
ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height=5)
feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(sce.all, group.by = "orig.ident", features=feats, pt.size=0, ncol=3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2
w=length(unique(sce.all$orig.ident))/2+5;w
ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)

p3=FeatureScatter(sce.all, "nCount_RNA", "nFeature_RNA", group.by="orig.ident", pt.size=0.5)
p3
ggsave(filename="Scatterplot.pdf",plot=p3)


if(F){
  selected_c <- WhichCells(sce.all, expression = nFeature_RNA > 500)
  selected_f <- rownames(sce.all)[Matrix::rowSums(sce.all@assays$RNA$counts > 0 ) > 3]
  sce.all.filt <- subset(sce.all, features = selected_f, cells = selected_c)
  dim(sce.all) 
  dim(sce.all.filt) 
}
sce.all.filt =  sce.all
C=subset(sce.all.filt,downsample=100)@assays$RNA$counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100

most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]

pdf("TOP50_most_expressed_gene.pdf",width=12)
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
        cex = 0.1, las = 1, 
        xlab = "% total count per cell", 
        col = (scales::hue_pal())(50)[50:1], 
        horizontal = TRUE)
dev.off()
rm(C)

selected_mito <- WhichCells(sce.all.filt, expression = percent_mito < 25)
selected_ribo <- WhichCells(sce.all.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(sce.all.filt, expression = percent_hb < 1 )
length(selected_hb)
length(selected_ribo)
length(selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_mito)
sce.all.filt <- subset(sce.all.filt, cells = selected_ribo)
sce.all.filt <- subset(sce.all.filt, cells = selected_hb)
dim(sce.all.filt)
table(sce.all.filt$orig.ident)
length(sce.all.filt$orig.ident)

feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
w=length(unique(sce.all.filt$orig.ident))/3+5;w 
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(sce.all.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
w=length(unique(sce.all.filt$orig.ident))/2+5;w 
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 

p2
sce.all.filt <- NormalizeData(sce.all.filt, 
                              normalization.method = "LogNormalize",
                              scale.factor = 1e4) 

sce.all.filt <- FindVariableFeatures(sce.all.filt)
p4 <- VariableFeaturePlot(sce.all.filt) 
p4

sce.all.filt <- ScaleData(sce.all.filt)

sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))
VizDimLoadings(sce.all.filt, dims = 1:4, reduction = "pca")
p <- VizDimLoadings(sce.all.filt, dims = 1:4, reduction = "pca")
ggsave("dim_loadings_plot.pdf", plot = p, width = 8, height = 10)

DimPlot(sce.all.filt, reduction = "pca") + NoLegend()
p <-DimPlot(sce.all.filt, reduction = "pca") + NoLegend()
ggsave("pca-2.pdf", plot = p, width = 8, height = 6)

DimHeatmap(sce.all.filt, dims = 1:12, cells = 500, balanced = TRUE)
p <- DimHeatmap(sce.all.filt, dims = 1:12, cells = 500, balanced = TRUE)

seuratObj <- RunHarmony(sce.all.filt, "orig.ident")
names(seuratObj@reductions)

seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=F ) 
seuratObj <- RunTSNE(seuratObj, dims = 1:15, 
                     reduction = "harmony")
DimPlot(seuratObj,reduction = "tsne",label=F ) 


sce.all.filt=seuratObj
sce.all.filt <- FindNeighbors(sce.all.filt, reduction = "harmony",
                              dims = 1:15) 

sce.all.filt.all=sce.all.filt

for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all.filt.all=FindClusters(sce.all.filt.all, #graph.name = "CCA_snn", 
                                resolution = res, algorithm = 1)
}
colnames(sce.all.filt.all@meta.data)
apply(sce.all.filt.all@meta.data[,grep("RNA_snn",colnames(sce.all.filt.all@meta.data))],2,table)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.01") + 
                   ggtitle("louvain_0.01"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.1") + 
                   ggtitle("louvain_0.1"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.2") + 
                   ggtitle("louvain_0.2"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_low.pdf",width = 14)

p1_dim=plot_grid(ncol = 3, DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.8") + 
                   ggtitle("louvain_0.8"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.1") + 
                   ggtitle("louvain_1"), DimPlot(sce.all.filt.all, reduction = "umap", group.by = "RNA_snn_res.0.3") + 
                   ggtitle("louvain_0.3"))
ggsave(plot=p1_dim, filename="Dimplot_diff_resolution_high.pdf",width = 18,height = 8)

p2_tree=clustree(sce.all.filt.all@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="Tree_diff_resolution.pdf",height = 8)
table(sce.all.filt.all@active.ident) 

saveRDS(sce.all.filt.all,"sce.all.filt.all.rds")
saveRDS(sce.all,"sce.all-2.rds")

