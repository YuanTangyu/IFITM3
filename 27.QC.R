
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(clustree)
library(harmony)

logFCfilter=0.585       
adjPvalFilter=0.05   

workDir=""
setwd(workDir)

dirs=list.dirs(workDir)
dirs_sample=dirs[-1]
names(dirs_sample)=gsub(".+\\/(.+)", "\\1", dirs_sample)
counts <- Read10X(data.dir = dirs_sample)
pbmc = CreateSeuratObject(counts, min.cells=3, min.features=100)

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pdf(file="01.featureViolin.pdf", width=10, height=6.5)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
dev.off()
pbmc=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 15)    

pdf(file="01.featureCor.pdf", width=13, height=7)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pbmc <- NormalizeData(object=pbmc, normalization.method="LogNormalize", scale.factor=10000)
pbmc <- FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=1500)
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf", width=10, height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()


pbmc=ScaleData(pbmc)       
pbmc=RunPCA(object= pbmc, npcs=20, pc.genes=VariableFeatures(object=pbmc))     
pbmc=RunHarmony(pbmc, "orig.ident")

pdf(file="02.pcaGene.pdf", width=10, height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures=20)
dev.off()

pdf(file="02.PCA.pdf", width=7.5, height=5)
DimPlot(object=pbmc, reduction="pca")
dev.off()

pdf(file="02.pcaHeatmap.pdf", width=10, height=8)
DimHeatmap(object=pbmc, dims=1:4, cells=500, balanced=TRUE, nfeatures=30, ncol=2)
dev.off()

pbmc <- JackStraw(object=pbmc, num.replicate=100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object=pbmc, dims=1:20)
dev.off()

pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)    

pbmc <- FindClusters(pbmc, resolution=seq(0.5, 1.2, by=0.1))
pbmc <- FindClusters(object = pbmc, resolution=0.6)
pdf(file="03.cluster.pdf", width=7, height=6)
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)          
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    
dev.off()
write.table(pbmc$seurat_clusters,file="03.Cluster.txt",quote=F,sep="\t",col.names=F)

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(file="03.clusterHeatmap.pdf",width=15, height=15)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()
pbmc_for_SingleR <- GetAssayData(pbmc, layer="data")
clusters<-pbmc@meta.data$seurat_clusters
ref1=get(load("ref_Human_all.RData"))
singler=SingleR(test=pbmc_for_SingleR, ref =ref1,
                labels=ref1$label.main, clusters = clusters)

clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
clusterAnn$labels=gsub("_", " ", clusterAnn$labels)
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
cellAnn=c()
for(i in 1:length(pbmc$seurat_clusters)){
  index=pbmc$seurat_clusters[i]
  cellAnn=c(cellAnn, clusterAnn[index,2])
}
cellAnnOut=cbind(names(pbmc$seurat_clusters), cellAnn)
colnames(cellAnnOut)=c("id", "labels")
write.table(cellAnnOut, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

newLabels=gsub("_", " ", singler$labels)
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pdf(file="04.cellAnn.pdf", width=7.5, height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)           
dev.off()
Type=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
names(Type)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=Type, col.name="Type")
pdf(file="04.group.cellAnn.pdf", width=11, height=6)
TSNEPlot(object = pbmc, pt.size = 1, label = TRUE, split.by="Type") 
dev.off()
save(pbmc, cellAnn, file="Seurat.Rdata")

monocle.matrix=as.matrix(pbmc@assays$RNA$data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name=row.names(monocle.matrix), row.names=row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

clusterAnn=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)

pdf(file="06.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State")
dev.off()
pdf(file="06.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
pdf(file="06.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
pdf(file="06.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()


