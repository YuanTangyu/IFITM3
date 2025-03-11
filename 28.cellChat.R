
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)

setwd("")

load("Seurat.Rdata")
expMatrix=as.matrix(pbmc@assays$RNA$data)
meta=as.data.frame(cellAnn)
colnames(meta)[1]="labels"
row.names(meta)=names(pbmc$seurat_clusters)

cellchat <- createCellChat(object = expMatrix, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use="labels")
groupSize <- as.numeric(table(cellchat@idents))     

CellChatDB <- CellChatDB.human      
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

pdf(file="COMM01.DatabaseCategory.pdf", width=7, height=5)
showDatabaseCategory(CellChatDB)
dev.off()

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)  
cellchat <- identifyOverExpressedInteractions(cellchat)    
cellchat <- projectData(cellchat, PPI.human)  

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)
df.net=subsetCommunication(cellchat)
write.table(file="COMM02.Comm.network.xls", df.net, sep="\t", row.names=F, quote=F)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
pdf(file="COMM03.cellNetworkCount.pdf", width=7, height=6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf(file="COMM04.cellNetworkWeight.pdf", width=7, height=6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
dev.off()

pdf(file="COMM05.singleCell.pdf", width=9, height=7.5)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,3), mgp=c(0,0,0), xpd=TRUE)
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.8,title.name=cel)
}
dev.off()

for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  pdf(file=paste0("COMM05.", cel, ".pdf"), width=6.5, height=5.5)
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.8,title.name=cel)
  dev.off()
}

pdf(file="COMM06.bubble.pdf", width=9.5, height=6)
netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x=45)
dev.off()


