
library(devtools)
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)
library(patchwork)
library(future)
library(doFuture)
plan(strategy = "multisession") 

options(stringsAsFactors = FALSE)

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


cellchat@netP$pathways

pathways.show <- c("CCL") 

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair
vertex.receiver = seq(1,4) 
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

pathways.show.all <- cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

plotGeneExpression(cellchat, signaling = "CCL")
plotGeneExpression(cellchat, signaling = "CCL", enriched.only = FALSE)


library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
library(ggalluvial)
netAnalysis_dot(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")
saveRDS(cellchat, file = "cellchat.rds")
sessionInfo()




