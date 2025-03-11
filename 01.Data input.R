rm(list=ls())
options(stringsAsFactors=F)
library(Seurat)
library(data.table)
library(dplyr)

setwd(" ")
dir='GSE214695_RAW/'
fs=list.files('GSE214695_RAW/','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify=T)[,1]

lapply(unique(samples),function(x){
  #x=unique(samples)[1]
  y=fs[grepl(x,fs)]
  folder=paste0("GSE214695_RAW/",paste(str_split(y[1],'_',simplify=T)[,1:2],collapse="_"))
  dir.create(folder,recursive=T)
  file.rename(paste0("GSE214695_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
  file.rename(paste0("GSE214695_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0("GSE214695_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
})

dir='GSE214695_RAW/'
samples=list.files(dir)
samples
sceList=lapply(samples,function(pro){
  #pro=samples[1]
  print(pro)
  tmp=Read10X(file.path(dir,pro))
  if(length(tmp)==2){
    ct=tmp[[1]]
  }else{ct=tmp}
  sce=CreateSeuratObject(counts=ct,
                          project=pro,
                          min.cells=5,
                          min.features=300)
  return(sce)
})
View(sceList)


do.call(rbind,lapply(sceList,dim))
sce.all=merge(x=sceList[[1]],
              y=sceList[-1],
              add.cell.ids=samples)
names(sce.all@assays$RNA@layers)

sce.all[["RNA"]]$counts
#Alternateaccessor function with the same result
LayerData(sce.all, assay = "RNA", layer ="counts")
sce.all
sce.all <- JoinLayers(sce.all)
sce.all

dim(sce.all[["RNA"]]$counts )
as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 
length(sce.all$orig.ident)

library(stringr)
phe=sce.all@meta.data
table(phe$orig.ident)
View(phe)

phe$group=str_split(phe$orig.ident,'[_]',simplify=T)[,2]

phe$sample=phe$orig.ident
phe$sample=gsub("GSM6614348_HC-1","HC",phe$sample)
phe$sample=gsub("GSM6614349_HC-2","HC",phe$sample)
phe$sample=gsub("GSM6614350_HC-3","HC",phe$sample)
phe$sample=gsub("GSM6614351_HC-4","HC",phe$sample)
phe$sample=gsub("GSM6614352_HC-5","HC",phe$sample)
phe$sample=gsub("GSM6614353_HC-6","HC",phe$sample)

phe$sample=gsub("GSM6614360_CD-1","CD",phe$sample)
phe$sample=gsub("GSM6614361_CD-2","CD",phe$sample)
phe$sample=gsub("GSM6614362_HC-3","CD",phe$sample)
phe$sample=gsub("GSM6614363_HC-4","CD",phe$sample)
phe$sample=gsub("GSM6614364_HC-5","CD",phe$sample)
phe$sample=gsub("GSM6614365_HC-6","CD",phe$sample)
sce.all@meta.data=phe
View(phe)

phe$patient = phe$orig.ident
table(phe$patient)
phe$patient = gsub("GSM6614348_HC-1", "Patient1", phe$patient)
phe$patient = gsub("GSM6614349_HC-2", "Patient2", phe$patient)
phe$patient = gsub("GSM6614350_HC-3", "Patient3", phe$patient)
phe$patient = gsub("GSM6614351_HC-4", "Patient4", phe$patient)
phe$patient = gsub("GSM6614352_HC-5", "Patient5", phe$patient)
phe$patient = gsub("GSM6614353_HC-6", "Patient6", phe$patient)

phe$patient = gsub("GSM6614360_CD-1", "Patient7", phe$patient)
phe$patient = gsub("GSM6614361_CD-2", "Patient8", phe$patient)
phe$patient = gsub("GSM6614362_CD-3", "Patient9", phe$patient)
phe$patient = gsub("GSM6614363_CD-4", "Patient10", phe$patient)
phe$patient = gsub("GSM6614364_CD-5", "Patient11", phe$patient)
phe$patient = gsub("GSM6614365_CD-6", "Patient12", phe$patient)
sce.all@meta.data=phe
View(phe)


phe$Gender = phe$orig.ident
table(phe$Gender)
phe$Gender = gsub("GSM6614348_HC-1", "Male", phe$Gender)
phe$Gender = gsub("GSM6614349_HC-2", "Male", phe$Gender)
phe$Gender = gsub("GSM6614350_HC-3", "Female", phe$Gender)
phe$Gender = gsub("GSM6614351_HC-4", "Male", phe$Gender)
phe$Gender = gsub("GSM6614352_HC-5", "Male", phe$Gender)
phe$Gender = gsub("GSM6614353_HC-6", "Female", phe$Gender)

phe$Gender = gsub("GSM6614360_CD-1", "Female", phe$Gender)
phe$Gender = gsub("GSM6614361_CD-2", "Male", phe$Gender)
phe$Gender = gsub("GSM6614362_CD-3", "Female", phe$Gender)
phe$Gender = gsub("GSM6614363_CD-4", "Male", phe$Gender)
phe$Gender = gsub("GSM6614364_CD-5", "Male", phe$Gender)
phe$Gender = gsub("GSM6614365_CD-6", "Male", phe$Gender)
sce.all@meta.data=phe
View(phe)

phe$Tissue = phe$orig.ident
table(phe$Tissue)
phe$Tissue = gsub("GSM6614348_HC-1", "sigma", phe$Tissue)
phe$Tissue = gsub("GSM6614349_HC-2", "sigma", phe$Tissue)
phe$Tissue = gsub("GSM6614350_HC-3", "sigma", phe$Tissue)
phe$Tissue = gsub("GSM6614351_HC-4", "sigma", phe$Tissue)
phe$Tissue = gsub("GSM6614352_HC-5", "sigma", phe$Tissue)
phe$Tissue = gsub("GSM6614353_HC-6", "sigma", phe$Tissue)

phe$Tissue = gsub("GSM6614360_CD-1", "ascending colon", phe$Tissue)
phe$Tissue = gsub("GSM6614361_CD-2", "sigma", phe$Tissue)
phe$Tissue = gsub("GSM6614362_CD-3", "descending colon", phe$Tissue)
phe$Tissue = gsub("GSM6614363_CD-4", "sigma", phe$Tissue)
phe$Tissue = gsub("GSM6614364_CD-5", "sigma", phe$Tissue)
phe$Tissue = gsub("GSM6614365_CD-6", "sigma", phe$Tissue)
sce.all@meta.data=phe
View(phe)


saveRDS(sce.all,"sce.all.rds")
