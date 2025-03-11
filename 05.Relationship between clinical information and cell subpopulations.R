rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(plotrix)
library(ggsci)
library(celldex)
library(singleseqgset)
library(devtools)

setwd("D:\\CD+Z\\03.单细胞分析\\05")
sce.all=readRDS("6pbmc_SingleR.rds")
sce.all

head(sce.all@meta.data)
table(sce.all$SingleR)
mynames<-table(sce.all$SingleR) %>% names()
myratio <- table(sce.all$SingleR) %>% as.numeric()
pielabel <- paste0(mynames," (", round(myratio/sum(myratio)*100,2), "%)")

cols <-c('#E64A35','#4DBBD4','#01A187','#6BD66B','#3C5588','#F29F80',
         '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF','#CE3D33',
         '#739B57','#EFE685','#446983','#BB6239','#5DB1DC','#7F2268','#800202','#D8D8CD'
)
pie(myratio, labels=pielabel,
    radius = 1.0,clockwise=T,
    main = "celltype",col = cols)


library(tidyr)
library(reshape2)
tb=table(sce.all$Gender, sce.all$SingleR)
head(tb)
library (gplots)
balloonplot(tb)
bar_data <- as.data.frame(tb)
bar_per <- bar_data %>%
  group_by(Var1) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per)
#write.csv(bar_per,file = "celltype_by_group_percent.csv")
col=c("#C0E2FD","#FEC0C1","#CDC6FF","#FDC0F7","#F3D8F1",    "#D6EBBF","#E1CAF7","#BFDCE2","#F8F0BE","#BEEFBF","#F8C9C8","#C0E2D2","#E9BFCD","#E3E3E3")
colnames(bar_per)

library(ggthemes)
p1 = ggplot(bar_per, aes(x = percent, y = Var1)) +
  geom_bar(aes(fill = Var2) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "top",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = " ", fill = NULL)+labs(x = 'Relative proportion(%)')+
  scale_fill_manual(values=col)+
  theme_few()+
  theme(plot.title = element_text(size=12,hjust=0.5))

p1
