
library(GOplot)
setwd("")
ego=read.table("GO.txt", header = T,sep="\t",check.names=F)    
go=data.frame(Category = "All",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)

id.fc <- read.table("id.txt", header = T,sep="\t",check.names=F)
genelist <- data.frame(ID = id.fc$id, logFC = id.fc$logFC)
row.names(genelist)=genelist[,1]

circ <- circle_dat(go, genelist)
termNum = 3                                  
geneNum = nrow(genelist)            

chord <- chord_dat(circ, genelist[1:geneNum,], go$Term[1:termNum])
pdf(file="circ.pdf",width = 13,height = 13)
GOChord(chord, 
        space = 0.001,         
        gene.order = 'logFC',   
        gene.space = 0.25,     
        gene.size = 5,         
        border.size = 0.1,     
        process.label = 10.5) 
dev.off()

termCol <- c("#223D6C","#D20A13","#FFD121","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")
pdf(file="cluster.pdf",width = 13,height = 13)
GOCluster(circ.gsym, 
          go$Term[1:termNum], 
          lfc.space = 0.2,                  
          lfc.width = 1,                
          term.col = termCol[1:termNum],    
          term.space = 0.2,             
          term.width = 1)                  
dev.off()          
