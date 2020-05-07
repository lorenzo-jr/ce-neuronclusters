library(ggplot2)
source("functions.R")
library(readr)

# aax1971_Table_S14 <- read_delim("aax1971_Table_S14.csv",
# "\t", escape_double = FALSE, trim_ws = TRUE)
# CellTypes=read.csv("Neurons",sep = "\t",header = F,stringsAsFactors = F)
# CellTypes=CellTypes$V1

load("ct_tbl")
IDsFinal <- read_csv("IDsFinal.csv")
IDsFinal <- na.omit(as.character(IDsFinal$x))
aax1971_Table_S14=ct_tbl[ct_tbl$cell.bin %in% IDsFinal,]
CellTypes<-IDsFinal

aax1971_Table_S14=aax1971_Table_S14[aax1971_Table_S14$adjusted.tpm.estimate>0,]
aax1971_Table_S14=aax1971_Table_S14[order(aax1971_Table_S14$raw.tpm.estimate,decreasing = T),]

AllMarkers92_4second <- read_csv("AllMarkers92.4WBid.csv")
AllMarkers92_4second<-as.data.frame(AllMarkers92_4second)
AllMarkers92_4second<-AllMarkers92_4second[AllMarkers92_4second$avg_logFC>0,]
genesdetected=length(unique(c((as.character(AllMarkers92_4second$WB_id)),as.character(aax1971_Table_S14$gene.id))))

# aax1971_Table_S14=aax1971_Table_S14[aax1971_Table_S14$cell.type=="FLP_PVD",]

matrix=matrix(nrow = length(CellTypes), ncol =  length(unique(AllMarkers92_4second$cluster)))
colnames(matrix)<-unique(AllMarkers92_4second$cluster)
rownames(matrix)<-CellTypes

n2overlap=200
# AllMarkers92_4second<-AllMarkers92_4second[AllMarkers92_4second$p_val_adj<0.05,]
for (CellT in CellTypes) {
  temp=aax1971_Table_S14[aax1971_Table_S14$cell.bin==CellT,]
  results=data.frame()
  for (Cluster in unique(AllMarkers92_4second$cluster)) {
    # n2overlap=length(AllMarkers92_4second[AllMarkers92_4second$cluster==Cluster,"X1"])
    overlap=length(intersect(AllMarkers92_4second[AllMarkers92_4second$cluster==Cluster,]$WB_id[1:n2overlap],
                             temp$gene.id[1:n2overlap]))
    pval=overlap.pvalue(n2overlap,n2overlap,overlap,N=genesdetected)
    results=rbind(results,
                  data.frame(Cluster=Cluster,
                             AssignedIDs=AllMarkers92_4second[AllMarkers92_4second$cluster==Cluster,]$clusterID[1],
                             Overlap=overlap,
                             pvalue=pval))
  }
  results=cbind(results,data.frame(
    FDR=p.adjust(results$pvalue,method = "fdr"),
    Bonferroni=p.adjust(results$pvalue,method = "bonferroni")
  ))
  rownames(results)<-results$Cluster
  matrix[CellT,]=results[colnames(matrix),"Bonferroni"]
  
}
# colnames(matrix)<-paste(results[colnames(matrix),]$Cluster,results[colnames(matrix),]$AssignedIDs)

write.csv(matrix,"resultMatrix.csv")


results=results[order(results$FDR,decreasing = F),]
