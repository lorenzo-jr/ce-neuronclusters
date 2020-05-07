library(Seurat)
library(ggplot2)
library(readr)

#try t.test
tryTtest <- function (sample,alternative = "greater") {
  return(tryCatch(t.test(sample,alternative = alternative), error=function(sample) list(p.value=1)))
}


getW<-function(gene,expression){
  if (gene %in% rownames(expression)){
    return(1-sum(expression[gene,]>0)/length(expression[gene,]))
  }
  else{
    return(NA)
  }
}

getGenes<-function(neuron,brainatlas=cured_brain_atlas_PL_and_OH,genelist=names(wheigt)){
  temp=brainatlas[genelist,]
  return(rownames(temp[temp[,neuron]>0,]))
}

### Load needed datasets
cured_brain_atlas_PL_and_OH <- read.table("../Datasets/cured brain atlas PL and OH-CircularOtherMarkers.csv",
                                          sep = ";",
                                          row.names=1,
                                          header = TRUE)


Atlas_simbol <- read.table("../Datasets/Atlas_simbol-Circular.csv",
                           sep = ",",
                           stringsAsFactors=FALSE,
                           header = TRUE)

row.names(Atlas_simbol)<-Atlas_simbol$ID

load("../Datasets/rawdata.Neurons.rda")
load("../Datasets/id_simbol.rda")


id_simbol=data.frame(gene_id=as.character(id_simbol$gene_id),
                     symbol=as.character(id_simbol$symbol), 
                     row.names =row.names(id_simbol),
                     stringsAsFactors = FALSE)


#Change names as in expression atlas
for (geneID in Atlas_simbol$ID) {
  if (geneID%in%rownames(id_simbol)) {
    if (id_simbol[geneID,2]!=Atlas_simbol[geneID,2]) {
      print(paste(geneID, id_simbol[geneID,2],Atlas_simbol[geneID,2],sep = " "))   
      id_simbol[geneID,2]=Atlas_simbol[geneID,2]
    }
    
  }
  else{
    print(paste("not in data",geneID, Atlas_simbol[geneID,2],sep = " "))
  }
}

load("CustomIterationClusters3dit.rda")

#Refining last assignements
CustomIterationClusters = CustomIterationClustersThird[!CustomIterationClustersThird == "27.NA"]# to discard 27.4: only 6 cells
CustomIterationClusters = CustomIterationClusters[!CustomIterationClusters == "6.0"]# other not neuronal: Rectal gland?
CustomIterationClusters = CustomIterationClusters[!CustomIterationClusters == "51.1"]# other not neuronal: germline

dfFinalAssignement <- data.frame(cells = colnames(rawdata.neuronsALL),
              cluster="NA",stringsAsFactors = FALSE)
rownames(dfFinalAssignement)<-dfFinalAssignement$cells
dfFinalAssignement[names(CustomIterationClusters),"cluster"]=as.character(CustomIterationClusters)

#Removing doublets previously analyzed
load("../DoubletAnalysis/doublets.rda")
FinalAssignement <- factor(dfFinalAssignement$cluster)
names(FinalAssignement) <- dfFinalAssignement$cells

FinalAssignement <- FinalAssignement[!names(FinalAssignement) %in% doublets]
rawdata.neuronsALL<-rawdata.neuronsALL[,names(FinalAssignement)]
row.names(rawdata.neuronsALL)<-id_simbol[,2]


#CREATING SEURAT OBEJECT
CeNeurons <- CreateSeuratObject(counts = rawdata.neuronsALL, min.cells = 0, min.features = 0, project = "Cao Neurons")

CeNeurons <- NormalizeData(CeNeurons, normalization.method = "LogNormalize",scale.factor = 10000)

CeNeurons <- FindVariableFeatures(CeNeurons, selection.method = "vst", nfeatures = 10000)

all.genes <- rownames(CeNeurons)
CeNeurons <- ScaleData(CeNeurons, features = all.genes,do.scale = TRUE,do.center = TRUE)

CeNeurons <- RunPCA(CeNeurons, npcs = 100)


#Stouffer's method calculations
CeNeurons.assay=GetAssay(CeNeurons)
wheigt=sapply(row.names(cured_brain_atlas_PL_and_OH), getW,CeNeurons.assay@data)
wheigt=na.omit(wheigt)
scaled=CeNeurons.assay@scale.data[row.names(CeNeurons.assay@scale.data)%in% names(wheigt),]
wheigt=wheigt[rownames(scaled)]
sqWheigtNeuron=numeric()
for (neuron in colnames(cured_brain_atlas_PL_and_OH)) {
  NeuronalIDgenes=c(getGenes(neuron))
  sqWheigtNeuron[neuron]=sum(wheigt[NeuronalIDgenes]*wheigt[NeuronalIDgenes])
}
sqWheigtNeuron=sqrt(sqWheigtNeuron)

scaled=scaled*wheigt

#Best N of PCs dimensions from first screening
dims=92
#Best resolution parameter from first screening
SeuratResol=4
#Create PDF with all plots
pdf(paste("./BoxPlotClusterMedianZscore_nPC_",paste(dims,"_Resol_",SeuratResol),".pdf")
    ,paper="A4",
    height = 0.01, 
    width = 0.05, 
    useDingbats=FALSE)
par(las=2)
par(mar = c(6,5,2,1))
par(mfrow=c(4,1))
par(cex=0.3)

#Constructs a Shared Nearest Neighbor (SNN) Graph
CeNeurons <- FindNeighbors(CeNeurons, reduction = "pca", dims = 1:dims)
# Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering.
# Using original Louvain 
CeNeurons <- FindClusters(CeNeurons, resolution = SeuratResol) 
# Run t-SNE dimensionality reduction on selected features. 
CeNeurons <- RunTSNE(CeNeurons, dims = 1:dims, tsne.method = "Rtsne", nthreads = 4, max_iter = 2000)

# load Final Assignement
Idents(CeNeurons)<-FinalAssignement
CeNeurons$seurat_clusters<-FinalAssignement


p1= DimPlot(CeNeurons, reduction = "tsne", pt.size = 0.5,label = TRUE,order = (length(levels(CeNeurons$seurat_clusters))-1):0) + 
  ggtitle(label = paste("Selected PCs:",dims,". Resolution:",SeuratResol)) +
  theme(legend.position="bottom")

p2= ElbowPlot(CeNeurons, ndims = length(CeNeurons[["pca"]]))+ 
  geom_vline(xintercept=dims, colour = "red")

# plot t-SNE and Elbow
plot(p1)
plot(p2)

# This calculates a Stouffer combined z-score for each cell using
# scaled values of gen expression corresponding to the different 
# neuronal identities
databoxplot<-data.frame(row.names = names(CeNeurons$seurat_clusters))
for (neuron in colnames(cured_brain_atlas_PL_and_OH)) {
  NeuronTest=scaled[c(getGenes(neuron)),]
  NeuronTest=apply(NeuronTest, 2, sum)
  NeuronTest=NeuronTest/sqWheigtNeuron[neuron]#Stouffer's method calculations
  NeuronTest=NeuronTest[rownames(databoxplot)]
  databoxplot<-cbind(databoxplot,NeuronTest)
  names(databoxplot)[names(databoxplot) == "NeuronTest"]<-neuron
}

# plot all boxplots corresponding to obtained clusters
for (cluster in levels(CeNeurons$seurat_clusters)) {
  NeuronOrder=levels(reorder(colnames(databoxplot), 
                             apply(databoxplot[CeNeurons$seurat_clusters==cluster,], 2, median)))
  
  data2test=databoxplot[CeNeurons$seurat_clusters==cluster,rev(NeuronOrder)]
  
  
  ptreshold=apply(data2test, 2, FUN = tryTtest,alternative = "greater")
  ptreshold <- sapply(ptreshold, '[[', 'p.value')
  pvalues=rep(1,length(NeuronOrder))
  for (i in 1:(length(NeuronOrder)-1)) {
    test=try(t.test(data2test[i],data2test[i+1],alternative = "greater"),silent = TRUE)#manage t.test error
    if (!is.list(test)) {
      break
    }
    pvalues[i]=test$p.value
    if (test$p.value<0.05) {
      break
    }
  }
  
  FDRtreshold=p.adjust(ptreshold,method = "fdr")
  names(pvalues)<-rev(NeuronOrder)
 
  
  #plot boxplot

  boxplot(data2test,
          ylim = c(-2, 22),
          outline = FALSE,
          notch = TRUE,
          ylab="Combined z-score",col = c(c("white","red")[((pvalues<1)*(FDRtreshold<0.05))+1]))
  
  legend("topleft", paste("cluster ",cluster,". nº cells: ",sum(CeNeurons$seurat_clusters==cluster),". nPC: ",dims,". Resolution: ",SeuratResol), bty="n",cex = 2)
  
  
  
  
}
dev.off()

# #save second iteration clusters
# CeNeurons_2st_seurat_clusters=CeNeurons$seurat_clusters
# save(file="7-FinalPlotAndClusteringMarkers/CeNeurons_2st_seurat_clusters.rda",CeNeurons_2st_seurat_clusters)


# Write final expression table
AverageExpression92.4second=AverageExpression(CeNeurons)
write.csv(AverageExpression92.4second,file ="./AverageExpression92.4.csv")

# Write final markers table
AllMarkers92.4second=FindAllMarkers(CeNeurons)
write.csv(AllMarkers92.4second,file ="./AllMarkers92.4.csv")

####################################################
# Add anotation data to Markers table
library(readxl)
CLUSTER_ID <- read_excel("../Datasets/ClusteringfinalAssignment.xlsx")

rownames(id_simbol)<-id_simbol$symbol

AllMarkers92_4 <- read_csv("./AllMarkers92.4.csv",
                                 col_types = cols(cluster = col_character()))

AllMarkers92_4[is.na(AllMarkers92_4$cluster),]$cluster <- "Non-neuronal"

uniqueMarkers=AllMarkers92_4
uniqueMarkers=uniqueMarkers[uniqueMarkers$avg_logFC>0,]
duplicatedMarkers=table(uniqueMarkers$gene)
duplicatedMarkers=names(duplicatedMarkers[duplicatedMarkers>1])

uniqueMarkers=uniqueMarkers[!uniqueMarkers$gene %in% duplicatedMarkers,]

AllMarkers92_4=cbind(AllMarkers92_4,clusterID="NA",stringsAsFactors = FALSE)
AllMarkers92_4=cbind(AllMarkers92_4,Unique="",stringsAsFactors = FALSE)
AllMarkers92_4$Unique[AllMarkers92_4$gene %in% uniqueMarkers$gene]="Unique"

AllMarkers92_4[AllMarkers92_4$cluster=="Non-neuronal",]$clusterID <- "Non-neuronal"

CLUSTER_ID=cbind(CLUSTER_ID,paste(CLUSTER_ID$Cluster,CLUSTER_ID$`Manual assignment to Neuron classes`),stringsAsFactors = FALSE)

for (ID in CLUSTER_ID$`Cluster`) {
  AllMarkers92_4[AllMarkers92_4$cluster==ID,9]=as.character(CLUSTER_ID[CLUSTER_ID$`Cluster`==ID,4])
}

for (symbol in id_simbol$symbol) {
  AllMarkers92_4[AllMarkers92_4$gene==symbol,"X1"]=as.character(id_simbol[symbol,]$gene_id)
}
AllMarkers92_4[AllMarkers92_4$gene=="nlp-49","X1"]="WBGene00019160"
AllMarkers92_4[AllMarkers92_4$gene=="nlp-56","X1"]="WBGene00013334"
AllMarkers92_4[AllMarkers92_4$gene=="Y111B2A.8","X1"]="WBGene00013732"

colnames(AllMarkers92_4) <- c("WB_id","p_val","avg_logFC","pct.1","pct.2","p_val_adj","cluster","gene","clusterID","Unique_Marquer")
write.csv(AllMarkers92_4,"./AllMarkers92.4WBid.csv")


#average counts per cluster
data <- data.frame()
for (cluster in levels(CeNeurons$seurat_clusters)) {
  cells <- CeNeurons$seurat_clusters %in% cluster
  if (sum(cells) >0) {
    data <- rbind(data, 
                  data.frame(cluster = cluster, 
                             average_RNAcount = sum(CeNeurons$nCount_RNA[cells])/sum(cells),
                             average_GENEcount = sum(CeNeurons$nFeature_RNA[cells])/sum(cells))
    )
  }
}
write.csv(data,"./AverageRNAcount.csv")

