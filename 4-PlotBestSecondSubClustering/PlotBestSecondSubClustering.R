source("InitializationScript.R")
load("../2-PlotBestFirstClustering/CeNeurons_1st_seurat_clusters.rda")

#Load parameters to assign clusters
ParameterChoose <- read_delim("../Datasets/ParameterChoose.csv",
                              ";", escape_double = FALSE, trim_ws = TRUE)

#Get parameters to assign clusters
ParametersRecluster <- ParameterChoose[!ParameterChoose$dim == -1,1:3]

colnames(ParametersRecluster) <- c("parent_cluster","dim","resol")

ClusToTest <- ParametersRecluster$parent_cluster

databoxplot <- data.frame(row.names = names(CeNeurons_1st_seurat_clusters))
#neurons z-score
for (neuron in colnames(cured_brain_atlas_PL_and_OH)) {
  NeuronTest = scaled[c(getGenes(neuron)),]
  NeuronTest = apply(NeuronTest, 2, sum)
  NeuronTest = NeuronTest/sqWheigtNeuron[neuron]#Stouffer's method calculations
  NeuronTest = NeuronTest[rownames(databoxplot)]
  databoxplot <- cbind(databoxplot,NeuronTest)
  names(databoxplot)[names(databoxplot) == "NeuronTest"]<-neuron
}
SubclusterIDs = data.frame()

PCs=92 #Maximum PCs to test

# generate individual plots for each sbucluster using slected parameters 
for (parentcluster in ClusToTest[order(as.numeric(ClusToTest))]) {
  lsitaClusterFinal = list()
  cells = names(CeNeurons_1st_seurat_clusters[CeNeurons_1st_seurat_clusters == parentcluster])
  CeNeuronSubcluster <- rawdata.neuronsALL[,cells]
  
  CeNeuronSubcluster <- CreateSeuratObject(counts = CeNeuronSubcluster, min.cells = 0, min.features = 0, project = "Cao Neurons")
  
  CeNeuronSubcluster <- NormalizeData(CeNeuronSubcluster, normalization.method = "LogNormalize",scale.factor = 10000)
  
  CeNeuronSubcluster <- FindVariableFeatures(CeNeuronSubcluster, selection.method = "vst", nfeatures = 10000)
  
  all.genes <- rownames(CeNeuronSubcluster)
  CeNeuronSubcluster <- ScaleData(CeNeuronSubcluster, features = all.genes,do.scale = TRUE,do.center = TRUE)
  
  dim=ParametersRecluster[ParametersRecluster$parent_cluster==parentcluster,]$dim
  
  npcs=length(cells)-1
  if (npcs > PCs) {
    npcs = PCs
  }
  CeNeuronSubcluster <- RunPCA(CeNeuronSubcluster, npcs = npcs)
  
  dims = ParametersRecluster[ParametersRecluster$parent_cluster == parentcluster,]$dim
  SeuratResol = ParametersRecluster[ParametersRecluster$parent_cluster == parentcluster,]$resol
  
  CeNeuronSubcluster <- FindNeighbors(CeNeuronSubcluster, reduction = "pca", dims = 1:dims)
  CeNeuronSubcluster <- FindClusters(CeNeuronSubcluster, resolution = SeuratResol)
  
  SubclusterIDs <- rbind(SubclusterIDs,
                      PDFseurat2(CeNeuronSubcluster,dims = dims,SeuratResol = SeuratResol,clustName = parentcluster, path = "./"))
  
}

#Get parameters to assign clusters
MergeAndExtract <- ParameterChoose[!ParameterChoose$dim == -1,c(1,4,5,6)]

FinalSubclusterIDs = cbind(SubclusterIDs, finalClusters=NA)

#use selected parameters to assign clusters
for (cluster in MergeAndExtract$cluster) {
  if (MergeAndExtract[MergeAndExtract$cluster==cluster,]$OnlyExtract == "0") {
    FinalSubclusterIDs[FinalSubclusterIDs$ParentCluster == cluster,"finalClusters"] = "0"
  } else if (!is.na(MergeAndExtract[MergeAndExtract$cluster==cluster,]$Merge)) {
    
    SubclustertoMerge = unlist(#merge subcluster 
      strsplit(MergeAndExtract[MergeAndExtract$cluster==cluster,]$Merge,
               split = ",")
      )
    SubclustertoExtract = unlist(#only extract 
      strsplit(MergeAndExtract[MergeAndExtract$cluster==cluster,]$OnlyExtract,
               split = ",")
    )
    
    FinalSubclusterIDs[FinalSubclusterIDs$Subcluster %in% SubclustertoMerge & FinalSubclusterIDs$ParentCluster==cluster,
                       "finalClusters"] = paste(SubclustertoMerge,collapse = "_")
    
    for (subcluster in SubclustertoExtract) {#only extract 
      FinalSubclusterIDs[FinalSubclusterIDs$Subcluster == subcluster & FinalSubclusterIDs$ParentCluster==cluster,
                         "finalClusters"] =  subcluster
    }

  } else{
    SubclustertoExtract = unlist(#only extract 
      strsplit(MergeAndExtract[MergeAndExtract$cluster==cluster,]$OnlyExtract,
               split = ",")
    )
    for (subcluster in SubclustertoExtract) {#only extract 
      FinalSubclusterIDs[FinalSubclusterIDs$Subcluster == subcluster & FinalSubclusterIDs$ParentCluster == cluster,
                         "finalClusters"] = subcluster
    }
    
  }
  
  
  
}
FinalSubclusterIDs=cbind(FinalSubclusterIDs,paste(FinalSubclusterIDs[,2],FinalSubclusterIDs[,6], sep = "."))

# save data from final assignement of clusters
CustomIterationClusters = factor(FinalSubclusterIDs[,7])
names(CustomIterationClusters) <- FinalSubclusterIDs[,1]
save(CustomIterationClusters,file = "CustomIterationClusters.rda")

