source("InitializationScript.R")
load("../4-PlotBestSecondSubClustering/CustomIterationClusters.rda")

#Load parameters to assign clusters
ParameterChoose <- read_delim("../Datasets/ParameterChoose3iteration.csv",
                              ";", escape_double = FALSE, trim_ws = TRUE,col_types = "cddccc")

#Get parameters to assign clusters
ParametersRecluster <- ParameterChoose[!ParameterChoose$dim == -1,1:3]

colnames(ParametersRecluster) <- c("parent_cluster","dim","resol")

ClusToTest <- ParametersRecluster$parent_cluster

databoxplot <- data.frame(row.names = names(CustomIterationClusters))
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
  cells = names(CustomIterationClusters[CustomIterationClusters == parentcluster])
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
#merged names
FinalIDs <- paste(FinalSubclusterIDs[,2],FinalSubclusterIDs[,6], sep = ".")

# Merge 10.0.1 (AVH?) with 10.0.4
FinalIDs[FinalIDs == "10.0.1"] <- "10.0.1_4" 
FinalIDs[FinalIDs == "10.0.4"] <- "10.0.1_4" 

# Merge 10.0.2 (BDU?) with 10.0.3
FinalIDs[FinalIDs == "10.0.2"] <- "10.0.2_3" 
FinalIDs[FinalIDs == "10.0.3"] <- "10.0.2_3" 

FinalSubclusterIDs=cbind(FinalSubclusterIDs,FinalIDs)
# combine and save data for final assignement of clusters
CustomIterationClustersThird <- factor(FinalSubclusterIDs[,7])
names(CustomIterationClustersThird) <- FinalSubclusterIDs[,1]

clusterID <- CustomIterationClusters[!(names(CustomIterationClusters) %in% names(CustomIterationClustersThird))] #Previos iteration cluster
combine <- as.data.frame(clusterID)
clusterID <- as.data.frame(CustomIterationClustersThird) #Third iteration clusters
colnames(clusterID) <- c("clusterID")
combine <- rbind(combine,clusterID)
CustomIterationClustersThird <- combine$clusterID
names(CustomIterationClustersThird) <- row.names(combine)

ClusterPseudotime <- CustomIterationClustersThird #for psudotime
#merge 16.0 and 16.1 both are FLP
CustomIterationClustersThird[CustomIterationClustersThird == "16.1"] <- "16.0"
#merge 0 and 2
levels(CustomIterationClustersThird)[1]<-"0_2"
CustomIterationClustersThird[CustomIterationClustersThird=="0.1"] <- "0_2"
CustomIterationClustersThird[CustomIterationClustersThird=="2.0"] <- "0_2"

save(CustomIterationClustersThird,file = "CustomIterationClusters3dit.rda")

#Removing doublets previously analyzed for Pseudotime
load("../DoubletAnalysis/doublets.rda")
ClusterPseudotime <- ClusterPseudotime[!names(ClusterPseudotime) %in% doublets]

#Clusters for Pseudotime
ClusterPseudotime <- ClusterPseudotime[ClusterPseudotime %in% c("3.0",
                                                                "3.1",
                                                                "4.0",
                                                                "4.1",
                                                                "16.0",
                                                                "16.1",
                                                                "17.0",
                                                                "17.1",
                                                                "0.0",
                                                                "0.1",
                                                                "1.0",
                                                                "1.1",
                                                                "2.0",
                                                                "33.0",
                                                                "42.0",
                                                                "42.1",
                                                                "56.0",
                                                                "8.0.0",
                                                                "8.0.1",
                                                                "8.1",
                                                                "8.2",
                                                                "8.3",
                                                                "10.0.0",
                                                                "10.0.1_4",
                                                                "10.0.2_3",
                                                                "10.0.5",
                                                                "10.1.0",
                                                                "10.1.1",
                                                                "10.2",
                                                                "56.0")]
RawCountsPseudotime <- rawdata.neuronsALL[,names(ClusterPseudotime)]
save(ClusterPseudotime,file = "../PseudoTime/ClusterPseudotime.rda")
save(RawCountsPseudotime,file = "../PseudoTime/RawCountsPseudotime.rda")


