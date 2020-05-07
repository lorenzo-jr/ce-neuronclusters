
load("neuronsCao.rda")

load("../../6-PlotBestThirdSubClustering/CustomIterationClusters3dit.rda")
#Refining last assignements
CustomIterationClustersThird = CustomIterationClustersThird[!CustomIterationClustersThird == "27.NA"]# to discard 27.4: only 6 cells
CustomIterationClustersThird = CustomIterationClustersThird[!CustomIterationClustersThird == "6.0"]# other not neuronal: Rectal gland?
CustomIterationClustersThird = CustomIterationClustersThird[!CustomIterationClustersThird == "51.1"]# other not neuronal: germline

LorenzoClusters <- as.data.frame(CustomIterationClustersThird)
LorenzoClusters <- cbind(LorenzoClusters,cell = rownames(LorenzoClusters))

mat <- matrix(nrow = length(unique(LorenzoClusters$CustomIterationClustersThird)))
for (clust1 in unique(neuronsCao$cluster)) {
  clusCao <- neuronsCao[neuronsCao$cluster == clust1,]$cell
  colcomp <- c()
  for (clust2 in unique(LorenzoClusters$CustomIterationClustersThird)) {
    
    clusLorenzo <- LorenzoClusters[LorenzoClusters$CustomIterationClustersThird == clust2,]$cell
    colcomp <- c(colcomp,length(intersect(clusCao,clusLorenzo)))
  }
  mat <- cbind(mat,colcomp) 
  
}
rownames(mat) <- unique(LorenzoClusters$CustomIterationClustersThird)
colnames(mat) <- c("na",unique(neuronsCao$cluster))
write.csv(mat,file = "CellOveralpMatrx.csv")


#new cells come from Cao all tissues cds clusters: 5  8 13 15 17 21
# but mostly all cells were finally labeled as non-neuronal
difcells <- setdiff(names(CustomIterationClustersThird),neuronsCao$cell)

print(CustomIterationClustersThird[names(CustomIterationClustersThird) %in% difcells])
