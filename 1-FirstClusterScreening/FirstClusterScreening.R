#Load functions and datasets needed
source("InitializationScript.R")

####################################################################################

## PARAMETERS TO SCREEN IN SEURAT CLUSTERING (PCA DIMENSIONS AND RESOLUTION)  
PCs = 40:100
Resols = 1:8

IterationData = data.frame()

# iterate changing number of dimensions and resolution parameters
for (dims in PCs) {
  for (SeuratResol in Resols) {
    
    #Constructs a Shared Nearest Neighbor (SNN) Graph
    clusterScreen = FindNeighbors(CeNeurons, reduction = "pca", dims = 1:dims)
    
    # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering.
    # Using original Louvain 
    clusterScreen = FindClusters(clusterScreen, resolution = SeuratResol)  
    
    # This calculates a Stouffer combined z-score for each cell using
    # scaled values of gen expression corresponding to the different 
    # neuronal identities
    databoxplot = data.frame(row.names = names(clusterScreen$seurat_clusters))
    for (neuron in colnames(cured_brain_atlas_PL_and_OH)) { #for neuron in expression atlas 
      NeuronTest = scaled[c(getGenes(neuron)),] #get expression of genes corresponding to each neuron identity
      NeuronTest = apply(NeuronTest, 2, sum)
      NeuronTest = NeuronTest / sqWheigtNeuron[neuron]# Second part of Stouffer's method calculation
      NeuronTest = NeuronTest[rownames(databoxplot)]
      databoxplot = cbind(databoxplot,NeuronTest) #attach the calculated neuron score to databoxplot
      names(databoxplot)[names(databoxplot) == "NeuronTest"] = neuron
    }
    
    # number of clusters produced in the iteration
    n_clusters = length(levels(clusterScreen$seurat_clusters))
    
    for (cluster in levels(clusterScreen$seurat_clusters)) {
      cellsInCluster = names(clusterScreen$seurat_clusters[clusterScreen$seurat_clusters == cluster])
      NeuronOrder = levels(reorder(colnames(databoxplot), apply(databoxplot[cellsInCluster,], 2, median)))
      data2test = databoxplot[cellsInCluster,rev(NeuronOrder)]
      
      
      # A sequential t-test from the highest to lowest ranked neuron classes
      # applied iteratively to the successive pair
      pvalues = rep(1, length(NeuronOrder))
      names(pvalues) = colnames(data2test)
      for (i in 1:(length(NeuronOrder) - 1)) {
        # Try two sample t-test. 
        test = try(t.test(data2test[i], data2test[i+1], alternative = "greater"), silent = TRUE)#manage t.test error
        
        if (!is.list(test)) {# If it fails set p-value = 1 (this happend in small clusters when cells have the same score)
          break
        }
        pvalues[i] = test$p.value
        if (test$p.value < 0.05) {
          break
        }
      }
      
      # one sample t-test with mu=0 (in the scaled data zero is the mean)
      ptreshold = apply(data2test, 2, FUN = tryTtest,alternative = "greater")
      ptreshold = sapply(ptreshold, '[[', 'p.value')
      # FDR adjustment of the p-value
      FDRtreshold = p.adjust(ptreshold,method = "fdr")
      
      # Name of automatically detected identities
      totalIdentities = intersect(names(FDRtreshold[FDRtreshold < 0.05]),names(pvalues[pvalues < 1]))
      cluster_id = paste(totalIdentities,collapse = ", ")
      
      # number of identities atumatically detected
      n_ident = length(totalIdentities)
      
      # number of cells in cluster
      CellNumber = length(cellsInCluster)
      
      # In each iteration for each cluster fill the iteration data
      IterationData = rbind(IterationData,
                           data.frame(cluster = cluster,# cluster id
                                      dim = dims, # number of dimensions used
                                      resol = SeuratResol, # resolution used
                                      n_cells = CellNumber, # number of cells in cluster
                                      n_ids = n_ident, # number of identities atumatically detected
                                      cluster_ids = cluster_id, # list of identities atumatically detected
                                      n_clusters = n_clusters # number of clusters produced in the iteration
                                      ))
      }
    }
  }



### WRITE RESULTS
duplicated = matrix(nrow = length(1:100), ncol = length(1:8)) #Number of clusters sharing identities 
nclusters = matrix(nrow = length(1:100), ncol = length(1:8)) #Total number of clusters
n1idclusters = matrix(nrow = length(1:100), ncol = length(1:8)) #Number of clusters with only 1 identity 

oneIdClusters = NA
duplicatedIds = NA
IterationData = cbind(IterationData, oneIdClusters)
IterationData = cbind(IterationData, duplicatedIds)
for (dims in PCs) {
  for (SeuratResol in Resols){
    Test = IterationData[IterationData$dim == dims & IterationData$resol == SeuratResol,]
    oneIdClusters = length(Test[Test$n_ids == 1,]$n_ids)
    duplicatedIds = sum(duplicated(Test[Test$n_ids == 1,]$cluster_ids))
    IterationData[IterationData$dim == dims & IterationData$resol == SeuratResol,"oneIdClusters"] = oneIdClusters
    IterationData[IterationData$dim == dims & IterationData$resol == SeuratResol,"duplicatedIds"] = duplicatedIds
    duplicated[dims,SeuratResol] = duplicatedIds
    nclusters[dims,SeuratResol] = IterationData[IterationData$dim == dims & IterationData$resol == SeuratResol,"n_clusters"][1]
    n1idclusters[dims,SeuratResol] = oneIdClusters
  }
}

write.table(duplicated[40:100,], file = "duplicated,csv",row.names = paste("PCs",40:100,sep = "_"), col.names = paste("resol",1:8,sep = "_"))
write.table(nclusters[40:100,], file = "nclusters.csv",row.names = paste("PCs",40:100,sep = "_"), col.names = paste("resol",1:8,sep = "_"))
write.table(n1idclusters[40:100,], file = "n1idclusters.csv",row.names = paste("PCs",40:100,sep = "_"), col.names = paste("resol",1:8,sep = "_"))

