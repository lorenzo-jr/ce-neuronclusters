library(Seurat)
library(ggplot2)
library(DoubletFinder)

#generates multiple plots
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#### load datasets
Atlas_simbol <- read.table("../../Datasets/Atlas_simbol-Circular.csv",
                           sep = ",",
                           stringsAsFactors=FALSE,
                           header = TRUE)
row.names(Atlas_simbol) <- Atlas_simbol$ID

load("../../Datasets/rawdata.Neurons.rda")
load("../../Datasets/id_simbol.rda")


id_simbol = data.frame(gene_id = as.character(id_simbol$gene_id),
                     symbol = as.character(id_simbol$symbol), 
                     row.names = row.names(id_simbol),
                     stringsAsFactors = FALSE)

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
#### end load datasets

#### Creating Seurat Object
PCs=92
Resols=4
NeuronalClustering <- rawdata.neuronsALL
row.names(NeuronalClustering)<-id_simbol[,2]
load("../doublets.rda")
NeuronalClustering <- NeuronalClustering[,!colnames(NeuronalClustering) %in% doublets]


NeuronalClustering <- CreateSeuratObject(counts = NeuronalClustering, min.cells = 0, min.features = 0, project = "Recluster Cao et al. Neurons")

NeuronalClustering <- NormalizeData(NeuronalClustering, normalization.method = "LogNormalize",scale.factor = 10000)

NeuronalClustering <- FindVariableFeatures(NeuronalClustering, selection.method = "vst", nfeatures = 10000)

all.genes <- rownames(NeuronalClustering)
NeuronalClustering <- ScaleData(NeuronalClustering, features = all.genes,do.scale = TRUE,do.center = TRUE)

NeuronalClustering <- RunPCA(NeuronalClustering, npcs = 100)

NeuronalClustering <- FindNeighbors(NeuronalClustering, reduction = "pca", dims = 1:PCs)
NeuronalClustering <- FindClusters(NeuronalClustering, resolution = Resols)  
NeuronalClustering <- RunTSNE(NeuronalClustering, dims = 1:PCs, tsne.method = "Rtsne", nthreads = 4, max_iter = 2000)

clustersNoDoublets <- NeuronalClustering$seurat_clusters

#clusters withDoublets
#### Creating Seurat Object
PCs=92
Resols=4
NeuronalClustering <- rawdata.neuronsALL
row.names(NeuronalClustering)<-id_simbol[,2]

NeuronalClustering <- CreateSeuratObject(counts = NeuronalClustering, min.cells = 0, min.features = 0, project = "Recluster Cao et al. Neurons")

NeuronalClustering <- NormalizeData(NeuronalClustering, normalization.method = "LogNormalize",scale.factor = 10000)

NeuronalClustering <- FindVariableFeatures(NeuronalClustering, selection.method = "vst", nfeatures = 10000)

all.genes <- rownames(NeuronalClustering)
NeuronalClustering <- ScaleData(NeuronalClustering, features = all.genes,do.scale = TRUE,do.center = TRUE)

NeuronalClustering <- RunPCA(NeuronalClustering, npcs = 100)

NeuronalClustering <- FindNeighbors(NeuronalClustering, reduction = "pca", dims = 1:PCs)
NeuronalClustering <- FindClusters(NeuronalClustering, resolution = Resols)  
NeuronalClustering <- RunTSNE(NeuronalClustering, dims = 1:PCs, tsne.method = "Rtsne", nthreads = 4, max_iter = 2000)

clustersWithDoublets <- NeuronalClustering$seurat_clusters

clustersWithDoublets <- clustersWithDoublets[!colnames(NeuronalClustering) %in% doublets]

#compare clusterings
CompareMatrix <- matrix(data = NA,nrow = length(levels(clustersWithDoublets)), ncol = length(levels(clustersNoDoublets)))
rownames(CompareMatrix)<-levels(clustersWithDoublets)
colnames(CompareMatrix)<-levels(clustersNoDoublets)
for (iDoub in levels(clustersWithDoublets)) {
  iDoubNames = names(clustersWithDoublets[clustersWithDoublets == iDoub])
  for (iNoDoub in levels(clustersNoDoublets)) {
    iNoDoubNames = names(clustersNoDoublets[clustersNoDoublets == iNoDoub])
    testUnion = length(union(iDoubNames,iNoDoubNames))
    testIntersect = length(intersect(iDoubNames,iNoDoubNames))
    # print(testUnion)
    # print(testIntersect/testUnion)
    CompareMatrix[iDoub,iNoDoub] <- testIntersect/testUnion
  }
  
}
write.csv(CompareMatrix,file = "CompareMatrix.csv")


