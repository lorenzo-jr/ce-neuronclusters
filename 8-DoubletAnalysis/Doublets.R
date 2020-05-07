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
Atlas_simbol <- read.table("../Datasets/Atlas_simbol-Circular.csv",
                           sep = ",",
                           stringsAsFactors=FALSE,
                           header = TRUE)
row.names(Atlas_simbol) <- Atlas_simbol$ID

load("../Datasets/rawdata.Neurons.rda")
load("../Datasets/id_simbol.rda")


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

NeuronalClustering <- CreateSeuratObject(counts = NeuronalClustering, min.cells = 0, min.features = 0, project = "Recluster Cao et al. Neurons")

NeuronalClustering <- NormalizeData(NeuronalClustering, normalization.method = "LogNormalize",scale.factor = 10000)

NeuronalClustering <- FindVariableFeatures(NeuronalClustering, selection.method = "vst", nfeatures = 10000)

all.genes <- rownames(NeuronalClustering)
NeuronalClustering <- ScaleData(NeuronalClustering, features = all.genes,do.scale = TRUE,do.center = TRUE)

NeuronalClustering <- RunPCA(NeuronalClustering, npcs = 100)

NeuronalClustering <- FindNeighbors(NeuronalClustering, reduction = "pca", dims = 1:PCs)
NeuronalClustering <- FindClusters(NeuronalClustering, resolution = Resols)  
NeuronalClustering <- RunTSNE(NeuronalClustering, dims = 1:PCs, tsne.method = "Rtsne", nthreads = 4, max_iter = 2000)



## DoubletFinder script

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep_v3(NeuronalClustering, PCs = 1:92, sct = FALSE,num.cores = 4)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn$pK<-as.numeric(as.character(bcmvn$pK))
pK1=bcmvn$pK[bcmvn$BCmetric==max(bcmvn$BCmetric)]

p1=ggplot(data=bcmvn, aes(x=pK, y=BCmetric, group=2)) +
  geom_line(color="blue")+
  geom_point()+
  geom_vline(xintercept=pK1, linetype="dashed", color = "red")+
  labs(title="pK Selection",x="pK", y = "BCmvn")+
  theme_classic()
  

################## 3.1
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <-NeuronalClustering@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.031*length(NeuronalClustering$nFeature_RNA))  ## Assuming 3.1% doublet formation rate * number of cells (length(NeuronalClustering$nFeature_RNA) == 9078)
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) #nExp adjusted according to the estimated proportion of homotypic doublets

## Run DoubletFinder ----------------------------------------------------------------
FindDoublet1 <- doubletFinder_v3(NeuronalClustering, PCs = 1:92, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
FindDoublet <- doubletFinder_v3(NeuronalClustering, PCs = 1:92, pN = 0.25, pK = 0.3, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)


DF.classifications=FindDoublet@meta.data[7]

doublets=row.names(DF.classifications)[DF.classifications=="Doublet"]
Singlet=row.names(DF.classifications)[DF.classifications=="Singlet"]
FindDoublet <- SetIdent(FindDoublet, cells = doublets, value = "Doublets")
FindDoublet <- SetIdent(FindDoublet, cells = Singlet, value = "Singlet")

p3= DimPlot(FindDoublet, reduction = "tsne", pt.size = 1) + 
  ggtitle(label = "Assuming 3.1% doublet formation rate")

save(doublets,file = "doublets.rda")
pdf("DoubletPlot.pdf",paper = "a4",height = 10,width = 7)
multiplot(p1,p2)
dev.off()

#Prepare to Write Doublets.csv Table
load("../Datasets/rawdata.Neurons.rda")
load("../4-PlotBestSecondSubClustering/CustomIterationClusters.rda") #load cluster identities
#Refining last assignements
CustomIterationClusters = CustomIterationClusters[!CustomIterationClusters == "27.NA"]# to discard 27.4: only 6 cells
CustomIterationClusters = CustomIterationClusters[!CustomIterationClusters == "6.0"]# other not neuronal: Rectal gland?
CustomIterationClusters = CustomIterationClusters[!CustomIterationClusters == "51.1"]# other not neuronal: germline
dfFinalAssignement <- data.frame(cells = colnames(rawdata.neuronsALL),
                                 cluster="NA",stringsAsFactors = FALSE)
rownames(dfFinalAssignement)<-dfFinalAssignement$cells
dfFinalAssignement[names(CustomIterationClusters),"cluster"]=as.character(CustomIterationClusters)
FinalAssignement <- factor(dfFinalAssignement$cluster)
names(FinalAssignement) <- dfFinalAssignement$cells

# Write results in a table
write.csv(cbind(doublets=table(FinalAssignement[doublets]),total_cells=table(FinalAssignement)),"Doublets.csv")

#Plot comparing Neuronal vs Non-Neuronal
Idents(NeuronalClustering, cells = names(NeuronalClustering$seurat_clusters))<-"Non-Neuronal" 
Idents(NeuronalClustering, cells = names(CustomIterationClusters))<-"Neuronal" 
Idents(NeuronalClustering, cells = doublets)<-"Doublets" 

pdf("DoubletPlotNeuronal.pdf",paper = "a4",height = 5,width = 7)
TSNEPlot(NeuronalClustering, pt.size =1, do.return = T) 
dev.off()

