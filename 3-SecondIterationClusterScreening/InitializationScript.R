library(Seurat)
library(readr)
library(ggplot2)

# Multiple plot function
multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Function to generate PDF with dimensional reduction plots and boxplots
PDFseurat2 = function(SeuratObject=SeuratObject,dims=dims,SeuratResol=SeuratResol,clustName=clustName,colorVector=c(0),path=""){
  pdf(paste(path,"BoxPlotClusterMeanZscore_Cluster_",clustName,"_nPC_",paste(dims,"_Resol_",SeuratResol),".pdf"),
      paper="A4",
      height = 0.01, 
      width = 0.05,
      useDingbats=FALSE)
  par(las=2)
  par(mar = c(6,5,2,1))
  par(mfrow=c(4,1))
  par(cex=0.3)
  
  if (dims>0) {
    SeuratObject = RunUMAP(SeuratObject, dims = 1:dims, min.dist = 0.01)
    
    p1=DimPlot(SeuratObject, reduction = "umap", pt.size = 0.5,label = TRUE, order = (length(levels(SeuratObject$seurat_clusters))-1):0) + 
      ggtitle(label = paste("Parent: cluster",clustName,"Selected PCs:",dims,". Resolution:",SeuratResol)) +
      theme(legend.position="bottom")
    
    p2=ElbowPlot(SeuratObject, ndims = length(SeuratObject[["pca"]]))+ 
      geom_vline(xintercept = dims, colour = "red")
    
    multiplot(p1, p2)
  }
  
  SubIdentities = data.frame()
  clusterOrder = 0:(length(levels(SeuratObject$seurat_clusters)) - 1)
  for (cluster in clusterOrder) {
    cellsInCluster = names(SeuratObject$seurat_clusters[SeuratObject$seurat_clusters == cluster])
    NeuronOrder = levels(reorder(colnames(databoxplot), apply(databoxplot[cellsInCluster,], 2, median)))    ###
    data2test = databoxplot[cellsInCluster, rev(NeuronOrder)]

    ptreshold = apply(data2test, 2, FUN = tryTtest,alternative = "greater")
    ptreshold = sapply(ptreshold, '[[', 'p.value')
    pvalues = rep(1,length(NeuronOrder))
    names(pvalues) = colnames(data2test)
    for (i in 1:(length(NeuronOrder)-1)) {
      test = try(t.test(data2test[i], data2test[i+1], alternative = "greater"),silent = TRUE)#manage t.test error
      if (!is.list(test)) {
        break
      }
      pvalues[i] = test$p.value
      if (test$p.value < 0.05) {
        break
      }
    }
    
    FDRtreshold = p.adjust(ptreshold, method = "fdr")
    
    boxplot(data2test,
            ylim = c(-2, 22),
            outline = FALSE,
            notch = TRUE,
            ylab="Mean z-score", col = c(c("white", "red")[((pvalues < 1 ) * (FDRtreshold < 0.05)) + 1]))
    
    legend("topleft", 
           paste("cluster ",
                 cluster,
                 ". nº cells: ",
                 sum(SeuratObject$seurat_clusters==cluster),
                 ". nPC: ",
                 dims,
                 ". Resolution: ",
                 SeuratResol),
           bty="n",
           cex = 2)
    
    SubIdentities<-rbind(SubIdentities,
                         data.frame(CellsInSubcluster = cellsInCluster,
                                    ParentCluster = clustName,
                                    Subcluster = cluster,
                                    IDsSublcuster = paste(names(FDRtreshold[(pvalues < 1) & (FDRtreshold<0.05)]), collapse = ","),
                                    N_IDs = length(names(FDRtreshold[(pvalues<1) & (FDRtreshold<0.05)]))
                                    
                         ))
    
    
    
  }
  dev.off()
  return(SubIdentities)
  
  
}

#try t.test
tryTtest = function (sample,alternative = "greater") {
  return(tryCatch(t.test(sample, alternative = alternative), error = function(sample) list(p.value = 1)))
}

#calculate the weights
getW = function(gene,expression){
  if (gene %in% rownames(expression)){
    return(1 - sum(expression[gene, ] > 0) / length(expression[gene, ]))
  }
  else{
    return(NA)
  }
  
}

#Function that get genes expresed in a defined neuron as reported in the brain atlas
getGenes = function(neuron, brainatlas = cured_brain_atlas_PL_and_OH, genelist = names(wheigt)){
  temp = brainatlas[genelist,]
  return(rownames(temp[temp[, neuron] > 0, ]))
}

#read brain atlas table
cured_brain_atlas_PL_and_OH = read.table("../Datasets/cured brain atlas PL and OH-Circular.csv",
                                          sep = ";",
                                          row.names=1,
                                          header = TRUE)
#load a table with used gen symbols in the brain atlas table and the equivalent and Wormbase IDs
Atlas_simbol = read.table("../Datasets/Atlas_simbol-Circular.csv",
                           sep = ",",
                           stringsAsFactors=FALSE,
                           header = TRUE)
row.names(Atlas_simbol) = Atlas_simbol$ID

#load single cell gen expression count
load("../Datasets/rawdata.Neurons.rda")

#load a table with used gen symbols in rawdata.Neurons and the equivalent Wormbase IDs
load("../Datasets/id_simbol.rda")
id_simbol = data.frame(gene_id = as.character(id_simbol$gene_id),
                     symbol = as.character(id_simbol$symbol), 
                     row.names = row.names(id_simbol),
                     stringsAsFactors = FALSE)

#This change rawdata.Neurons symbols to brain atlas gen symbols
for (geneID in Atlas_simbol$ID) {
  if (geneID%in%rownames(id_simbol)) {
    if (id_simbol[geneID, 2]!=Atlas_simbol[geneID, 2]) {
      print(paste(geneID, id_simbol[geneID, 2], Atlas_simbol[geneID, 2],sep = " "))   
      id_simbol[geneID, 2] = Atlas_simbol[geneID, 2]
    }
    
  }
  else{
    print(paste("not in data", geneID, Atlas_simbol[geneID, 2], sep = " "))
  }
  
}
row.names(rawdata.neuronsALL) = id_simbol[,2]


#Creates Seurat Object
CeNeurons = CreateSeuratObject(counts = rawdata.neuronsALL, min.cells = 0, min.features = 0, project = "Cao Neurons")

CeNeurons = NormalizeData(CeNeurons, normalization.method = "LogNormalize", scale.factor = 10000)

CeNeurons = FindVariableFeatures(CeNeurons, selection.method = "vst", nfeatures = 10000)

all.genes = rownames(CeNeurons)
CeNeurons = ScaleData(CeNeurons, features = all.genes,do.scale = TRUE,do.center = TRUE)

CeNeurons = RunPCA(CeNeurons, npcs = 100)


#First part of Stouffer's method calculations
CeNeurons.assay = GetAssay(CeNeurons)
wheigt = sapply(row.names(cured_brain_atlas_PL_and_OH), getW,CeNeurons.assay@data)
wheigt = na.omit(wheigt)
scaled = CeNeurons.assay@scale.data[row.names(CeNeurons.assay@scale.data) %in% names(wheigt),]
wheigt = wheigt[rownames(scaled)]
sqWheigtNeuron = numeric()
for (neuron in colnames(cured_brain_atlas_PL_and_OH)) {
  NeuronalIDgenes = c(getGenes(neuron))
  sqWheigtNeuron[neuron] = sum(wheigt[NeuronalIDgenes] * wheigt[NeuronalIDgenes])
}
sqWheigtNeuron = sqrt(sqWheigtNeuron)

scaled = scaled * wheigt 
