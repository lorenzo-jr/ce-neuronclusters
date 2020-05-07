source("InitializationScript.R")

####################################################################################

#Best N of PCs dimensions from first screening
dims=92 

#Best resolution parameter from first screening
SeuratResol=4

#Create PDF with all plots
pdf(paste("BoxPlotClusterMedianZscore_nPC_",paste(dims,"_Resol_",SeuratResol),".pdf"),paper="A4",height = 0.01, width = 0.05)
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
  
  # A sequential t-test from the highest to lowest ranked neuron
  # classes  applied iteratively to the successive pair
  pvalues=rep(1,length(NeuronOrder))
  names(pvalues)<-rev(NeuronOrder)
  for (i in 1:(length(NeuronOrder)-1)) {
    # Try two sample t-test. 
    test=try(t.test(data2test[i],data2test[i+1],alternative = "greater"),silent = TRUE)#manage t.test error
    if (!is.list(test)) {
      break
    }
    pvalues[i]=test$p.value
    if (test$p.value<0.05) {
      break
    }
  }
  
  # one sample t-test with mu=0 (in the scaled data zero is the mean)
  ptreshold=apply(data2test, 2, FUN = tryTtest,alternative = "greater")
  ptreshold <- sapply(ptreshold, '[[', 'p.value')
  
  # FDR adjustment of the p-value
  FDRtreshold=p.adjust(ptreshold,method = "fdr")
  
  #plot boxplot
  boxplot(data2test,
          ylim = c(-2, 22),
          outline = FALSE,
          notch = TRUE,
          ylab="Mean z-score",col = c(c("white","red")[((pvalues<1)*(FDRtreshold<0.05))+1])) #Apply treshold conditions to boxplot colors
 
  
  legend("topleft", paste("cluster ",cluster,". nº cells: ",sum(CeNeurons$seurat_clusters==cluster),". nPC: ",dims,". Resolution: ",SeuratResol), bty="n",cex = 2)
  
  
  
  
}
dev.off()

#save first iteration clusters
CeNeurons_1st_seurat_clusters=CeNeurons$seurat_clusters
save(file="CeNeurons_1st_seurat_clusters.rda",CeNeurons_1st_seurat_clusters)


#Feature plot 'ric-4', 'unc-31','unc-104'
pdf(file = "ric-4_unc-31_unc-104.pdf",paper = "A4")
FeaturePlot(object = CeNeurons, features = c('ric-4', 'unc-31','unc-104'))
dev.off()

listacells = list()
listacells[["unc-25"]] = WhichCells(CeNeurons,expression = `unc-25` > 0)
listacells[["unc-46"]] = WhichCells(CeNeurons,expression = `unc-46` > 0)
listacells[["eat-4"]] = WhichCells(CeNeurons,expression = `eat-4` > 0)
listacells[["cha-1"]] = WhichCells(CeNeurons,expression = `cha-1` > 0)
listacells[["unc-17"]] = WhichCells(CeNeurons,expression = `unc-17` > 0)

combi = combn(c("unc-25", "unc-46", "eat-4","cha-1", "unc-17"), 2, FUN = NULL, simplify = F)
NeuronalClustering <- CeNeurons

#check coexpression of markers
for (genepair in combi) {
  Idents(NeuronalClustering) <- NA
  Idents(NeuronalClustering, cells = listacells[[genepair[1]]]) <- genepair[1]
  Idents(NeuronalClustering, cells = listacells[[genepair[2]]]) <- genepair[2]
  Idents(NeuronalClustering, cells = intersect(listacells[[genepair[1]]],listacells[[genepair[2]]]) ) <-"Coexpress" 
  x <- TSNEPlot(NeuronalClustering, pt.size =0.1)
  ggsave(filename = paste("",paste(genepair,collapse = "_"),".pdf",sep = ""), x, width=7, height=4)


}

