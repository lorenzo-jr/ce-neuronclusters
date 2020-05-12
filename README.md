# Seurat clustering and ploting scrpits used in the article "Combining single-cell RNA-sequencing with a molecular atlas unveils new markers for Caenorhabditis elegans neuron classes".

Authors: Ramiro Lorenzo, Michiho Onizuka,Matthieu Defrance, Patrick Laurent

Datasets used for the analysis are in Datasets folder.

The scripts are organized as the workflow clustering analysis.

1-FirstClusterScreening:
Clustering screening using Seurat by changing PCA dimensions (between 40 and 100) and resolution (from 1 to 8). Louvain algorithm was used for clustering.

2-PlotBestFirstClustering:
Plot tSNE space and boxplots corresponding to each cluster using best parameters obtained from 1-FirstClusterScreening. Some neuronal markers are also plotted at this step on tSNE space.

3-SecondIterationClusterScreening:
Second round of clustering on the 64 parent clusters obtained in 2-PlotBestFirstClustering. Single-cell profiles coming from parent clusters were independently re-clustered. Multiple re-clustering trials were generated for resolution values from 0.1 to 3 and for PCs values from 3 to 92.

4-PlotBestSecondSubClustering:
Single UMAP space plots for each individual parent clusters and refinement of new clusters by manual curation.

5-ThirdIterationClusterScreening:
Third round of clustering on parent clusters 8.0, 10.0 and 10.1 obtained in 4-PlotBestSecondSubClustering. Single-cell profiles
coming from parent clusters were independently re-clustered. Multiple re-clustering trials were generated for resolution values from 0.1 to 3 and for PCs values from 3 to 92.

6-PlotBestThirdSubClustering:
Single UMAP space plot for each individual clusters and refinement of new clusters by manual curation.

7-FinalPlotAndClusteringMarkers:
Final results plots, expression and markers analysis.

8-DoubletAnalysis:
Analysis of cell doublets using DoubletFinder

9-Psudotime:
Pseudotime analysis using 5 sets of clusters 0_2, 3, 4, 16 and 42.1-56.

10-ClusteringComparison:
Cluster wise comparison against Cao where only cell names were compared and Packer where marker genes were compared between clusters. Other statistics are also found in this folder.

Datasets:
Folder containing datasets used for the iterative clustering procedures
