Script files:

- InitializationScript.R: used to load library and functions and to initialize seurat object

- PlotBestFirstClustering.R: this script is used for plotting in tSNE space using best clustering parameters found  in previous step (folder 1-FirstClusteringScreening; PCs = 92; Resolution = 4). Also plots in tSNE spaces highlighting some neuronal markers are produced at this step.


Output files from PlotBestFirstClustering.R:

- BoxPlotClusterMedianZscore_nPC_92_Resol_4.pdf: tSNE reduced space for the clustering using PCs: 92 and resolution: 4. Also includes for each cluster the corresponding boxplot with automatic assignment of neuronal identities (in red)

- CeNeurons_1st_seurat_clusters.rda: Cell clusters from the first cluster iteration

- ric-4_unc-31_unc-104.pdf: Cell expressing the neuronal markers ric-4, unc-31 and unc-104 plotted in tSNE reduced space

- Other PDF files: Combinatorial expression for neuronal markers unc-25, unc-46, eat-4, cha-1 and unc-17 in tSNE reduced space. 
