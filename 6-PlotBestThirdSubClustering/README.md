Script files:

- InitializationScript.R: used to load library and functions and to initialize Seurat object

- PlotBestSecondSubClustering.R: this script is used for piloting in UMAP space using best clustering parameters chose  in previous step. It uses manually selected parameters found on ParameterChoose.csv

Output files from PlotBestSecondSubClustering.R:

- BoxPlotClusterMeanZscore_Cluster_(8.0, 10.0 or 10.1)_nPC_x_Resol_x.pdf : UMAP reduced space for second iteration of clustering. Also includes for each subcluster the corresponding boxplot with automatic assignment of neuronal identities (in red)

- CustomIterationClusters3dit: Cell clusters results from the third cluster iteration
