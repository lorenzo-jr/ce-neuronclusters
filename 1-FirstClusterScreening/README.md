Script files:

- InitializationScript.R: used to load library and functions and to initialize seurat object

- FirstClusterScreening.R: this script iterates by changing the numbers of principal components  (between 40 and 100)  and resolution parameters (from 1 to 8) used for Louvain clustering. Finally generates tree tables with summary results for each condition (output files bellow).

Output files from FirstClusteringScreening.R:

- duplicated.csv: Table indicating number of clusters sharing same identities 

- n1idclusters.csv: Table indicating total number of clusters with single identity

- nclusters.csv: Total number of clusters

Summary in ClusteringOptimization.xlsx. Maximum number of clusters with one non duplicated neuronal identity was 33. The lowest of number PCs yielding 33 identities was 92 (resolution 4).
