Script files:

- InitializationScript.R: used to load library and functions and to initialize Seurat object

- SecondIterationClusterScreening.R: this script iterates by changing the numbers of principal components  (between 3 to 92)  and resolution parameters (from 0.1 to 3) used for Louvain clustering.

- 1IDfreq.R: This script counts neuronal identity occurrence during clustering screening. Only clusterings condition giving not repeated identities are tacked into account.

- IterationDataAnalysis.py: This script counts neuronal identities occurring during clustering screening and saves conditions that generate those identities. Also includes conditions giving repeated ids.

Output files:

- IterationData.rda and IterationData.csv: Dataframe and table with the results of the second clustering iteration (SecondIterationClusterScreening.R outputs)

- Ident_Clusters.txt: Neuronal identity occurrence count in subclusters generated during screening of clustering conditions (1IDfreq.R output). 

- ClusterAnalysis.csv: List of neuronal identities obtained by different clustering conditions (IterationDataAnalysis.py output)
