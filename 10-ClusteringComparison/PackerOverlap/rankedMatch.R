library(readxl)
library(readr)
pvalues <- read_csv("resultMatrix.csv")
pvalues<-as.data.frame(pvalues)
rownames(pvalues)<-pvalues$X1
pvalues<-pvalues[,2:length(pvalues)]


CLUSTER_ID <- read_excel("/home/rama/Desktop/C elegans Transcriptome/RevisionNAR/GitHub/ScriptsModify2iter/Datasets/ClusteringfinalAssignment.xlsx")
CLUSTER_ID<-as.data.frame(CLUSTER_ID)
row.names(CLUSTER_ID)<-CLUSTER_ID$Cluster


result=list()
neurons=list()
resultsdata=data.frame(stringsAsFactors = FALSE)
for (cluster in colnames(pvalues)) {#
  x=pvalues[order(pvalues[,cluster]),cluster]
  names(x)<-rownames(pvalues)[order(pvalues[,cluster])]
  resultsdata=rbind(resultsdata,
    data.frame(Parick=cluster,
               Packer1=names(x[1]),
               pvalue1=as.double(pvalues[names(x[1]),cluster]),
               Packer2=names(x[2]),
               pvalue2=as.double(pvalues[names(x[2]),cluster]),
               Packer3=names(x[3]),
               pvalue3=as.double(pvalues[names(x[3]),cluster]),
               Packer4=names(x[4]),
               pvalue4=as.double(pvalues[names(x[4]),cluster]),
               Packer5=names(x[5]),
               pvalue5=as.double(pvalues[names(x[5]),cluster]),
               Packer6=names(x[6]),
               pvalue6=as.double(pvalues[names(x[6]),cluster])
               )
    )
  if (is.na(sum(resultsdata$pvalue1))) {
    break
  }
  
  print(paste(cluster,names(x[2])))
    x=list(x)
  names(x)<-cluster
  result=c(result,x)

}

rownames(resultsdata)<-resultsdata$Parick

resultsdata=resultsdata[CLUSTER_ID$Cluster,]
resultsdata=cbind(resultsdata,CLUSTER_ID[rownames(resultsdata),]$`Manual assignment to Neuron classes`)

write.csv(resultsdata,"rankedMatch.csv")



