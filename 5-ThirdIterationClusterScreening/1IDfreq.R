load("IterationData.rda")
IterationData=unique(IterationData)



oneIdClusters=NA
duplicatedIds=NA
IterationData<-cbind(IterationData,oneIdClusters)
IterationData<-cbind(IterationData,duplicatedIds)

reslist=c(NA*1:30)
names(reslist)<-1:30/10

for (parentcluster in unique(IterationData$parent_cluster)) { # fill columns oneIdClusters and duplicatedIds, coming from the subclusters in the same set



  for (dims in unique(IterationData[IterationData$parent_cluster==parentcluster,"dim"])) {

    nclusters=reslist
    oneId=reslist
    dupID=reslist


    for (SeuratResol in unique(IterationData[IterationData$parent_cluster==parentcluster,"resol"])){
      subset2test=IterationData$dim==dims & IterationData$resol==SeuratResol & IterationData$parent_cluster==parentcluster
      temp=IterationData[subset2test,]
      oneIdClusters=length(temp[temp$n_ident==1,]$n_ident)
      duplicatedIds=sum(duplicated(temp[temp$n_ident==1,]$cluster_id))

      
      IterationData[subset2test,"oneIdClusters"]<-oneIdClusters
      IterationData[subset2test,"duplicatedIds"]<-duplicatedIds

      nclusters[paste(SeuratResol)]=IterationData[subset2test,"n_clusters"][1]
      oneId[paste(SeuratResol)]=oneIdClusters
      dupID[paste(SeuratResol)]=duplicatedIds



    }

  }

}


ResultData=IterationData[IterationData$n_ident==1,]#filtering to get only 1 id subclusters
ResultData=ResultData[ResultData$duplicatedIds==0,]#filtering to get only subclusters coming from clustering conditions without duplicated ID assignements

TestData=data.frame()
for (parentcluster in unique(IterationData$parent_cluster)) {#Test data is a summary of 1IDs subcluster
  temp=ResultData[ResultData$parent_cluster==parentcluster,]
  for (dim in unique(temp$dim)) {
    for (resol in unique(temp[temp$dim==dim,]$resol)) {
      
      UniqueIDs=as.character(temp[temp$dim==dim & temp$resol==resol,]$cluster_id)
      n_UIs=length(UniqueIDs)
      UniqueIDs=paste(UniqueIDs[order(UniqueIDs)],collapse = ",")
      n_clust=temp[temp$dim==dim & temp$resol==resol,]$n_clusters[1]
      
      TestData<-rbind(TestData,
                      data.frame(parentcluster=parentcluster,
                                 dim=dim,
                                 resol=resol,
                                 UniqueIDs=UniqueIDs,
                                 stringsAsFactors = FALSE,
                                 n_clust=n_clust,
                                 n_UIs=n_UIs
                      )
      )
    }
  }
  
}




TestData=unique(TestData)

unlink("Ident_Clusters.txt")
fileConn<-file("Ident_Clusters.txt","a")
for (parentcluster in unique(IterationData$parent_cluster)) { #this write a file with the occurrence of 1ID clusters
  print(
    
    paste("parent cluster:", parentcluster, "n trials:",
          length(unique(IterationData[IterationData$parent_cluster==parentcluster,c("dim","resol")])[,1]))
    
  )
  writeLines(
    paste("parent cluster:", 
          parentcluster, 
          "n trials:",
          length(unique(IterationData[IterationData$parent_cluster==parentcluster,c("dim","resol")])[,1])),
    fileConn
  )
  freq=table(as.character(ResultData[ResultData$parent_cluster==parentcluster,]$cluster_id))
  freq=freq[order(freq,decreasing = TRUE)]
  writeLines("Individual Ocurrence",fileConn)
  writeLines(paste(as.character(names(freq)),collapse = "\t"),fileConn)
  writeLines(paste(as.character(freq),collapse = "\t"),fileConn)         
  print(freq)
  freq=table(TestData[TestData$parentcluster==parentcluster,]$UniqueIDs)
  freq=freq[order(freq,decreasing = TRUE)]
  writeLines("Combined Ocurrence",fileConn)
  writeLines(paste(as.character(names(freq)),collapse = "\t"),fileConn)
  writeLines(paste(as.character(freq),collapse = "\t"),fileConn)   
  writeLines(" ",fileConn)  
  print(freq[order(freq,decreasing = TRUE)])
}
close(fileConn)


