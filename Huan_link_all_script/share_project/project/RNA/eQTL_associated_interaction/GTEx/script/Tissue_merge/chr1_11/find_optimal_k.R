library(irlba)
library(data.table)
# library(KKLClustering)
library(NbClust)
# data = fread("/share/swap/huan/6mers.csv",data.table=FALSE)
# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_11/")
data = fread("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/6mers.csv", data.table=FALSE)

rownames(data) = data[,1]
data = data[,-1]

res<-NbClust(data, distance = "euclidean", min.nc=2, max.nc=10, 
            method = "kmeans", index = "all")

# NbClust(data = data,diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15, method = NULL, index = "all", alphaBeale = 0.1)

# res<-NbClust(data, distance = "euclidean", min.nc=2, max.nc=10, 
#             method = NULL, index = "all")

res<-NbClust(data, distance = "euclidean", min.nc=2, max.nc=10, 
            method = "kmeans", index = "all")

            
res<-NbClust(data, distance = "euclidean", min.nc=2, max.nc=10, 
             index = "all")