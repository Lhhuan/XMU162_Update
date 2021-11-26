library(irlba)
library(data.table)
library(KKLClustering)
library(ggplot2)
data = fread("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/kmer/6/6mers.csv",data.table=FALSE)
# data = fread("/share/data0/GTEx/test_kmer/chr1_22_6kmers.csv", data.table=FALSE)
rownames(data) = data[,1]
data = data[,-1]

#-------------------------------------------------------------------------------
#                                                                              
#                                      PCA                                 
#                                                                              
#-------------------------------------------------------------------------------
Sys.time() 
pca = prcomp_irlba(data, n = 50)
pca_result = pca$x
rownames(pca_result) = rownames(data)
Sys.time() 
#-------------------------------------------------------------------------------
#                                                                              
#                                   clustering                                 
#                                                                              
#-------------------------------------------------------------------------------

result = kkl(pca_result,
             outlier_q = 0,
             down_n = 2000,
             knn_range = seq(5,200,5),
             iter = 50,
             compute_index =  "Calinski_Harabasz",
             assess_index = "Calinski_Harabasz",
             cores = 5,
             seed = 723)
Sys.time() 

table(result$cluster_df$cluster)
a <-result$cluster_df
a$name <-gsub(">","",a$name)
write.table(a,"01_kkl_result.txt",sep="\t",quote=F,col.names=T,row.names=F)

#-------------------------------------------------------------------------------
#                                                                              
#                                和seekr的结果对比                                 
#                                                                              
#-------------------------------------------------------------------------------

bed <- as.data.frame(read.table("/share/data0/GTEx/test_kmer/chr1_11_communities_5.bed.gz",
								header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
rn <- sapply(1:nrow(bed), function(i){
	x=bed[i,]
	paste0(">",x[1],":",x[2],"-",x[3])
})
all.equal(rn, rownames(data))
rownames(bed) = rn

aricode::ARI(bed$V4, result$cluster_df$cluster)
# 0.5479268
#------------------------tsne
pca = prcomp_irlba(data, n = 50)
pca_result = pca$x
tsne = Rtsne::Rtsne(pca_result, pca = FALSE,check_duplicates = FALSE)
df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = result$cluster_df$cluster)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = cluster), size = 0.8, alpha = 0.3)
ggsave("aa.png",p,dpi=72)