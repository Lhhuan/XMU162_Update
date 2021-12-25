library(irlba)
library(data.table)
library(KKLClustering)
library(ggplot2)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")

org<-read.table("../../../output/Tissue_merge/Cis_eQTL/02_hotspot_mark_annotatiob_fraction.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org$hotspot <-paste(org$chr,org$start,org$end,sep="_")
org <-org[,c("Marker","overlap_fraction","hotspot")]
org1 <-org%>%group_by(Marker,hotspot)%>%summarise(sum_overlap_fraction =sum(overlap_fraction))%>%data.frame()
# org1%>%group_by(Marker)%>%summarise(count=n())
org2 <-dcast(org1[,c(2,1,3)],hotspot~Marker)
rownames(org2)<-org2$hotspot
data <-org2[,-1]

#-------------------------------------------------------------------------------
#                                                                              
#                                      PCA                                 
#                                                                              
#-------------------------------------------------------------------------------
Sys.time() 
pca = prcomp_irlba(data, n = 9)
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

#------------------------tsne
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1), axis.title.y = element_text(size = 10), #,size=1.2
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_blank(), plot.title = element_text(hjust = 0.5)) 



tsne = Rtsne::Rtsne(pca_result, pca = FALSE,check_duplicates = FALSE)
df = data.frame(tsne1 = tsne$Y[,1], tsne2 = tsne$Y[,2], cluster = result$cluster_df$cluster)
p =ggplot(df) + geom_point(aes(x = tsne1, y = tsne2, colour = factor(cluster)), size = 0.4, alpha = 0.3)+ggtitle("KKL")+theme_bw() + p_theme
ggsave("05_KKL_tsne_cluster.png",p,dpi=150,width=5.5,height=5)

#--------------------umap
rumap <- uwot::umap(pca_result, n_neighbors = 30)

df <-data.frame(UMAP_1=rumap[,1],UMAP_2=rumap[,2],cluster = result$cluster_df$cluster)
p =ggplot(df) + geom_point(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster)), size = 0.3, alpha = 0.3)+ggtitle("KKL")+theme_bw() + p_theme
ggsave("05_KKL_Umap_cluster.png",p,dpi=150,width=5.5,height=5)