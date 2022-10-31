library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)
library(reshape2)
library(R.utils)
library(ggsci)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/gene/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
Cluster <-read.table("../../11_whole_genome_leiden_pca3_k50_resolution1e-04.txt",header = T,sep = "\t") %>% as.data.frame()
tmp1 <-strsplit(Cluster$hotspot,":")
tmp2 <- do.call(rbind,tmp1)%>%data.frame()
Cluster$CHR <- tmp2$X1
cluster_chr <-Cluster %>%group_by(cluster,CHR)%>%summarise(count=n())%>%data.frame()
cluster_n <- Cluster%>%group_by(cluster)%>%summarise(cluster_count=n())%>%data.frame()
cluster_chr <-left_join(cluster_chr,cluster_n,by="cluster")
cluster_chr$ratio <- cluster_chr$count/cluster_chr$cluster_count


all_gene$hotspot <-paste0(all_gene$V1,":",all_gene$V2,"-",all_gene$V3)
all_gene$eqtl <- paste0(all_gene$V4,":",all_gene$V5,"-",all_gene$V6)
all_gene$length <- all_gene$V3 - all_gene$V2
all_gene<-left_join(all_gene,Cluster[,c("cluster","hotspot")],by="hotspot")
colnames(all_gene)[7]<-"ENSEMBL"
all_gene<-all_gene[,c("cluster","hotspot","eqtl","length","ENSEMBL")]
all_gene_count <-all_gene%>%group_by(hotspot,cluster,length)%>%summarise(eqtl_gene_count=n())%>%data.frame()
all_gene_count$adjust_eqtl_gene_count <-all_gene_count$eqtl_gene_count/all_gene_count$length*1000
all_gene_count1 <-filter(all_gene_count,cluster %in%c(1:6))

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

library(clusterProfiler)
library(org.Hs.eg.db)

data = all_gene[,"ENSEMBL"]%>%unique()
data = as.vector(data)
annots <- select(org.Hs.eg.db, keys=data, 
                columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
annots <-annots[!duplicated(annots$ENSEMBL),]
all_gene <-left_join(all_gene,annots,by="ENSEMBL")
all_gene1 <-unique(all_gene[,c("cluster","ENSEMBL","SYMBOL","ENTREZID")])

#======
bbb<-list()
for(i in 1:6){
    itmp <-filter(all_gene1,cluster==i)
    bbb[[i]]<-unique(itmp$ENSEMBL)
    names(bbb)[i] <-paste0("C",i)
    print(i)
    # print(nrow(filter(all_gene1,cluster==i)%>%dplyr::select(ENSEMBL)%>%unique()))
}
library(venn)
pdf("13_egene_venn.pdf",width=3.1,height=3.1)
venn(bbb, ilab = TRUE, zcolor = "style")
dev.off()

rs <-data.frame()
for(j in c(1:length(bbb))){
    b <-bbb
    all <-bbb[[j]]
    CC <- names(bbb[j])
    CC <-gsub("C","",CC)
    b[[j]]=NULL
    names(b)=NULL
    other <-unlist(b)
    specific <- setdiff(all,other)%>%data.frame()
    colnames(specific) <-"ENSG"
    specific$cluster <-j 
    rs <-bind_rows(rs,specific)
}

gene_cluster_n <-all_gene1%>%group_by(ENSEMBL)%>%summarise(count=n())%>%data.frame()
share_gene <-filter(gene_cluster_n,count==6)
share_gene$cluster <-"Share"
colnames(share_gene)[1] <-"ENSG"
share_gene <- share_gene[,c("ENSG","cluster")]
rs$cluster <-as.character(rs$cluster)
rs <- bind_rows(rs,share_gene)
other <-setdiff(unique(gene_cluster_n$ENSEMBL),rs$ENSG)
rs<-bind_rows(rs,data.frame(ENSG=other,cluster="Other"))
# random <- data.frame()
# for(i in 1:6){
#     c_specific <-filter(rs,cluster==i)
#     set.seed(123)
#     tmp3 <-c_specific[sample(1:nrow(c_specific), 20, replace = FALSE),]
#     random <- bind_rows(random,tmp3)
# }
# set.seed(123)
# tmp3 <-share_gene[sample(1:nrow(share_gene), 20, replace = FALSE),]
# random$cluster <- as.character(random$cluster)
# random <- bind_rows(random,tmp3[,c(1,3)])
#==============
gene_exp <-read.table("/share/data0/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",header = T,sep = "\t",skip=2) %>% as.data.frame()
gene_exp$Name<-gsub("\\..*","",gene_exp$Name)
colnames(gene_exp) <-gsub("\\.","_",colnames(gene_exp))
colnames(gene_exp) <-gsub("___","_",colnames(gene_exp))
colnames(gene_exp) <-gsub("__","_",colnames(gene_exp))
colnames(gene_exp)[1] <-"ENSG"
gene_exp$sd <- apply(gene_exp[,3:ncol(gene_exp)],1,sd)
gene_exp <-gene_exp[order(-gene_exp$sd),]
dat <-gene_exp[!duplicated(gene_exp$ENSG),]
fgene <-inner_join(rs,dat,rs,by="ENSG")
fgene <-fgene[order(-fgene$sd),]
fgene$Name <-paste(fgene$cluster,fgene$ENSG,sep="-")
rownames(fgene)<-fgene$Name

#==================
# fgene <- filter(fgene,cluster!="Other")
library(pheatmap)
color2 = colorRampPalette(c('#c6dbef','#6baed6','#2171b5'))(50)
pdf("13_1_gene_expression_tissue_heatmap_top500.pdf",width=15,height=50)
pheatmap(as.matrix(fgene[1:500,4:57]),cluster_rows = TRUE, cluster_cols = FALSE,angle_col = 90,color= color2,display_numbers = FALSE,show_rownames = T,cellwidth= 10)
dev.off()
log_fgene <- log(fgene[,4:57]+0.0001)
pdf("13_1_gene_expression_tissue_heatmap_top500_log.pdf",width=15,height=55)
pheatmap(as.matrix(log_fgene[1:500,]),cluster_rows = TRUE, cluster_cols = FALSE,angle_col = 90,color= color2,display_numbers = FALSE,show_rownames = T,cellwidth= 10)
dev.off()
#====================
pca<-prcomp(log_fgene[1:5000,])#
pc <-pca$x
df <-data.frame(pc[,1:2],Name=rownames(pc))
p <- ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(size = 1, alpha = 0.6) + 
  theme_bw()
ggsave("./egene_cluster/log/13_1_top_gene_pca.png", p, width = 6, height = 5.3)

library(kernlab)
pca_list <-list()
pca_list[[1]]<- kpca(~.,data=log_fgene[1:2000,], kernel="rbfdot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)
pca_list[[2]]<- kpca(~.,data=log_fgene[1:2000,], kernel="laplacedot", kpar = list(sigma = 0.01))
pca_list[[3]]<- kpca(~.,data=log_fgene[1:2000,], kernel="tanhdot", kpar = list(scale = 0.01, offset = 1))
pca_list[[4]]<- kpca(~.,data=log_fgene[1:2000,], kernel="besseldot", kpar = list(sigma = 1, order = 1, degree = 1))
pca_list[[5]]<- kpca(~.,data=log_fgene[1:2000,], kernel="anovadot", kpar = list(degree = .1))
pca_list[[6]]<- kpca(~.,data=log_fgene[1:2000,], kernel="vanilladot", kpar = list())
pca_list[[7]]<- kpca(~.,data=log_fgene[1:2000,], kernel="polydot", kpar = list(degree = 2, scale = 1, offset = 5))
pca_list[[8]] <-kpca(~.,data=log_fgene[1:2000,], kernel='splinedot', kpar = list())

names(pca_list)=c("rbfdot","laplacedot","tanhdot","besseldot","anovadot","vanilladot","polydot",'splinedot')
# save(pca_list,file = "08_kpca.Rdata")
plot <-function(i=NULL){
    aa=as.data.frame(pca_list[[i]]@rotated)
    colnames(aa)[1:2] <- c("PC1","PC2")
    p =ggplot(aa) + geom_point(aes(x = PC1, y = PC2), size = 0.1, alpha = 0.3)+ggtitle(names(pca_list)[i])+theme_bw() + p_theme
    return(p)
}

plist <-lapply(c(1:8),plot)
pdf(paste0("./egene_cluster/log/13_1_kpca_sampling_result.pdf"),width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()
#===============================kpca
num<-round(nrow(fgene)*0.3)
library(kernlab)
pca1<- kpca(~.,data=log_fgene[1:num,], kernel="anovadot",kpar = list(sigma = 0.1),features = 0, th = 1e-4)
aa=as.data.frame(pca1@rotated)
colnames(aa)[1:2] <- c("PC1","PC2")

p <- ggplot(aa, aes(x = PC1, y = PC2)) + 
  geom_point(size = 1, alpha = 0.6) + 
  theme_bw()
ggsave(paste0("./egene_cluster/log/13_1_top_gene_kpca_anovadot_top_",num,".png"), p, width = 6, height = 5.3)
#======================================================scree plot
library(ggplot2)
library(igraph)
library(dplyr)
library(ggsci)
variances <- data.frame(pca1@eig)
variances$Dimensions <- 1:nrow(variances)
variances$y <- variances$pca1.eig/sum(variances$pca1.eig) *100
variances$sum_y <-cumsum(variances$pca1.eig)/sum(variances$pca1.eig)*100
head(which(variances$sum_y>60))

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), 
                                                #   axis.title.y = element_text(size = 10),
                                                #   axis.title.x = element_text(size = 10),
                                                  axis.line = element_line(colour = "black"),legend.position ="none",plot.title = element_text(hjust = 0.5))


p1 <-ggplot(data = variances[1:30,], mapping = aes(x =Dimensions, y = y)) + 
    geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.9)+
    p_theme+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(breaks = seq(1, 30, by = 1))+
    labs(y="Percentage of explained variances",title="Scree plot")
ggsave(paste0("./egene_cluster/log/13_1_gene_expression_top_",num,"_kPCA_screeplot_whole_genome.pdf"),p1,dpi=300,width=6,height=5) 
#===================================umap,tsne
pcn = 2
# reduction
dat = pca1@rotated[,1:pcn]
rownames(dat)=1:nrow(dat)

set.seed(123)
umap <- uwot::umap(dat)
save(umap,file=paste0("./egene_cluster/log/13_1_top_",num,"_gene_expression_",pcn,"kpca_umap.Rdata"))
set.seed(123)
tsne = Rtsne::Rtsne(dat, initial_dims = pcn, pca = FALSE)
save(tsne,file=paste0("./egene_cluster/log/13_1_top_",num,"_gene_expression_",pcn,"kpca_tsne.Rdata"))
#================================graph 
k = 30
set.seed(1)
dat_knn = RcppHNSW::hnsw_knn(X = dat, k = k)
dat_edge = reshape2::melt(dat_knn$idx)
dat_graph = igraph::graph_from_edgelist(as.matrix(dat_edge[,c(1,3)]), directed = FALSE)
c = igraph::cluster_louvain(dat_graph)

# # plot
# # df = data.frame(dat[,1:2], umap$layout, tsne$Y, as.factor(c$membership) )
df = data.frame(dat[,1:2], umap, tsne$Y, as.factor(c$membership) )
colnames(df) = c(paste0("PC_",1:2), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
df$hotspot <- rownames(as.data.frame(pca1@xmatrix))[as.numeric(rownames(dat))]

p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 0.8, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("./egene_cluster/log/13_1_top_",num,"_gene_expression_louvain_pca",pcn,"_k",k,"_tsne.png"), p, width = 6, height = 5.3)

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 0.8, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("./egene_cluster/log/13_1_top_",num,"_gene_expression_louvain_pca",pcn,"_k",k,"_umap.png"), p, width = 6, height = 5)

#================leiden
resolution=0.0022
set.seed(40)#k=30
g <-cluster_leiden(dat_graph, resolution_parameter=resolution)
#plot
df = data.frame(dat, umap, tsne$Y, as.factor(g$membership) )
colnames(df) = c(paste0("PC_",1:2), paste0("UMAP_",1:2), paste0("tSNE_",1:2), "cluster")
#----
p <- ggplot(df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("./egene_cluster/log/13_1_top_",num,"_gene_expression_leiden_pca",pcn,"_k",k,"_resolution",resolution,"_tsne.png"), p, width = 6, height = 5.3)
#-----
p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_d3("category20") + 
  theme_bw()
ggsave(paste0("./egene_cluster/log/13_1_top_",num,"_gene_expression_leiden_pca",pcn,"_k",k,"_resolution",resolution,"_umap.png"), p, width = 6, height = 5.3)

df$all_gene <- rownames(pca1@xmatrix)
tmp1 <-strsplit(df$all_gene,"-")
tmp2 <- do.call(rbind,tmp1)%>%data.frame()
colnames(tmp2) <- c("hotspot_cluster","ENSG")
df <- bind_cols(df,tmp2)
df <-df%>%dplyr::select(-all_gene)
no_other <- filter(df,hotspot_cluster!="Other")
aricode::ARI(no_other$cluster,no_other$hotspot_cluster)
for (i in 1:6){
    cg <- filter(df,hotspot_cluster==i)
    print(aricode::ARI(cg$cluster,cg$hotspot_cluster))
}