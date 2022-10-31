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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/using_umap_clustering/re_2e_5/gene/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
Cluster <-read.table("../../../11_whole_genome_umap_leiden_pca5_k50_resolution2e-05.txt",header = T,sep = "\t") %>% as.data.frame()
all_gene$hotspot <-paste0(all_gene$V1,":",all_gene$V2,"-",all_gene$V3)
all_gene$eqtl <- paste0(all_gene$V4,":",all_gene$V5,"-",all_gene$V6)
all_gene$length <- all_gene$V3 - all_gene$V2
all_gene<-left_join(all_gene,Cluster[,c("cluster","hotspot")],by="hotspot")
colnames(all_gene)[7]<-"ENSEMBL"
all_gene<-all_gene[,c("cluster","hotspot","eqtl","length","ENSEMBL")]
all_gene_count <-all_gene%>%group_by(hotspot,cluster,length)%>%summarise(eqtl_gene_count=n())%>%data.frame()
all_gene_count$adjust_eqtl_gene_count <-all_gene_count$eqtl_gene_count/all_gene_count$length*1000


p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))

p1 <- ggplot(all_gene_count, aes(x=as.factor(cluster), y=log(adjust_eqtl_gene_count),group=cluster,fill = as.factor(cluster))) + 
    geom_boxplot()+
    labs(x="Cluster",y="Log(number of eQTL-eGene pairs per kb)")+
    scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    theme_bw()+
    ggtitle("eQTL-eGene pairs")+
    theme(legend.position="none")+
    p_theme
pdf("13_eqtl_gene_count_boxplot.pdf",height=4.5,width=4.5)
print(p1)
dev.off()


p <- ggplot(all_gene_count, aes(x=log(adjust_eqtl_gene_count),color=as.factor(cluster))) + 
geom_density(size=0.8)+
scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(number of eQTL-eGene pairs per kb)",color="Cluster") +
# coord_cartesian(xlim = c(0, 10))
ggtitle("eQTL-eGene pairs")+
p_theme
pdf("13_densityplot_eqtl_gene_count.pdf",height=4.5,width=4.9)
print(p)
dev.off()


library(clusterProfiler)
library(org.Hs.eg.db)

data = all_gene[,"ENSEMBL"]%>%unique()
data = as.vector(data)
annots <- select(org.Hs.eg.db, keys=data, 
                columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
annots <-annots[!duplicated(annots$ENSEMBL),]
all_gene <-left_join(all_gene,annots,by="ENSEMBL")

enrichment_and_plot <-function(i=NULL){
    dir_name =paste0("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/gene/",i)
    dir.create(dir_name)
    setwd(dir_name)
    cluster <-filter(all_gene,cluster==i)
    kegg <- enrichKEGG(cluster$ENTREZID, 
    organism = 'hsa', 
    keyType = 'kegg', 
    pvalueCutoff = 0.05,
    pAdjustMethod = 'fdr', 
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.05,
    use_internal_data = FALSE)
    kegg1 <-data.frame(kegg)
    write.table(kegg1,paste0(i,"_KEGG.txt"),col.names=T,row.names =F,quote=F,sep="\t")
    if(nrow(kegg)>0){
        pdf(paste0(i,"_KEGG.pdf"),width=5,height=5)
        p1 <-dotplot(kegg)
        print(p1)
        dev.off()
    }
    types <-c("BP","CC","MF")
    for(type in types){
        go <- enrichGO(cluster$ENTREZID, 
            OrgDb = 'org.Hs.eg.db',
            keyType = 'ENTREZID',
            ont=type,
            pvalueCutoff = 0.05,
            pAdjustMethod = 'fdr', 
            minGSSize = 10,
            maxGSSize = 500,
            qvalueCutoff = 0.05,   
            readable = FALSE,
            pool = FALSE)
        go1 <-data.frame(go)
        write.table(go1,paste0(i,"_GO_",type,".txt"),col.names=T,row.names =F,quote=F,sep="\t")
        if(nrow(go)>0){
            pdf(paste0(i,"_GO_",type,".pdf"))
            p <-dotplot(go)
            print(p)
            dev.off()
        }
        print(type)
    }
    print(i)
}

lapply(unique(all_gene$cluster),enrichment_and_plot)


test = bitr(x, #数据集
fromType="SYMBOL", #输入为SYMBOL格式
toType="ENSEMBL",  # 转为ENTERZID格式
OrgDb="org.Hs.eg.db")