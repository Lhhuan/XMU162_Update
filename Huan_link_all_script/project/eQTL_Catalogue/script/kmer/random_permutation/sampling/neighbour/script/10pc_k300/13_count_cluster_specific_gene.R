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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
Cluster <-read.table("../11_chr1_louvain_pca10_k300.txt",header = T,sep = "\t") %>% as.data.frame()
chr1_gene <- filter(all_gene,V1=="chr1")
chr1_gene$hotspot <-paste0(chr1_gene$V1,":",chr1_gene$V2,"-",chr1_gene$V3)
chr1_gene<-left_join(chr1_gene,Cluster[,c("cluster","hotspot")],by="hotspot")
colnames(chr1_gene)[4]<-"ENSEMBL"
chr1_gene<-chr1_gene[,c("cluster","hotspot","ENSEMBL")]
chr1_gene_count <-chr1_gene%>%group_by(hotspot,cluster)%>%summarise(gene_count=n())%>%data.frame()


p1 <- ggplot(chr1_gene_count, aes(x=cluster, y=gene_count,group=cluster)) + 
    geom_boxplot()+
    theme_bw()
pdf("13_gene_count_boxplot.pdf")
print(p1)
dev.off()


library(clusterProfiler)
library(org.Hs.eg.db)

data = chr1_gene[,"ENSEMBL"]%>%unique()
data = as.vector(data)
annots <- select(org.Hs.eg.db, keys=data, 
                columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
annots <-annots[!duplicated(annots$ENSEMBL),]
chr1_gene <-left_join(chr1_gene,annots,by="ENSEMBL")

enrichment_and_plot <-function(i=NULL){
    dir_name =paste0("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/enrichment/",i)
    dir.create(dir_name)
    setwd(dir_name)
    cluster <-filter(chr1_gene,cluster==i)
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

lapply(unique(chr1_gene$cluster),enrichment_and_plot)


test = bitr(x, #数据集
fromType="SYMBOL", #输入为SYMBOL格式
toType="ENSEMBL",  # 转为ENTERZID格式
OrgDb="org.Hs.eg.db")