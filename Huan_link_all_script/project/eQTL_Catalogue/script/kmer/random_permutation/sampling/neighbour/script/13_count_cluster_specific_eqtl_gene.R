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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/gene/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
whole_genome_gene <-read.table("../11_whole_genome_leiden_pca5_k30_resolution1e-04.txt",header = T,sep = "\t") %>% as.data.frame()
whole_genome_gene$hotspot <-paste0(whole_genome_gene$V1,":",whole_genome_gene$V2,"-",whole_genome_gene$V3)
whole_genome_gene$eqtl <- paste0(whole_genome_gene$V4,":",whole_genome_gene$V5,"-",whole_genome_gene$V6)
whole_genome_gene<-left_join(whole_genome_gene,Cluster[,c("cluster","hotspot")],by="hotspot")
colnames(whole_genome_gene)[7]<-"ENSEMBL"
whole_genome_gene<-whole_genome_gene[,c("cluster","hotspot","eqtl","ENSEMBL")]
whole_genome_gene_count <-whole_genome_gene%>%group_by(hotspot,cluster)%>%summarise(eqtl_gene_count=n())%>%data.frame()


p1 <- ggplot(whole_genome_gene_count, aes(x=as.factor(cluster), y=eqtl_gene_count,group=cluster,fill = as.factor(cluster))) + 
    geom_boxplot()+
    labs(x="Cluster",y="Number of eQTL-eGene pairs")+
    theme_bw()+
    theme(legend.position="none")
pdf("13_eqtl_gene_count_boxplot.pdf",height=6,width=6)
print(p1)
dev.off()


p <- ggplot(whole_genome_gene_count, aes(x=eqtl_gene_count,color=as.factor(cluster))) + 
geom_density()+
scale_color_d3("category20") + 
theme_bw()+
xlab("Number of eQTL-eGene pairs") #+
# coord_cartesian(xlim = c(0, 50))

pdf("13_densityplot_eqtl_gene_count.pdf",height=6,width=7)
print(p)
dev.off()


library(clusterProfiler)
library(org.Hs.eg.db)

data = whole_genome_gene[,"ENSEMBL"]%>%unique()
data = as.vector(data)
annots <- select(org.Hs.eg.db, keys=data, 
                columns=c("SYMBOL","ENTREZID"), keytype="ENSEMBL")
annots <-annots[!duplicated(annots$ENSEMBL),]
whole_genome_gene <-left_join(whole_genome_gene,annots,by="ENSEMBL")

enrichment_and_plot <-function(i=NULL){
    dir_name =paste0("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/09_chr/chr1/gene/",i)
    dir.create(dir_name)
    setwd(dir_name)
    cluster <-filter(whole_genome_gene,cluster==i)
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

lapply(unique(whole_genome_gene$cluster),enrichment_and_plot)


test = bitr(x, #数据集
fromType="SYMBOL", #输入为SYMBOL格式
toType="ENSEMBL",  # 转为ENTERZID格式
OrgDb="org.Hs.eg.db")