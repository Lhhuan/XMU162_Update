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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_2cluster/gene/")
all_gene <-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge_egene_eQTL_5e_8.bed.gz",header = F,sep = "\t") %>% as.data.frame()
Cluster <-read.table("../../11_whole_genome_leiden_pca3_k50_resolution2e-05.txt",header = T,sep = "\t") %>% as.data.frame()
all_gene$hotspot <-paste0(all_gene$V1,":",all_gene$V2,"-",all_gene$V3)
all_gene$eqtl <- paste0(all_gene$V4,":",all_gene$V5,"-",all_gene$V6)
all_gene$length <- all_gene$V3 - all_gene$V2
all_gene<-left_join(all_gene,Cluster[,c("cluster","hotspot")],by="hotspot")
colnames(all_gene)[7]<-"ENSEMBL"
all_gene<-all_gene[,c("cluster","hotspot","eqtl","length","ENSEMBL")]
all_gene_count <-all_gene%>%group_by(hotspot,cluster,length)%>%summarise(eqtl_gene_count=n())%>%data.frame()
all_gene_count$adjust_eqtl_gene_count <-all_gene_count$eqtl_gene_count/all_gene_count$length*1000
all_gene_count1 <-filter(all_gene_count,cluster %in%c(1:2))

p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5))



all_gene_count1$cluster <-factor(all_gene_count1$cluster,levels=c(1,2))
p1 <- ggplot(all_gene_count1, aes(x=as.factor(cluster), y=log(adjust_eqtl_gene_count),group=cluster,fill = as.factor(cluster))) + 
    geom_boxplot()+
    labs(x="Cluster",y="Log(number of pairs per kb)")+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B","#E277C1","#7F7F7F","#BCBC22","#15BECE"))+
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    # scale_fill_manual(values=c("#8C554B","#C22324","#1E77B4","#FF7F0E","#9567BD","#2CA02C"))+
    # scale_fill_manual(values=c("#C22324","#FF7F0E","#1E77B4","#9567BD","#8C554B","#2CA02C"))+
    scale_fill_manual(values=c("#1E77B4","#FF7F0E"))+
    theme_bw()+
    ggtitle("eQTL-eGene pairs")+
    theme(legend.position="none")+
    stat_compare_means(method = "wilcox.test",comparisons =list(c(1,2)),label="p.signif",method.args = list(alternative = "less"))+
    p_theme
pdf("13_eqtl_gene_count_boxplot.pdf",height=3.1,width=3)
print(p1)
dev.off()


p <- ggplot(all_gene_count1, aes(x=log(adjust_eqtl_gene_count),color=as.factor(cluster))) + 
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

#==============
enrichment_and_plot <-function(i=NULL){
    dir_name =paste0("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/gene/",i)
    dir.create(dir_name)
    setwd(dir_name)
    cluster <-filter(all_gene1,cluster==i)
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

library(venn)
# x <- list(First = 1:20, Second = 10:30, Third = sample(25:50, 15))
# venn(x, ilab = TRUE, zcolor = "style")

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/whole_genome/3pca_1e_4/gene/")
types<-c("GO_BP","GO_CC","GO_MF","KEGG")
for(type in  types){
    aaa<-list()
    for(i in 1:6){
        en <-read.table(paste0("./",i,"/",i,"_",type,".txt"),header = T,sep = "\t") %>% as.data.frame()
        if(nrow(en)>0){
        aaa[[length(aaa)+1]] <-en$Description
        names(aaa)[length(aaa)] <-paste0("C",i)
        print(i)
        }
    }
    pdf(paste0("13_",type,"_venn.pdf"),width=3.1,height=3.1)
    venn(aaa, ilab = TRUE, zcolor = "style")
    dev.off()
    for(j in c(1:length(aaa))){
        b <-aaa
        all <-aaa[[j]]
        CC <- names(aaa[j])
        CC <-gsub("C","",CC)
        b[[j]]=NULL
        names(b)=NULL
        other <-unlist(b)
        specific <- setdiff(all,other)%>%data.frame()
        write.table(specific,paste0("./",CC,"/",CC,"_",type,"_specific.txt"),col.names=F,row.names =F,quote=F,sep="\t")
    }
}

#===================================================================compare eQTL and hotspots

eqtl_count <-all_gene%>%group_by(eqtl)%>%summarise(eqtl_gene_count=n())%>%data.frame()
hotspot_eqtl_count <-all_gene[,c("hotspot","ENSEMBL")]%>%unique()%>%group_by(hotspot)%>%summarise(gene_count=n())%>%data.frame()
egene_cover_hotspot <-all_gene[,c("hotspot","ENSEMBL")]%>%unique()%>%group_by(ENSEMBL)%>%summarise(number_of_hotspot=n())%>%data.frame()

# x <- list(First = 1:20, Second = 10:30, Third = sample(25:50, 15))


hotspot1 <-filter(all_gene,hotspot=="chr1:865450-866545")

p <- ggplot(eqtl_count, aes(x=eqtl_gene_count)) + 
geom_density(size=0.8)+
# scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(Number of egene in eQTL)") +
# coord_cartesian(xlim = c(0, 10))
ggtitle("Number of egene in eQTL")+
p_theme
pdf("13_densityplot_gene_count_in_eQTL.pdf",height=4.5,width=4.9)
print(p)
dev.off()

#==========================
p <- ggplot(hotspot_eqtl_count, aes(x=gene_count)) + 
geom_density(size=0.8)+
# scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(Number of egene in hotspot)") +
# coord_cartesian(xlim = c(0, 10))
ggtitle("Number of egene in hotspot")+
p_theme
pdf("13_densityplot_gene_count_in_hotspot.pdf",height=4.5,width=4.9)
print(p)
dev.off()

h_eqtl_egene_number<-all_gene%>%group_by(hotspot,ENSEMBL)%>%summarise(number_of_eQTL=n())%>%data.frame()
eqtl_egene_number<-all_gene%>%group_by(ENSEMBL)%>%summarise(number_of_total_eQTL=n())%>%data.frame()

h_eqtl_egene_number <- left_join(h_eqtl_egene_number,eqtl_egene_number,by="ENSEMBL")

h_eqtl_egene_number$ratio <-h_eqtl_egene_number$number_of_eQTL/h_eqtl_egene_number$number_of_total_eQTL
ratios <-bind_rows()

egene_h <- h_eqtl_egene_number%>%group_by(ENSEMBL)%>%summarise(count=n())%>%data.frame()
egene_qtl <- all_gene[,c("eqtl","ENSEMBL")]%>%group_by(ENSEMBL)%>%summarise(count=n())%>%data.frame()

egene_h$Class<-"Hotspot"
egene_qtl$Class <-"eQTL"
egene <-bind_rows(egene_h,egene_qtl)

p <- ggplot(egene, aes(x=log(count),color=Class)) + 
geom_density(size=0.8)+
scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(Number of class in egene)") +
# coord_cartesian(xlim = c(0, 10))
ggtitle("Number of class in egene")+
p_theme
pdf("13_densityplot_gene_distribution_in_hotspot_eqtl.pdf",height=4.5,width=4.9)
print(p)
dev.off()

p <- ggplot(h_eqtl_egene_number, aes(x=ratio)) + 
geom_density(size=0.8)+
# scale_color_d3("category20") + 
theme_bw()+
labs(x="Log(Ratio of hotspot cover gene)") +
# coord_cartesian(xlim = c(0, 10))
ggtitle("Ratio of hotspot cover gene")+
p_theme
pdf("13_densityplot_hotsopt_cover_gene_ratio.pdf",height=4.5,width=4.9)
print(p)
dev.off()


test = bitr(x, #数据集
fromType="SYMBOL", #输入为SYMBOL格式
toType="ENSEMBL",  # 转为ENTERZID格式
OrgDb="org.Hs.eg.db")