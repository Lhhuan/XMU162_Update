library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(pheatmap)
library(reshape2)

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/5/kmer/extend/anno")
All <-read.table("../communities.bed.gz",header = F,sep = "\t") %>% as.data.frame()
CHROMATIN_Accessibility <-read.table("CHROMATIN_Accessibility.bed.gz",header = F,sep = "\t") %>% as.data.frame()
CTCF<-read.table("CTCF.bed.gz",header = F,sep = "\t") %>% as.data.frame()
TFBS<-read.table("TFBS.bed.gz",header = F,sep = "\t") %>% as.data.frame()

H3K27ac<-read.table("H3K27ac.bed.gz",header = F,sep = "\t") %>% as.data.frame() #active
H3K27me3<-read.table("H3K27me3.bed.gz",header = F,sep = "\t") %>% as.data.frame() #repressed
H3K36me3<-read.table("H3K36me3.bed.gz",header = F,sep = "\t") %>% as.data.frame() #transcribed
H3K4me1<-read.table("H3K4me1.bed.gz",header = F,sep = "\t") %>% as.data.frame() #enhancer
H3K4me3<-read.table("H3K4me3.bed.gz",header = F,sep = "\t") %>% as.data.frame() #promoter
H3K9ac<-read.table("H3K9ac.bed.gz",header = F,sep = "\t") %>% as.data.frame() #active
H3K9me3<-read.table("H3K9me3.bed.gz",header = F,sep = "\t") %>% as.data.frame() #gene silencing

# Active_histone <-bind_rows(H3K27ac,H3K36me3,H3K4me1,H3K4me3,H3K9ac)
# Silent_histone <-bind_rows(H3K27me3,H3K9me3)


CHROMATIN_Accessibility$Class<-"CHROMATIN Accessibility"
CTCF$Class<-"CTCF"
TFBS$Class<-"TFBS"
H3K27ac$Class<-"H3K27ac"
H3K36me3$Class<-"H3K36me3"
H3K4me1$Class<-"H3K4me1"
H3K4me3$Class<-"H3K4me3"
H3K9ac$Class<-"H3K9ac"
H3K27me3$Class<-"H3K27me3"
H3K9me3$Class<-"H3K9me3"

# cluster_number <-c(0,1,2,3,5,"All")

# count <-c(18739,18477,4617,4429,103,46369)
# cluster <-data.frame(cluster_number,count)
all_anno <-bind_rows(CHROMATIN_Accessibility,CTCF,TFBS,H3K27ac,H3K36me3,H3K4me1,H3K4me3,H3K9ac,H3K27me3,H3K9me3)
# all_anno <-bind_rows(CHROMATIN_Accessibility,TFBS,H3K27ac,H3K36me3,H3K4me1,H3K4me3,H3K9ac,H3K27me3,H3K9me3)
colnames(all_anno)[4] <-"Cluster_number"
all_anno$pos <-paste(all_anno$V1,all_anno$V2,all_anno$V3,sep="_")
anno_count <-all_anno%>%group_by(cluster_number,Class)%>%summarise(factor=n())%>%as.data.frame()

# anno_count_all <-filter(all_anno,cluster_number!=4)%>%group_by(Class)%>%summarise(factor=n())%>%as.data.frame()




# anno_count_all$cluster_number <-"All"
# anno_count_all <- anno_count_all[,c("cluster_number", "Class", "factor")]
anno_count$cluster_number <-as.character(anno_count$cluster_number)
# anno_count <-bind_rows(anno_count_all,anno_count)
final_data<-filter(anno_count,cluster_number!=5)
figure_used <-acast(final_data,Class~cluster_number)
res =chisq.test(figure_used,simulate.p.value = TRUE)
expected = res$expected
roe = as.data.frame(figure_used/expected)
roe$factor <-rownames(roe)
bbb <- melt(roe,id="factor")
colnames(bbb)[2] <-"Cluster"
#-------------------------
bbb$factor<-factor(bbb$factor,levels= c("CHROMATIN Accessibility","TFBS","CTCF","H3K27ac","H3K4me1","H3K4me3","H3K9ac","H3K36me3","H3K27me3","H3K9me3"))
# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/extend_kmer/5kmer/")

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/extend_kmer/5kmer/figure/")
pdf("08_roe.pdf",height = 8, width = 9 )
p1 <-ggplot(bbb,aes(x=Cluster,y=value,fill=factor))+geom_bar(stat="identity",position="dodge")+
    ggtitle("5kmer")+
    theme(plot.title = element_text(hjust = 0.5))+
    ylab("RO/E")


print(p1)
dev.off()


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/extend_kmer/5kmer/figure/")
save(roe,file="5kmer_roe.Rdata")
ROE <-data.frame(roe,stringrAsFactors =T)
ROE$factor <-rownames(ROE)

roe_used_plot <- melt(roe,id=)
aaa <- as.data.frame(t(roe))
aaa$cluster <- rownames(aaa)
bbb <- melt(aaa,id="cluster")



