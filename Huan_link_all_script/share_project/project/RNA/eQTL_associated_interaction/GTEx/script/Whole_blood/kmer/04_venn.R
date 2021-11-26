library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/factor_anno")
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

Active_histone <-bind_rows(H3K27ac,H3K36me3,H3K4me1,H3K4me3,H3K9ac)
Silent_histone <-bind_rows(H3K27me3,H3K9me3)


CHROMATIN_Accessibility$Clsss <-"CHROMATIN_Accessibility"
CTCF$Clsss <-"CTCF"
TFBS$Clsss <-"TFBS"
Active_histone$Clsss <-"Active_histone"
Silent_histone$Clsss <-"Silent_histone"
all_anno <-bind_rows(CHROMATIN_Accessibility,CTCF,TFBS,Active_histone,Silent_histone)

#-----------------
i=0
aaaaaa <-c(0,1,2,3,5)
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/")
for(i in aaaaaa){
#-------------pie
    a_All <-filter(All,V4==i)%>%dplyr::select(V1,V2,V3)
    a_All$pos <-paste(a_All$V1,a_All$V2,a_All$V3,sep="_")

    a_anno <-filter(all_anno,V4==i)%>%dplyr::select(V1,V2,V3)%>%unique()
    a_anno$pos <-paste(a_anno$V1,a_anno$V2,a_anno$V3,sep="_")
    no_anno <-anti_join(a_All,a_anno,by = "pos")
    a_anno$Class <-"Annotation"
    no_anno$Class <-"Other"
    all_used_Pie <-bind_rows(a_anno,no_anno)%>%dplyr::select(pos,Class)
    a<-table(all_used_Pie$Class)%>%as.data.frame()
    mycolor<-c("#93abd3","#9ad3bc")

    pdf(paste("./figure/Cluater",i,".pdf",sep=""))
    pie(a$Freq, cex=2,col = mycolor,labels = a$Freq, radius = 1,main=paste0("Cluster",i),cex.main=2)
        legend("topright", c("Annotation","Other"), cex = 1, fill = mycolor)
    dev.off()

    #----------veen
    a <-filter(CTCF,V4==i)%>%dplyr::select(V1,V2,V3)
    a_CTCF <-paste(a$V1,a$V2,a$V3,sep="_")

    a <-filter(TFBS,V4==i)%>%dplyr::select(V1,V2,V3)
    a_TFBS<-paste(a$V1,a$V2,a$V3,sep="_")

    a <-filter(Active_histone,V4==i)%>%dplyr::select(V1,V2,V3)
    a_Active_histone <-paste(a$V1,a$V2,a$V3,sep="_")

    a <-filter(Silent_histone,V4==i)%>%dplyr::select(V1,V2,V3)
    a_Silent_histone <-paste(a$V1,a$V2,a$V3,sep="_")

    # setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/figure/")

    library(VennDiagram)

    my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c")
    pdf(paste("./figure/Veen_Cluater",i,".pdf",sep=""))
    T<-venn.diagram(list(CTCF=a_CTCF,TFBS=a_TFBS,"Active histone modification"=a_Active_histone,"Silent histone modification"=a_Silent_histone),filename=NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=1.5,fontface = "bold",cat.cex=1.2,cat.fontfamily = "serif")

    grid.draw(T)
    dev.off()
    print(i)
}
i=4
    a <-filter(CTCF,V4==i)%>%dplyr::select(V1,V2,V3)
    a_CTCF <-paste(a$V1,a$V2,a$V3,sep="_")

    a <-filter(TFBS,V4==i)%>%dplyr::select(V1,V2,V3)
    a_TFBS<-paste(a$V1,a$V2,a$V3,sep="_")

    a <-filter(Active_histone,V4==i)%>%dplyr::select(V1,V2,V3)
    a_Active_histone <-paste(a$V1,a$V2,a$V3,sep="_")

    a <-filter(Silent_histone,V4==i)%>%dplyr::select(V1,V2,V3)
    a_Silent_histone <-paste(a$V1,a$V2,a$V3,sep="_")

    # setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/ALL/figure/")

    library(VennDiagram)

    my_c2=c(colors()[616], colors()[38], colors()[468],"#beca5c")
    pdf(paste("./figure/Veen_Cluater",i,".pdf",sep=""))
    T<-venn.diagram(list(CTCF=a_CTCF,TFBS=a_TFBS,"Active histone modification"=a_Active_histone,"Silent histone modification"=a_Silent_histone),filename=NULL,lwd=1,lty=2,col=my_c2 ,fill=my_c2,reverse=TRUE,cex=1.5,fontface = "bold",cat.cex=1.2,cat.fontfamily = "serif")

    grid.draw(T)
    dev.off()
    print(i)