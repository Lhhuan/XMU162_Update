library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

#---------------------------------------------------------
library(Hmisc)
# tissue <-"Whole_Blood"
# cutoff=0.176
# file_name<-paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_",cutoff,".bed.gz")




setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_11/6kmer/5_community/")



org<-read.table("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/communities_5.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)[1] <-"CHR"
colnames(org)[2] <-"start"
colnames(org)[3] <-"end"
colnames(org)[4] <-"Community"
org$hotspot_length <-org$end - org$start
# setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/figure/")
pdf("12_boxplot_distribution_extend_hotspot_length_log10.pdf",width=3.5, height=3.5)

# p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65,outlier_colour=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Hotspot")+ylab("log10(length of hotspot)")

p<-ggplot(org,aes(x=factor(Community), y=log10(hotspot_length),fill=Community))+geom_violin(width=0.65)+geom_boxplot(width=0.1,outlier.color=NA)+ theme(legend.position ="none")+p_theme+xlab("")+ylab("log10(length of hotspot)")+coord_cartesian(ylim = c(0, 6.5))

print(p)
dev.off()
#-----------------------------

plot<-function(i){
    setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/chr1_11/6kmer/",i,"_community/"))
    org<-read.table(paste0("/share/data0/QTLbase/huan/GTEx/Tissue_merge/Cis_eQTL/hotspot/interval_18/chr1_11/kmer/6/communities_",i,".bed.gz"),header = F,sep = "\t") %>% as.data.frame()
    colnames(org) <-c("CHR","start","end","Community")
    org$hotspot_length <-org$end - org$start
    # setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/figure/")
    pdf("12_boxplot_distribution_extend_hotspot_length_log10.pdf",width=3.5, height=3.5)

    # p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65,outlier_colour=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Hotspot")+ylab("log10(length of hotspot)")

    p<-ggplot(org,aes(x=factor(Community), y=log10(hotspot_length),fill=Community))+geom_violin(width=0.65)+geom_boxplot(width=0.1,outlier.color=NA)+ theme(legend.position ="none")+p_theme+xlab("")+ylab("log10(length of hotspot)")+coord_cartesian(ylim = c(0, 6.5))

    print(p)
    dev.off()
    #-----------------------------

}
lapply(c(3:9),plot)