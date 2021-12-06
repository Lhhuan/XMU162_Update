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


org<-read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176_extend_sorted_merge.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)[1] <-"CHR"
colnames(org)[2] <-"start"
colnames(org)[3] <-"end"
org$hotspot_length <-org$end - org$start
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/figure/")
pdf("041_boxplot_distribution_extend_hotspot_length_log10.pdf",width=3.5, height=3.5)

# p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65,outlier_colour=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Hotspot")+ylab("log10(length of hotspot)")

p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("log10(length of hotspot)")+coord_cartesian(ylim = c(0, 6.5))

print(p)
dev.off()

pdf("041_boxplot_distribution_extend_hotspot_length_log10_no_ylim.pdf",width=3.5, height=3.5)

# p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65,outlier_colour=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Hotspot")+ylab("log10(length of hotspot)")

p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("log10(length of hotspot)")

print(p)
dev.off()
# p <-ggplot(org,aes(x=log10(hotspot_length)))+geom_density(fill="#a3d2ca")


stat <-data.frame(mean = mean(org$hotspot_length),min=min(org$hotspot_length),max=max(org$hotspot_length),sd=sd(org$hotspot_length))%>%t()%>%data.frame()
colnames(stat) <-"length"
stat <-add_column(stat, statistic=rownames(stat), .before = "length")
# stat <-mutate(stat, statistic=rownames(stat), .before = 1)
# # write.table(stat,"")
write.table(stat,"041_distribution_of_extend_hotspot_length.txt",row.names = F, col.names = T,quote =F,sep="\t")



org<-read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org)[1] <-"CHR"
colnames(org)[2] <-"start"
colnames(org)[3] <-"end"
org$hotspot_length <-org$end - org$start
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Tissue_merge/figure/")
pdf("041_boxplot_distribution_original_hotspot_length_log10.pdf",width=3.5, height=3.5)

# p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65,outlier_colour=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Hotspot")+ylab("log10(length of hotspot)")

p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("log10(length of hotspot)")+coord_cartesian(ylim = c(0, 6.5))

print(p)
dev.off()

stat <-data.frame(mean = mean(org$hotspot_length),min=min(org$hotspot_length),max=max(org$hotspot_length),sd=sd(org$hotspot_length))%>%t()%>%data.frame()
colnames(stat) <-"length"
stat <-add_column(stat, statistic=rownames(stat), .before = "length")
# stat <-mutate(stat, statistic=rownames(stat), .before = 1)
# # write.table(stat,"")
write.table(stat,"041_distribution_of_hotspot_length.txt",row.names = F, col.names = T,quote =F,sep="\t")