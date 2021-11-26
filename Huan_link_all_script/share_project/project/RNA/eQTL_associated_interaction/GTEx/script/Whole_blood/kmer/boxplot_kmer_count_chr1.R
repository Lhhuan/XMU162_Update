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

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/chr1/")
org<-read.csv("6mers.csv",header = T,sep = ",") %>% as.data.frame()

colnames(org)[1] <-"hotspot"
org1 <-melt(org,id ="hotspot")
colnames(org1)[2] <-"kmer"
colnames(org1)[3] <-"count"
pdf("kmer_distrbution.pdf")
p<-ggplot(org1,aes(x=kmer,y=count))+geom_boxplot(aes(fill=kmer),outlier.colour = NA)+ theme(legend.position ="none") +scale_y_continuous(limits=0.5)
# +scale_y_continuous(limits=c(0,0.5))

print(p)
dev.off()


org<-read.csv("6mers_non_normalization_.csv",header = T,sep = ",") %>% as.data.frame()

colnames(org)[1] <-"hotspot"
org1 <-melt(org,id ="hotspot")
colnames(org1)[2] <-"kmer"
colnames(org1)[3] <-"count"
pdf("kmer_distrbution_non_normalized.pdf")
p<-ggplot(org1,aes(x=kmer,y=count))+geom_boxplot(aes(fill=kmer),outlier.colour = NA)+ theme(legend.position ="none")+scale_y_continuous(limits=100)
# +scale_y_continuous(limits=c(0,0.5))

print(p)
dev.off()

pdf("kmer_distrbution_non_normalized_no_limit.pdf")
p<-ggplot(org1,aes(x=kmer,y=count))+geom_boxplot(aes(fill=kmer),outlier.colour = NA)+ theme(legend.position ="none")
# +scale_y_continuous(limits=c(0,0.5))

print(p)
dev.off()

pdf("kmer_distrbution_non_normalized_limit0_100.pdf")
p<-ggplot(org1,aes(x=kmer,y=count))+geom_boxplot(aes(fill=kmer),outlier.colour = NA)+ theme(legend.position ="none")+scale_y_continuous(limits=c(0,100))

print(p)
dev.off()




org<-read.csv("6mers_log_non_normalization_uncentered.csv",header = T,sep = ",") %>% as.data.frame()

colnames(org)[1] <-"hotspot"
org1 <-melt(org,id ="hotspot")
colnames(org1)[2] <-"kmer"
colnames(org1)[3] <-"count"
pdf("kmer_distrbution_log_non_normalization_uncentered.pdf")
p<-ggplot(org1,aes(x=kmer,y=count))+geom_boxplot(aes(fill=kmer),outlier.colour = NA)+ theme(legend.position ="none")
# +scale_y_continuous(limits=100)
# +scale_y_continuous(limits=c(0,0.5))

print(p)
dev.off()

pdf("kmer_distrbution_log_non_normalization_uncentered_limit_0_2.5.pdf")
p<-ggplot(org1,aes(x=kmer,y=count))+geom_boxplot(aes(fill=kmer),outlier.colour = NA)+ theme(legend.position ="none")+scale_y_continuous(limits=c(0,2.5))

print(p)
dev.off()



# random<-read.table("original_random_histone_marker.txt.gz",header = T,sep = "\t") %>% as.data.frame()

# hotspot$random_number <- 1
# hotspot$Class <- "Hotspot"
# random$Class <- "Random expectaion"
# # hotspot$percentage <- hotspot$number/46369*100
# # random$percentage <-random$number/46369*100
# rs<- bind_rows(hotspot,random)

# unique_Marker <- unique(rs$Marker)



# p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                                                 panel.background = element_blank(), axis.title.y = element_text(size = 8),
#                                                 # axis.title.x = element_text(size = 10),
#                                                 axis.line = element_line(colour = "black"))



# p <-ggplot(aa,aes(x=Class,y=count))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
#     #scale_y_continuous(limits=c(0,5)) + 
#     theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 7, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Count of segments")+p_theme