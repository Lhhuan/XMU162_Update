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
library(parallel)
library(data.table)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/")
org1 <- read.table("../../../figure/hit_hospot_ratio_kmer.txt",header= T,sep = "\t") %>% as.data.frame()
# kmer_need <-filter(org1,ratio>0.5)
wilcox <-read.table("./output/figure/03_hotspot_10_neighbour_random_kmer_wilcox_twoside_fdr_sig_foldchange.txt",header= T,sep = "\t") %>% as.data.frame()
sigK <- filter(wilcox,Fold_change>=1.15 | Fold_change <=0.9)
kmer_need <-filter(org1,seq %in%sigK$seq)

ProcessBedGz<-function(i=NULL){
    setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/random/10_fold_neighbour/")
    print(i)
    file_name <-paste0(i,"_6mers_uc_us_no_log.csv.gz")
    org<-fread(file_name,header = T,sep = ",") %>% as.data.frame()
    rownames(org) <- org[,1]
    Sorg <-org[,which(colnames(org) %in% sigK$seq)]
    Sorg$hotspot <- rownames(Sorg)
    org2 <-reshape2::melt(Sorg,"hotspot")
    colnames(org2)[2] <-"seq"
    # sig_hotspot <-filter(org2, value>0 & seq %in% sigK$seq)
    sig_hotspot <-filter(org2, value>0)
    seq_count_in <- sig_hotspot%>%group_by(seq)%>%dplyr::summarise(count=n())%>%as.data.frame()
    # seq_count_in <-filter(seq_count,seq %in% org1$seq) 
    gc()
    # rs<-bind_rows(rs,sig_hotspot)
    print(paste(i,"end",sep=" "))
    seq_count_in$random_ID=i
    return(seq_count_in)
}

# all_s_h <-mclapply(c(1:100), ProcessBedGz, mc.cores = 20)
# all_s_h <-mclapply(c(1:1000), ProcessBedGz, mc.cores = 20)
# all_6kmers <-do.call(rbind,all_s_h)
# a <-lapply(c(1:1000), ProcessBedGz)
a <-lapply(c(1:1000), ProcessBedGz)
random_more0_count<-do.call(rbind,a)

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/")
hotf<-read.table("/home/huanhuan/project/eQTL_Catalogue/output/all_tissue_status/hotspot/Tissue_merge_segment_hotspot_cutoff_0.176_extend_18_snp_sorted_merge.bed.gz",header = F,sep = "\t")%>% as.data.frame()

random_more0_count$ratio <-random_more0_count$count/nrow(hotf)
random <-random_more0_count
save(random,file="06_random_kmer_need_test_ratio_1_1000.Rdata")
colnames(kmer_need)[2:3] <-c("hotspot_count","hotspot_ratio")
dat <-inner_join(random,kmer_need,by="seq")
#-------------hotspot higher than random
dat$diff <- dat$ratio -dat$hotspot_ratio 
dat$diff0 <-NA
dat[which(dat$diff>=0),"diff0"]=1
dat[which(dat$diff<0),"diff0"]=0
fdat_higher <-dat%>%group_by(seq)%>%summarise(n=sum(diff0))%>%data.frame()
fdat_higher$p_value <-fdat_higher$n/1000
fdat_higher$class <-"hotspot_higher"
#------------------------hotspot lower than random
dat$diff2 <- dat$hotspot_ratio - dat$ratio 
dat$diff02 <-NA
dat[which(dat$diff2>=0),"diff02"]=1
dat[which(dat$diff2<0),"diff02"]=0
fdat_lower <-dat%>%group_by(seq)%>%summarise(n=sum(diff02))%>%data.frame()
fdat_lower$p_value <-fdat_lower$n/1000
fdat_lower$class <-"hotspot_lower"

fdat <-bind_rows(fdat_higher,fdat_lower)
fdat$FDR <-p.adjust(fdat$p_value,method="fdr")

fdat0 <-filter(fdat,FDR<0.05)
save(fdat0,file="06_permutation_test_1000_sig_kmer.Rdata")
# neighbor <- fdat0 #1835
# #-----------------------------





# load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/emplambda_0/06_permutation_test_100_sig_kmer.Rdata")
# emp0<- fdat0

# overlap_neig_emp0 <-filter(neighbor,seq %in% emp0$seq)


# wilcox_sig <- read.table("./figure/03_hotspot_10_neighbour_random_kmer_wilcox_twoside.txt",header = T,sep = "\t") %>% as.data.frame()
# wilcox_sig <-filter(wilcox_sig,p_value<0.05)
# save(wilcox_sig,file="06_wilcox_sig_two_side.Rdata")
# wilcox_permutation_overlap <-filter(neighbor,seq %in% wilcox_sig$seq)
# wilcox_emp0_overlap <- filter(wilcox_sig,seq %in% emp0$seq)
# save(wilcox_permutation_overlap,file="06_wilcox_permutation_overlap_significant_twoside.Rdata")

# #------------------------------------------------------------------------------------------------------- 
# library(VennDiagram)

# # A = 1:150
# # B = c(121:170,300:320)
# pdf("./figure/06_venn_permutation_wilcox_two_side.pdf",width=5,height=5)
# T<-venn.diagram(list(Wilcox=wilcox_sig$seq,Permutation=neighbor$seq),
#                 filename=NULL,
#                 lwd=1,lty=2,
#                 col=c('red','green'),
#                 fill=c('red','green'),
#                 cat.col=c('black','black'),
#                 reverse=TRUE)
# grid.draw(T)
# dev.off()

# overlao_sig_dis <-filter(org1,seq %in%wilcox_permutation_overlap$seq)

# p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                                                 panel.background = element_blank(), axis.title.y = element_text(size = 8),
#                                                 # axis.title.x = element_text(size = 10),
#                                                 axis.line = element_line(colour = "black"))

# pdf("./figure/06_sig_permutation_wilcox_cover_hotspot.pdf",width=3.5, height=3.5)
# p<-ggplot(overlao_sig_dis,aes(x=1, y=ratio))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Kmer")+ylab("Fraction of segments")

# print(p)
# dev.off()



