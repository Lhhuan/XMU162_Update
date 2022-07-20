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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/")
load("06_wilcox_permutation_overlap_significant.Rdata")

permutation_sig <- wilcox_permutation_overlap
permutation_sig <-permutation_sig[order(permutation_sig$p_value),]
tsg <-permutation_sig[1:30,]
load("all_random_kmer_need_test_ratio_1_100.Rdata")
org1 <- read.table("../../../../figure/hit_hospot_ratio_kmer.txt",header = T,sep = "\t") %>% as.data.frame()
# kmer_need <-filter(org1,ratio>0.5)
kmer_need <-org1

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), 
                                                axis.title.y = element_text(size = 8,color = "black"),
                                                axis.line = element_line(colour = "black"),
                                                plot.title = element_text(hjust = 0.5,size = 8),
                                                axis.title.x=element_blank(),
                                                axis.text.x = element_text(size = 8, color = "black"),
                                                axis.text.y= element_text(size = 6.5))
#------------------------
random <-random[,-3]

unique_hotspot_kmer_need_test1 <-as.character(unique(tsg$seq))
# for(kmer in unique_hotspot_kmer_need_test1){
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure")
plot_dis <-function(kmer){
    # ProcessBedGz_test <-function(hotspot_kmer_need_test_value =NULL, all_random_kmer_need_test_value =NULL,kmer= NULL){
    hotspot <-filter(kmer_need,seq == kmer)
    random <-filter(random,seq == kmer)
    hotspot$class <- "Hotspot"
    random$class <- "Random"
    dat1 <-bind_rows(hotspot,random)
    title_name<-kmer
    p <-ggplot(dat1,aes(x=class,y=ratio))+geom_boxplot(aes(fill=class),width=0.3,outlier.colour = NA)+ 
    # scale_y_continuous(limits=c(0,10)) + 
    theme(legend.position ="none")+ggtitle(title_name) +ylab("Fraction of segments")+p_theme
    p1 <- add_pval(p,annotation = "****",pval_star = T)
    fig_name=paste0("./kmer_fig/",kmer,".pdf")
    pdf(fig_name,width=1.5, height=1.5)
    print(p1)
    dev.off()
    return(p1)
    print(kmer)
}

plist = lapply(unique_hotspot_kmer_need_test1,plot_dis)

pdf("./07_sig_permutation_overlap_kmer_distrbution.pdf",width=12, height=10)
CombinePlots(plist,ncol=6,nrow=5)
dev.off()
print(111)

pdf("./07_sig_permutation_overlap_kmer_distrbution_10.pdf",width=8, height=4.1)
gridExtra::marrangeGrob(plist,nrow=2,ncol=4)
dev.off()


pdf("./07_sig_permutation_overlap_kmer_distrbution_12.pdf",width=8, height=6.1)
gridExtra::marrangeGrob(plist,nrow=3,ncol=4)
dev.off()
print(111)





all_s_h <-mclapply(c(1:300), ProcessBedGz, mc.cores = 40)