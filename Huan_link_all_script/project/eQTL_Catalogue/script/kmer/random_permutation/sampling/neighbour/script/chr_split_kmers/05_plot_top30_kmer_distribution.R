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

#---------
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/")
dat <- read.table("./03_hotspot_10_neighbour_random_kmer_wilcox_twoside_fdr_sig_foldchange.txt",header = T,sep = "\t") %>% as.data.frame()
sdat <- filter(dat,Fold_change>=1.15 | Fold_change <=0.9)
sdat <-sdat[order(sdat$FDR),]
tsg <-sdat[1:30,]
save(tsg,file="05_wilcox_plot_kmer.Rdata")

org<-read.csv("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(org)[1] <-"hotspot"
org2 <-melt(org,"hotspot")
colnames(org2)[2] <-"seq"

sorg <-filter(org2,seq %in% tsg$seq)

all_random<-read.csv("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/random/10_fold_neighbour/1_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(all_random)[1] <-"hotspot"
all_random2 <-melt(all_random,"hotspot")
colnames(all_random2)[2] <-"seq"
srandom <- filter(all_random2,seq %in% tsg$seq)
#--------------------------------------------------------------------------------------------------
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/figure/")
unique_hotspot_kmer_need_test1 <-as.character(unique(tsg$seq))

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), 
                                                axis.title.y = element_text(size = 8,color = "black"),
                                                axis.line = element_line(colour = "black"),
                                                plot.title = element_text(hjust = 0.5,size = 8),
                                                axis.title.x=element_blank(),
                                                axis.text.x = element_text(size = 8, color = "black"),
                                                axis.text.y= element_text(size = 6.5))

# for(kmer in unique_hotspot_kmer_need_test1){
plot_dis <-function(kmer){
    # ProcessBedGz_test <-function(hotspot_kmer_need_test_value =NULL, all_random_kmer_need_test_value =NULL,kmer= NULL){
    hotspot <-filter(sorg,seq == kmer)
    random <-filter(srandom,seq == kmer)
    hotspot$class <- "Hotspot"
    random$class <- "Random"
    dat <-bind_rows(hotspot,random)
    title_name<-kmer
    p <-ggplot(dat,aes(x=class,y=log(value)))+geom_boxplot(aes(fill=class),width=0.3,outlier.colour = NA)+ 
    # scale_y_continuous(limits=c(0,10)) + 
    theme(legend.position ="none")+ggtitle(title_name) +ylab("Log(adjust count)")+p_theme
    p1 <- add_pval(p,annotation = "****",pval_star = T)
    fig_name=paste0("./kmer_fig/",kmer,".pdf")
    pdf(fig_name,width=1.5, height=1.5)
    print(p1)
    dev.off()
    return(p1)
    print(kmer)
}

plist = lapply(unique_hotspot_kmer_need_test1,plot_dis)

pdf("./05_sig_wilcox_kmer_distrbution.pdf",width=12, height=10)
CombinePlots(plist,ncol=6,nrow=5)
dev.off()
print(111)

pdf("./05_sig_wilcox_kmer_distrbution_8.pdf",width=8, height=4.1)
gridExtra::marrangeGrob(plist,nrow=2,ncol=4)
dev.off()
print(111)

pdf("./05_sig_wilcox_kmer_distrbution_12.pdf",width=8, height=6.1)
gridExtra::marrangeGrob(plist,nrow=3,ncol=4)
dev.off()
print(111)


pdf("./05_sig_wilcox_kmer_distrbution.pdf",width=12, height=10)
CombinePlots(plist,ncol=6,nrow=5)
dev.off()
print(111)




# load("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/emplambda_0/06_permutation_test_100_sig_kmer.Rdata")

# overlap <- filter(sdat,seq%in%fdat0$seq)




