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

load("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/figure/hotspot_kmer_need_test_value.Rdata")
#---------
setwd("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/random_permutation/sampling/emplambda_0/figure/")
org <- read.table(paste0("hotspot_emplambda_0_random_kmer_wilcox.txt"),header = T,sep = "\t") %>% as.data.frame()
org1 <-org[order(org$p_value),]
tsg <-org1[1:30,]
hotspot_kmer <-filter(hotspot_kmer_need_test_value,seq %in% tsg$seq)
load("all_random_kmer_need_test_value_1_100.Rdata")
random_kmer <-filter(all_random_kmer_need_test_value,seq %in% tsg$seq)

#------------------------
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), 
                                                axis.title.y = element_text(size = 3,color = "black"),
                                                axis.line = element_line(colour = "black"),
                                                plot.title = element_text(hjust = 0.5,size = 4),
                                                axis.title.x=element_blank(),
                                                axis.text.x = element_text(size = 3, color = "black"),
                                                axis.text.y= element_text(size = 3))
#------------------------
unique_hotspot_kmer_need_test1 <-as.character(unique(tsg$seq))
rs <-data.frame()
# for(kmer in unique_hotspot_kmer_need_test1){
plot_dis <-function(kmer){
    # ProcessBedGz_test <-function(hotspot_kmer_need_test_value =NULL, all_random_kmer_need_test_value =NULL,kmer= NULL){
    hotspot <-filter(hotspot_kmer,seq == kmer)
    random <-filter(random_kmer,seq == kmer)
    hotspot$class <- "Hotspot"
    random$class <- "Random"
    dat <-bind_rows(hotspot,random)
    title_name<-kmer
    p <-ggplot(dat,aes(x=class,y=value))+geom_boxplot(aes(fill=class),width=0.3,outlier.colour = NA)+ 
    scale_y_continuous(limits=c(0,10)) + 
    theme(legend.position ="none")+ggtitle(title_name) +ylab("Adjust count")+p_theme
    p1 <- add_pval(p,annotation = "****",pval_star = T)
    fig_name=paste0("./kmer_fig/",kmer,".pdf")
    pdf(fig_name,width=1.5, height=1.5)
    print(p1)
    dev.off()
    return(p1)
    print(kmer)
}

plist = lapply(unique_hotspot_kmer_need_test1,plot_dis)

pdf("./sig_kmer_distrbution.pdf",width=12, height=10)
CombinePlots(plist,ncol=6,nrow=5)
dev.off()
print(111)