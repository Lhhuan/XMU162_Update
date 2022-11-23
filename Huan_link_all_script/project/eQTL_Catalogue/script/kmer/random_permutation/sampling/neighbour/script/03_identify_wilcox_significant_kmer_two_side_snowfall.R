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
library(snowfall)
# setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/")
org<-read.csv("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(org)[1] <-"hotspot"
org2 <-melt(org,"hotspot")
colnames(org2)[2] <-"seq"

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/")
all_random<-read.csv("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/random/10_fold_neighbour/1_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(all_random)[1] <-"hotspot"
all_random2 <-melt(all_random,"hotspot")
colnames(all_random2)[2] <-"seq"

unique_hotspot_kmer_need_test<-as.character(unique(org2$seq))
# rs <-data.frame()
# for(kmer in unique_hotspot_kmer_need_test){
# for(i in c(1:length(unique_hotspot_kmer_need_test))){
sfInit(parallel = TRUE, cpus = 30)
kmer_test<-function(i){
    kmer = unique_hotspot_kmer_need_test[i]
    hotspot <-filter(org2,seq == kmer)
    hotspot_value <- hotspot$value
    random <-filter(all_random2 ,seq == kmer)
    random_value <- random$value
    re <- wilcox.test(hotspot_value, random_value,alternative = "two.sided")
    p_value = re$p.value
    #---------------------
    if(p_value <= 0.0001){
        annotation = "****"
    }else if(p_value <= 0.001){
        annotation = "***"
    }else if(p_value <= 0.01){
        annotation = "**"
    }else if(p_value <= 0.05){
        annotation = "*"
    }else{
        annotation = "ns"
    }
    #-------------
    W <-re$statistic
    names(W)=NULL
    w_re <-data.frame(seq = kmer,W=W,p_value = re$p.value,annotation=annotation)
    return(w_re)
    # rs <-bind_rows(rs,w_re)
    print(i)
}
#------------------------
sfLibrary(dplyr)
sfExport("unique_hotspot_kmer_need_test")
sfExport("org2")
sfExport("all_random2")
re1 <- sfLapply(1:length(unique_hotspot_kmer_need_test),kmer_test)

sfStop()


# write.table(rs,"./figure/hotspot_10_neighbour_random_kmer_wilcox_twoside.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(rs,"./figure/03_hotspot_10_neighbour_random_kmer_wilcox_twoside.txt",row.names = F, col.names = T,quote =F,sep="\t")

rs_FDR <-rs
rs_FDR$FDR <-p.adjust(rs_FDR$p_value,method="fdr")
rs_FDR_sig <-filter(rs_FDR,FDR<0.05)
write.table(rs_FDR_sig,"./figure/03_hotspot_10_neighbour_random_kmer_wilcox_twoside_fdr_sig.txt",row.names = F, col.names = T,quote =F,sep="\t")


fold_c <-data.frame()
# for(kmer in unique_hotspot_kmer_need_test){
for(i in c(1:length(unique(rs_FDR_sig$seq)))){
    kmer = unique(rs_FDR_sig$seq)[i]
    # ProcessBedGz_test <-function(hotspot_kmer_need_test_value =NULL, all_random_kmer_need_test_value =NULL,kmer= NULL){
    hotspot <-filter(org2,seq == kmer)
    hotspot_value <- hotspot$value
    random <-filter(all_random2 ,seq == kmer)
    random_value <- random$value
    fold_change =mean(hotspot_value)/mean(random_value)
    result <-data.frame(seq = kmer,Fold_change=fold_change)
    fold_c <-bind_rows(fold_c,result)
    print(i)
    # gc()
    # return(w_re)
}

final_re <-left_join(rs_FDR_sig,fold_c,by="seq")
write.table(final_re,"./figure/03_hotspot_10_neighbour_random_kmer_wilcox_twoside_fdr_sig_foldchange.txt",row.names = F, col.names = T,quote =F,sep="\t")

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))



pdf("./figure/03_boxplot_wilcox_sig_fold_change.pdf",width=3.5, height=3.5)
p<-ggplot(final_re,aes(x=1, y=Fold_change))+
geom_violin(fill="#a3d2ca",width=0.65)+
geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ 
p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position ="none")+
xlab("")+ylab("Fold change")+scale_y_continuous( limits=c(0.8, 1.45),breaks=seq(0.8,1.45,0.05))

print(p)
dev.off()

sigK <- filter(final_re,Fold_change>=1.15 | Fold_change <=0.9)