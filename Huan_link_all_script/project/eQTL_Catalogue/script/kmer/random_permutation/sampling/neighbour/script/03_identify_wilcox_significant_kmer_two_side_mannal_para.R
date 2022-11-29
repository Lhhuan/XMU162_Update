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
# library(data.table)
# setwd("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/")
org<-read.csv("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/hotspot/6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(org)[1] <-"hotspot"
org2 <-melt(org,"hotspot")
colnames(org2)[2] <-"seq"

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/output/")
all_random<-read.csv("/share/data0/QTLbase/huan/eQTL_Catalogue/kmer/random/10_fold_neighbour/100_6mers_uc_us_no_log.csv.gz",header = T,sep = ",") %>% as.data.frame()
colnames(all_random)[1] <-"hotspot"
all_random2 <-melt(all_random,"hotspot")
colnames(all_random2)[2] <-"seq"

unique_hotspot_kmer_need_test<-as.character(unique(org2$seq))

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

re1_1000 <- lapply(1:1000,kmer_test)
d_re1_1000 <- do.call(rbind,re1_1000)
save(d_re1_1000,file="03_d_re1_1000.Rdata")
#=========
re1001_2000 <- lapply(1001:2000,kmer_test)
d_re1001_2000 <- do.call(rbind,re1001_2000)
save(d_re1001_2000,file="03_d_re1001_2000.Rdata")
#=========
re2001_3000 <- lapply(2001:3000,kmer_test)
d_re2001_3000 <- do.call(rbind,re2001_3000)
save(d_re2001_3000,file="03_d_re2001_3000.Rdata")
#=========
re3001_4096 <- lapply(3001:4096,kmer_test)
d_re3001_4096 <- do.call(rbind,re3001_4096)
save(d_re3001_4096,file="03_d_re3001_4096.Rdata")
load("03_d_re1_1000.Rdata")
load("03_d_re2001_3000.Rdata")
load("03_d_re3001_4096.Rdata")

rs <- bind_rows(d_re1_1000,d_re1001_2000,d_re2001_3000,d_re3001_4096)
save(rs,file="03_hotspot_10_neighbour_random_kmer_wilcox_twoside.Rdata")
# write.table(rs,"./figure/hotspot_10_neighbour_random_kmer_wilcox_twoside.txt",row.names = F, col.names = T,quote =F,sep="\t")
write.table(rs,"./figure/03_hotspot_10_neighbour_random_kmer_wilcox_twoside.txt",row.names = F, col.names = T,quote =F,sep="\t")

rs_FDR <-rs
rs_FDR$FDR <-p.adjust(rs_FDR$p_value,method="fdr")
rs_FDR_sig <-filter(rs_FDR,FDR<0.05)
write.table(rs_FDR_sig,"./figure/03_hotspot_10_neighbour_random_kmer_wilcox_twoside_fdr_sig.txt",row.names = F, col.names = T,quote =F,sep="\t")


# fold_c <-data.frame()
# for(kmer in unique_hotspot_kmer_need_test){
# for(i in c(1:length(unique(rs_FDR_sig$seq)))){
fold_change <- function(i){
    kmer = unique(rs_FDR_sig$seq)[i]
    # ProcessBedGz_test <-function(hotspot_kmer_need_test_value =NULL, all_random_kmer_need_test_value =NULL,kmer= NULL){
    hotspot <-filter(org2,seq == kmer)
    hotspot_value <- hotspot$value
    random <-filter(all_random2 ,seq == kmer)
    random_value <- random$value
    fold_change =mean(hotspot_value)/mean(random_value)
    result <-data.frame(seq = kmer,Fold_change=fold_change)
    return(result)
}
f1_900 <-lapply(1:900,fold_change)
df1_900 <-do.call(rbind,f1_900)
save(df1_900,file="03_fold_change_1_900.Rdata")
f901_1800 <-lapply(901:1800,fold_change)
df901_1800 <-do.call(rbind,f901_1800)
save(df901_1800,file="03_fold_change_901_1800.Rdata")
f1801_2700<-lapply(1801:2700,fold_change)
df1801_2700<-do.call(rbind,f1801_2700)
save(df1801_2700,file="03_fold_change_1801_2700.Rdata")
f2701_3245 <-lapply(2701:length(unique(rs_FDR_sig$seq)),fold_change)
df2701_3245<-do.call(rbind,f2701_3245)
# df2701_3245 <-bind_rows(df2701_3000,df3001_3245)
save(df2701_3245,file="03_fold_change_2701_3245.Rdata")
fold_c <-bind_rows(df1_900,df901_1800,df1801_2700,df2701_3245)

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