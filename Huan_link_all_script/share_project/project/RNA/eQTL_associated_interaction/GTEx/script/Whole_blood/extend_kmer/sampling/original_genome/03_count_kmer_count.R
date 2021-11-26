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

setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/random/original_genome/fa")
# rs<-data.frame()
# for(i in c(1:1000)){
#     print(i)
#     file_name <-paste0(i,"_6mers_uc_us_no_log.csv.gz")
#     org<-read.csv(file_name,header = T,sep = ",") %>% as.data.frame()
#     colnames(org)[1] <-"hotspot"
#     org2 <-melt(org,"hotspot")
#     colnames(org2)[2] <-"seq"
#     sig_hotspot <-filter(org2, value>0)
#     rs<-bind_rows(rs,sig_hotspot)
#     # print(i)
# }

ProcessBedGz<-function(i=NULL){
    setwd("/share/data0/QTLbase/huan/GTEx/Whole_Blood/Cis_eQTL/hotspot_cis_eQTL/interval_18_filter/6/kmer/extend/random/original_genome/fa")
    print(i)
    file_name <-paste0(i,"_6mers_uc_us_no_log.csv.gz")
    org<-read.csv(file_name,header = T,sep = ",") %>% as.data.frame()
    colnames(org)[1] <-"hotspot"
    org2 <-melt(org,"hotspot")
    colnames(org2)[2] <-"seq"
    sig_hotspot <-filter(org2, value>0)
    seq_count <- sig_hotspot%>%group_by(seq)%>%dplyr::summarise(count=n())%>%as.data.frame()
    gc()
    # rs<-bind_rows(rs,sig_hotspot)
    # print(i)
    return(seq_count)
}

# all_s_h <-mclapply(c(1:10), ProcessBedGz, mc.cores = 20)
# all_s_h <-mclapply(c(1:1000), ProcessBedGz, mc.cores = 20)
# all_6kmers <-do.call(rbind,all_s_h)

a <-mclapply(c(1:1000), ProcessBedGz, mc.cores = 10)
b <-do.call(rbind,a)
final_seq_count <- b%>%group_by(seq)%>%dplyr::summarise(count=sum(count))%>%as.data.frame()

final_seq_count$ratio <-final_seq_count$count/(46369*1000)

seq_count <-final_seq_count
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 8),
                                                # axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))



#-------------------------------------

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/extend_kmer/sampling/original_genome/")

#----------------------------point plot
png("./figure/original_random_Kmer_occur_count.png")

p1 <-ggplot(seq_count,mapping = aes(x=seq ,y = count))+geom_point(size=0.1) +ylab("Count")+ 
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Kmer")+p_theme
print(p1)
dev.off()

png("./figure/original_random_Kmer_occur_ratio.png")
p2 <-ggplot(seq_count,mapping = aes(x=seq ,y = ratio))+geom_point(size=0.1) +ylab("Ratio")+ 
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Kmer") +p_theme

print(p2)
dev.off()

#------------summary boxplot
pdf("./figure/original_random_Kmer_occur_boxplot_ratio.pdf",width=3.5, height=3.5)
p<-ggplot(seq_count,aes(x=1, y=ratio))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.15)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("Kmer")+ylab("Ratio")

print(p)
dev.off()

#-------------------


rs <-data.frame()

cutoffs<-seq(0,0.25,0.01)

for(cutoff in cutoffs){
    aa <-filter(seq_count,ratio >cutoff)
    kmer_ratio =nrow(aa)/(4^6)
    kmer_count = nrow(aa)
    tmp <-data.frame(hotspot_ratio =cutoff,kmer_ratio = kmer_ratio,kmer_count=kmer_count)
    rs<-bind_rows(rs,tmp)
    print(cutoff)
}


pdf("./figure/original_random_hit_random_ratio_Kmer_ratio.pdf",width=3.2, height=3)

p1 <-ggplot(rs,mapping = aes(x=hotspot_ratio ,y = kmer_ratio))+geom_point(size=0.1) +ylab("Number of kmer")+ 
    xlab("Hit random ratio")+p_theme
print(p1)
dev.off()
#-------------------
pdf("./figure/original_random_hit_random_ratio_Kmer_count.pdf",width=3.2, height=3)

p1 <-ggplot(rs,mapping = aes(x=hotspot_ratio ,y = kmer_count))+geom_point(size=0.1) +ylab("Number of kmer")+ 
    xlab("Hit random ratio")+p_theme
print(p1)
dev.off()

#----------

save(seq_count,file="1000.Rdata")