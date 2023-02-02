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
library(R.utils)
library(mclust)
library(umap)
library(uwot)
library(kernlab)
library(data.table)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/output/")
posi <- fread("081_eqtl_warm_region_predict_hotspot_true.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(posi) <- c("e_chr","e_start","e_end","egene","pvalue","tissue1","tissue2","hchr","hstart","h_end")
posi1 <-posi%>%group_by(e_chr,e_start,e_end,hchr,hstart,h_end)%>%summarise(mim_p=min(pvalue))%>%data.frame()
# posi <-posi[,c("e_chr","e_start","e_end","hchr","hstart","h_end")]%>%unique()

nega <- fread("081_eqtl_warm_region_predict_hotspot_false.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(nega) <- c("e_chr","e_start","e_end","egene","pvalue","tissue1","tissue2","hchr","hstart","h_end")
nega1 <-nega%>%group_by(e_chr,e_start,e_end,hchr,hstart,h_end)%>%summarise(mim_p=min(pvalue))%>%data.frame()
#==========================

    p_theme<-theme(
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))


#===========
# cutoff=1e-3
rs <-data.frame()
# filter_cutoff <- function(){
for (cutoff in c(0.05,1e-2, 1e-3,1e-4,1e-5,5e-8)){
    posi2 <- filter(posi1,mim_p<cutoff)
    posi_n <- posi2%>%group_by(hchr,hstart,h_end)%>%summarise(number_of_eqtl=n())%>%data.frame()
    all_posi <- fread("081_warm_region_predict_hotspot_true.bed",header = F,sep = "\t") %>% as.data.frame()
    colnames(all_posi) <- c("hchr","hstart","h_end")
    all_posi_n <- left_join(all_posi,posi_n,by=c("hchr","hstart","h_end") )
    all_posi_n$number_of_eqtl[is.na(all_posi_n$number_of_eqtl)]<-0
    all_posi_n$class <- "True"
    #=====================
    nega2 <- filter(nega1,mim_p<cutoff)
    nega_n <- nega2%>%group_by(hchr,hstart,h_end)%>%summarise(number_of_eqtl=n())
    all_nega <- fread("081_warm_region_predict_hotspot_false.bed",header = F,sep = "\t") %>% as.data.frame()
    colnames(all_nega) <- c("hchr","hstart","h_end")
    all_nega_n <- left_join(all_nega,nega_n,by=c("hchr","hstart","h_end") )
    all_nega_n$number_of_eqtl[is.na(all_nega_n$number_of_eqtl)]<-0
    all_nega_n$class <- "False"
    #=================
    dat <- bind_rows(all_posi_n,all_nega_n)
    dat$length <- dat$h_end -dat$hstart
    dat$adjust_num <-dat$number_of_eqtl/dat$length*1000
    binom.test(nrow(posi_n),nrow(all_posi),p=nrow(nega_n)/nrow(all_nega),alternative="greater")
    dat$class <-factor(dat$class,levels=c("True","False"))
    # p1 <- ggplot(dat, aes(x=as.factor(class), y=log(adjust_num+0.00000001),group=class,fill = as.factor(class))) + 
    p1 <- ggplot(dat, aes(x=as.factor(class), y=number_of_eqtl,group=class,fill = as.factor(class))) + 
    # p1 <- ggplot(dat, aes(x=as.factor(class), y=adjust_num,group=class,fill = as.factor(class))) + 
    geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values=c("#A593E0","#F68657"))+
    theme_bw()+
    # labs(x=NULL,y="log(number of eqtl per kb)")+
    labs(x=NULL,y="Number of eQTL)")+
    ggtitle(cutoff)+
    theme(legend.position ="none")+
    # stat_compare_means(method = "wilcox.test",label="p.signif")+
    stat_compare_means(method = "wilcox.test",comparisons =list(c("True","False")),label="p.signif",method.args = list(alternative = "greater"))+
    p_theme #+
    # ylim(-5,4)

    pdf(paste0("./figure/082_predicted_hotspot_eQTL_distribution_",cutoff,".pdf"),height=2.5,width=2.6)
    # pdf(paste0("./figure/082_predicted_hotspot_eQTL_distribution_",cutoff,".pdf"),height=2.5,width=2.6)
    print(p1)
    dev.off()
    cover_eqtl=c(nrow(posi_n),nrow(nega_n))
    total=c(nrow(all_posi),nrow(all_nega))
    class=c("True","False")
    prop_t <- data.frame(cover_eqtl=c(nrow(posi_n),nrow(nega_n)),total=c(nrow(all_posi),nrow(all_nega)),class=c("True","False"))
    prop_t$prop <- prop_t$cover_eqtl/prop_t$total
    prop_t$cutoff <- cutoff 
    rs <-bind_rows(rs,prop_t)
    print(cutoff)
}

rs$class <-factor(rs$class,levels=c("True","False"))
rs$cutoff <-as.factor(rs$cutoff)
p1 <- ggplot(rs,mapping=aes(x=cutoff,y=prop,fill=class))+geom_bar(stat = 'identity',position=position_dodge(0.9))+ #
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    scale_fill_manual(values=c("#A593E0","#F68657"))+
    p_theme+
    labs(x="",y="proportion of segments",fill="")
pdf("./figure/082_segment_proportion.pdf",height=4,width=5)
print(p1)
dev.off()