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
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/Blood/output/")
org <- fread("eqtl_warm_hotspot_region_win4518_large_than6.bed.gz",header = F,sep = "\t") %>% as.data.frame()
predict <- fread("08_warm_region_predict_key_info.txt",header = T,sep = "\t") %>% as.data.frame()
posi <- filter(predict,predict_class==1)
nega <-filter(predict,predict_class==0)

write.table(posi[,1:3],"10_eqtl_warm_region_positive.bed",row.names = F, col.names = F,quote =F,sep="\t")
write.table(nega[,1:3],"10_eqtl_warm_region_negative.bed",row.names = F, col.names = F,quote =F,sep="\t")



colnames(org) <- c("e_chr","e_start","e_end","pvalue","egene","hchr","hstart","h_end")
org$hotspot <-paste0(org$hchr,":",org$hstart,"-",org$h_end)
org1 <- left_join(predict[,c("hotspot","predict_class")],org,by="hotspot")
org2 <-org1%>%group_by(e_chr,e_start,e_end,hchr,hstart,h_end)%>%summarise(mim_p=min(pvalue))%>%data.frame()

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

rs <-data.frame()
for (cutoff in c(0.05,1e-2, 1e-3,1e-4,1e-5,5e-8)){
    org3 <- filter(org2,mim_p<cutoff)
    org_n <- org3%>%group_by(hchr,hstart,h_end)%>%summarise(number_of_eqtl=n())%>%data.frame()
    colnames(predict)[1:3] <- c("hchr","hstart","h_end")
    predict_n <- left_join(predict[,c(1:3,5)],org_n,by=c("hchr","hstart","h_end") )
    predict_n$number_of_eqtl[is.na(predict_n$number_of_eqtl)]<-0
    predict_n$predict_class <- gsub(0,"False",predict_n$predict_class)
    predict_n$predict_class <- gsub(1,"True",predict_n$predict_class)
    predict_n$length <- predict_n$h_end -predict_n$hstart
    predict_n$adjust_num <-predict_n$number_of_eqtl/predict_n$length*1000
    predict_n$class <-factor(predict_n$predict_class,levels=c("True","False"))
    p1 <- ggplot(predict_n, aes(x=as.factor(class), y=number_of_eqtl,group=class,fill = as.factor(class))) + 
    # p1 <- ggplot(dat, aes(x=as.factor(class), y=adjust_num,group=class,fill = as.factor(class))) + 
    geom_boxplot(outlier.shape = NA)+
    scale_fill_manual(values=c("#A593E0","#F68657"))+
    theme_bw()+
    # labs(x=NULL,y="log(number of eqtl per kb)")+
    labs(x=NULL,y="Number of eQTL")+
    ggtitle(cutoff)+
    theme(legend.position ="none")+
    stat_compare_means(method = "wilcox.test",comparisons =list(c("True","False")),label="p.signif",method.args = list(alternative = "greater"))+
    p_theme 
    pdf(paste0("./figure/10_predicted_hotspot_eQTL_distribution_",cutoff,".pdf"),height=2.5,width=2.6)
    print(p1)
    dev.off()
    #============================
    all_posi <-filter(predict_n,class=="True")
    all_nega <- filter(predict_n,class=="False")
    posi_n <- filter(predict_n,class=="True" & number_of_eqtl >0)
    nega_n <- filter(predict_n,class=="False" & number_of_eqtl >0)
    cover_eqtl=c(nrow(posi_n),nrow(nega_n))
    total=c(nrow(all_posi),nrow(all_nega))
    class=c("True","False")
    prop_t <- data.frame(cover_eqtl=c(nrow(posi_n),nrow(nega_n)),total=c(nrow(all_posi),nrow(all_nega)),class=c("True","False"))
    prop_t$prop <- prop_t$cover_eqtl/prop_t$total
    prop_t$cutoff <- cutoff 
    rs <-bind_rows(rs,prop_t)
    print(cutoff)
    binom.test(nrow(posi_n),nrow(all_posi),p=nrow(nega_n)/nrow(all_nega),alternative="greater")
}

rs$class <-factor(rs$class,levels=c("True","False"))
rs$cutoff <-as.factor(rs$cutoff)
p1 <- ggplot(rs,mapping=aes(x=cutoff,y=prop,fill=class))+geom_bar(stat = 'identity',position=position_dodge(0.9))+ #
    # scale_fill_manual(values=c("#1E77B4","#FF7F0E","#2CA02C","#C22324","#9567BD","#8C554B"))+
    scale_fill_manual(values=c("#A593E0","#F68657"))+
    p_theme+
    labs(x="",y="proportion of segments",fill="")
pdf("./figure/10_segment_proportion.pdf",height=4,width=5)
print(p1)
dev.off()
