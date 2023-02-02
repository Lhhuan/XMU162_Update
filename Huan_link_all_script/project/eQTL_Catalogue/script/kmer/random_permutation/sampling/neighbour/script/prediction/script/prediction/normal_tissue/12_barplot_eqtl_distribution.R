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
setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/normal_tissue/")

f1 <- read.table("./output/02_tissue_level_Marker_source.txt",header = T,sep = "\t") %>% as.data.frame()
f1 <-filter(f1,marker!="HISTONE_MARK_AND_VARIANT")
tissues <-sort(unique(f1$refine_tissue_label2))
p_theme<-theme(
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(), 
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5,size=11))

plist <- list()
all_rs <-data.frame()
for (i in 1:length(tissues)){
# for (i in 1:2){
    tissue = tissues[i]
    path = paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/",tissue)
    setwd(path)
    predict <- fread("predicted_key_info.txt.gz",header =T,sep = "\t") %>% as.data.frame()
    predict_n <- data.frame(table(predict$predict_class))
    org <- fread("11_predicted_region_eqtl.bed.gz",header =F,sep = "\t") %>% as.data.frame()
    # colnames(org) <- c("e_chr","e_start","e_end","pvalue","egene","hchr","hstart","h_end")
    colnames(org) <- c("e_chr","e_start","e_end","pvalue","ENSG","hchr","hstart","h_end")
    org1 <-org%>%group_by(e_chr,e_start,e_end,hchr,hstart,h_end)%>%summarise(mim_p=min(pvalue))%>%data.frame()  
    org1$hotspot <-paste0(org1$hchr,":",org1$hstart,"-",org1$h_end)
    org2 <- left_join(org1,predict[,c("hotspot","predict_class")],by="hotspot")
    rs <-data.frame()
    for (cutoff in c(0.05,1e-2, 1e-3,1e-4,1e-5,5e-8)){
        org3 <- filter(org2,mim_p<cutoff)%>%select(hotspot,predict_class)%>%unique()
        org4 <-data.frame(table(org3$predict_class))
        org5 <- left_join(org4,predict_n,by="Var1")
        colnames(org5) <-c("class","N_cutoff","N_total")
        org5$class <-as.character(org5$class)
        org5[4,] <- data.frame("1_2",org5[2,2]+org5[3,2],org5[2,3]+org5[3,3])
        org5$cutoff <- cutoff
        rs <-bind_rows(rs,org5)
        print(cutoff)
    }
    rs$tissue <- tissue
    # rs$class <- paste0("C",rs$class)
    rs$cutoff <-as.factor(rs$cutoff)
    rs$prop <- rs$N_cutoff/rs$N_total
    all_rs <-bind_rows(all_rs,rs)
    #=============================
    new_dir <- paste0("/share/data0/QTLbase/huan/eQTL_Catalogue/prediction/normal_tissue/",tissue,"/eQTL/figure")
    if(dir.exists(new_dir)){
        print(c(new_dir, "exists"))
    }else{
        dir.create(new_dir, recursive = TRUE)
    }
    setwd(new_dir)
    # rs$class <-factor(rs$class,levels=c("C0","C1_2","C1","C2"))
    tissue <-capitalize(tissue)
    rs$class <-factor(rs$class,levels=c("0","1_2","1","2"))
    plist[[i]] <- ggplot(rs,mapping=aes(x=cutoff,y=prop,fill=class))+geom_bar(stat = 'identity',position=position_dodge(0.9))+ #
    scale_fill_manual(values=c("#A593E0","#F68657","#84B1ED","#F6B352"))+
    p_theme+
    ggtitle(tissue)+
    labs(x="",y="proportion of segments",fill="")
    pdf("12_eqtl_segment_proportion.pdf",height=4,width=5)
    print(plist[[i]])
    dev.off()
    #=============================
    print(c(i,tissue))
}


setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/normal_tissue/output/figure/")
write.table(all_rs,"../12_All_tissue_eQTL_proportion.txt",row.names = F, col.names = T,quote =F,sep="\t")
pdf("12_ALL_tissue_eqtl_segment_proportion.pdf",width=15, height=8)
p2<-gridExtra::marrangeGrob(plist,nrow=3,ncol=4)
print(p2)
dev.off()
print("finish")
