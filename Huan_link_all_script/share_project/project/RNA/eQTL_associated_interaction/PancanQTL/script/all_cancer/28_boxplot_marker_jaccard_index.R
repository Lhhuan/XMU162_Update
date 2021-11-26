library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

# cancers <-c("BRCA","KIRC","KIRP","LAML","LIHC","LUAD","LUSC","OV","PAAD","PRAD","STAD","THCA")
cancers <-c("BRCA","KIRC","KIRP","LAML","LIHC","LUAD","OV","PAAD","PRAD","STAD","THCA")

c_plot <-function(cancer =NULL){
# for(cancer in cancers){
    print(cancer)
    path = paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/",cancer,"/Cis_eQTL/anno/")
    setwd(path)
    print(path)
    cancer_file <- paste0(cancer,"_segment_hotspot_cutoff_0.176_marker_jaccard_index_label.txt.gz")
    org <- read.table(cancer_file,header = T,sep = "\t") %>% as.data.frame()

    tri_markers <-c("CHROMATIN_Accessibility","TFBS","H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
    unique_marker <-unique(org$marker)
    overlap_marker <- unique_marker[unique_marker%in%tri_markers]

    p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                    panel.background = element_blank(), axis.title.y = element_text(size = 8),
                                                    # axis.title.x = element_text(size = 10),
                                                    axis.line = element_line(colour = "black"))
    my_comparisons <- list(c("Specific","Others"))
    # my_comparisons <- list(c("Others","specific"))
    large <-function(state){
        aa <-filter(org,marker==state)
        title_name <-str_replace(state,"CHROMATIN_Accessibility", "CA")
        # title_name<-state
        p <-ggplot(aa,aes(x=Class,y=jaccard_index,fill = Class))+ geom_violin()+
            ggtitle(title_name)+
            theme(legend.position ="none",plot.title = element_text(hjust = 0.5),axis.title.x= element_blank())+
            ylab("Jaccard index")+
            p_theme+coord_cartesian(ylim = c(0, 1.1))+
            stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))+
            scale_x_discrete(limits=rev(levels(aa$Class)))
    }
    
    plist = lapply(overlap_marker,large)
    pdf("large_typical_mark_jaccard_index.pdf",width=8.3, height=6.5)
    pp <-CombinePlots(plist,ncol=4,nrow=3)
    print(pp)
    dev.off()

    plist = lapply(unique_marker,large)
    pdf("large_all_mark_jaccard_index.pdf",width=11, height=11)
    # print(CombinePlots(plist,ncol=5,nrow=5))
    ppp <-CombinePlots(plist,ncol=5,nrow=5)
    print(ppp)
    dev.off()

    less <-function(state){
    aa <-filter(org,marker==state)
    title_name <-str_replace(state,"CHROMATIN_Accessibility", "CA")
    # title_name<-state
    p <-ggplot(aa,aes(x=Class,y=jaccard_index,fill = Class))+ geom_violin()+
        ggtitle(title_name)+
        theme(legend.position ="none",plot.title = element_text(hjust = 0.5),axis.title.x= element_blank())+
        ylab("Jaccard index")+
        p_theme+coord_cartesian(ylim = c(0, 1.1))+
        stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "less"))+
        scale_x_discrete(limits=rev(levels(aa$Class)))
    }
    #----------
    plist = lapply(overlap_marker,less)
    pdf("less_typical_mark_jaccard_index.pdf",width=8.3, height=6.5)
    pp <-CombinePlots(plist,ncol=4,nrow=3)
    print(pp)
    dev.off()

    plist = lapply(unique_marker,less)
    pdf("less_all_mark_jaccard_index.pdf",width=11, height=11)
    # print(CombinePlots(plist,ncol=5,nrow=5))
    ppp <-CombinePlots(plist,ncol=5,nrow=5)
    print(ppp)
    dev.off()
    #----------
}

lapply(cancers,c_plot)


#------------------try log transform

cancer ="BRCA"

    print(cancer)
    path = paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/PancanQTL/output/",cancer,"/Cis_eQTL/anno/")
    setwd(path)
    print(path)
    cancer_file <- paste0(cancer,"_segment_hotspot_cutoff_0.176_marker_jaccard_index_label.txt.gz")
    org <- read.table(cancer_file,header = T,sep = "\t") %>% as.data.frame()

    tri_markers <-c("CHROMATIN_Accessibility","TFBS","H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3","H3K27me3","H3K9me3")
    unique_marker <-unique(org$marker)
    overlap_marker <- unique_marker[unique_marker%in%tri_markers]
    my_comparisons <- list(c("Specific","Others"))
    less <-function(state){
        aa <-filter(org,marker==state)
        title_name <-str_replace(state,"CHROMATIN_Accessibility", "CA")
        # title_name<-state
        p <-ggplot(aa,aes(x=Class,y=log10(jaccard_index+0.001),fill = Class))+ geom_violin()+
            ggtitle(title_name)+
            theme(legend.position ="none",plot.title = element_text(hjust = 0.5),axis.title.x= element_blank())+
            ylab("Log(Jaccard index)")+
            p_theme+coord_cartesian(ylim = c(0, 1.1))+
            stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "less"))+
            scale_x_discrete(limits=rev(levels(aa$Class)))
    }

    plist = lapply(overlap_marker,less)
    pdf("log_less_typical_mark_jaccard_index.pdf",width=8.3, height=6.5)
    pp <-CombinePlots(plist,ncol=4,nrow=3)
    print(pp)
    dev.off()



# hotspot<-read.table("hotspot_cutoff_0.176_histone_marker_jaccard_index.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# # random<-fread("0_0.176_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t")
# random<-read.table("0_0.176_jaccard_index_histone_marker_1000.txt.gz",header = T,sep = "\t") %>% as.data.frame()


# active_mark <-c("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3")
# repressed_mark <-c("H3K27me3","H3K9me3")

# hotspot$Random_number <- 1
# hotspot$Class <- "hotspot"
# random$Class <- "random"

# p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#                                                 panel.background = element_blank(), axis.title.y = element_text(size = 8),
#                                                 # axis.title.x = element_text(size = 10),
#                                                 axis.line = element_line(colour = "black"))

# rs<- bind_rows(hotspot,random)

# # figure_list <-list()
# my_comparisons <- list(c("hotspot","random"))

# large <-function(state){
#     aa <-filter(rs,Marker==state)
#     title_name<-state

#     p <-ggplot(aa,aes(x=Class,y=jaacard_index,fill = Class))+ geom_violin(outlier.colour = NA)+ theme(legend.position ="none") +ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5))+ylab("Jaccard index")+p_theme+coord_cartesian(ylim = c(0, 1.1))+
#     stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", method.args = list(alternative = "greater"))
# }
 
# plist = lapply(active_mark,large)

# pdf("./figure/0_0.176/compare_0_0.176_active_mark_jaccard_index.pdf",width=6.3, height=4.5)
# CombinePlots(plist,ncol=3,nrow=2)
# dev.off()

# plist = lapply(repressed_mark,large)

# pdf("./figure/0_0.176/compare_0_0.176_repressed_mark_jaccard_index.pdf",width=4, height=2.25)
# CombinePlots(plist,ncol=2)
# dev.off()
# print("aaa")


# #----
