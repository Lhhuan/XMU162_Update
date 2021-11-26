library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Whole_Blood/Cis_eQTL/enrichment/interval_18/ALL/")
hotspot<-read.table("hotspot_cutoff_0.176_histone_marker.txt",header = T,sep = "\t") %>% as.data.frame()
random<-read.table("original_random_histone_marker.txt.gz",header = T,sep = "\t") %>% as.data.frame()

hotspot$random_number <- 1
hotspot$Class <- "Hotspot"
random$Class <- "Random expectaion"
# hotspot$percentage <- hotspot$number/46369*100
# random$percentage <-random$number/46369*100
rs<- bind_rows(hotspot,random)

unique_Marker <- unique(rs$Marker)



p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 8),
                                                # axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))


#---------------
#-------------"1_TssA"
symnum.args <- list(cutpoints = c(0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

active_mark <-c("H3K27ac","H3K9ac","H3K36me3","H3K4me1","H3K4me3")
repressed_mark <-c("H3K27me3","H3K9me3")


large <-function(state){
    aa <-filter(rs,Marker==state)
    Hotspot <-filter(aa,Class=="Hotspot")
    Random_expectaion <-filter(aa,Class!="Hotspot")
    observed <-Hotspot$count
    used_count_p <- nrow(filter(Random_expectaion,count >observed))
    p_value = used_count_p/nrow(Random_expectaion)
    print(p_value)

    title_name<-state
    state <-str_replace(state,"/","_")
    #pdf(paste0("./figure/percentage/15_state_",state,".pdf"),width=4, height=4) 
    p <-ggplot(aa,aes(x=Class,y=count))+geom_boxplot(aes(fill=Class),width=0.3,outlier.colour = NA)+ 
    #scale_y_continuous(limits=c(0,5)) + 
    theme(legend.position ="none")+ggtitle(title_name) +theme(plot.title = element_text(hjust = 0.5),axis.title.x=element_blank(),axis.text.x = element_text(size = 7, color = "black"),axis.text.y = element_text(color = "black")) +ylab("Count of segments")+p_theme

    #add_pval(p,annotation = p_value)

    if(p_value <= 0.0001){
        add_pval(p,annotation = "****",pval_star = T)
    }else if(p_value <= 0.001){
        add_pval(p,annotation = "***",pval_star = T)
    }else if(p_value <= 0.01){
        add_pval(p,annotation = "**",pval_star = T)
    }else if(p_value <= 0.05){
        add_pval(p,annotation = "*",pval_star = T)
    }else{
        add_pval(p,annotation = "ns",pval_star = T)
    }
    
}
 
plist = lapply(active_mark,large)

pdf("./figure/count/compare_orginal_random_active_mark_count.pdf",width=8.3, height=4)
CombinePlots(plist,ncol=4,nrow=2)
dev.off()

plist = lapply(repressed_mark,large)

pdf("./figure/count/compare_orginal_random_repressed_mark_count.pdf",width=4, height=2)
CombinePlots(plist,ncol=2)
dev.off()

