library(ggplot2)
library(Rcpp)
# library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gcookbook)
library(gridExtra)
library(Hmisc)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/ALL_eQTL/cis_trans/interval_18/annotation_out/factors_score/")


types= c("cis_1MB","cis_10MB","trans_1MB","trans_10MB")
figure_list <-list ()
box_list <-list()
i=1
for(type in types){
  org_hotspot <-read.table(paste("06_annotation_merge_factors_score_hotspot_",type,".txt.gz",sep=""),header = T,sep = "\t") %>% as.data.frame()
  org_non_hotspot <-read.table(paste("06_annotation_merge_factors_score_non_hotspot_",type,".txt.gz",sep=""),header = T,sep = "\t") %>% as.data.frame()
  org_hotspot$Class <-"YES"
  org_non_hotspot$Class <-"NO"
  all<-bind_rows(org_non_hotspot,org_hotspot)
  top_title <- str_replace(type,"_"," ")
  top_title <-capitalize(top_title)
  figure_list[[i]] <- ggplot(all,aes(x=Class,y=score,fill=Class))+geom_violin(width=0.8)+geom_boxplot(width=0.08)+scale_y_continuous(name = "Score") + scale_x_discrete(name = "Class")+
  ggtitle(top_title)+theme(plot.title = element_text(hjust = 0.5))
  #--------
  box_list[[i]] <- ggplot(all,aes(x=Class,y=score,fill=Class))+geom_boxplot(width=0.8)+scale_y_continuous(name = "Score") + scale_x_discrete(name = "Class")+
  ggtitle(top_title)+theme(plot.title = element_text(hjust = 0.5))

  print(i)
  i=i+1
}
top_title= "Interval 18"
pdf("boxplot_violin_plot_cutoff_7.3_interval_18_cis_trans.pdf",height = 8,width = 8) 
p1<-marrangeGrob(figure_list,nrow=2,ncol=2,top = top_title)  
print(p1)
dev.off()
pdf("boxplot_plot_cutoff_7.3_interval_18_cis_trans.pdf",height = 8,width = 8) 
p1<-marrangeGrob(box_list,nrow=2,ncol=2,top = top_title)  
print(p1)
dev.off()
