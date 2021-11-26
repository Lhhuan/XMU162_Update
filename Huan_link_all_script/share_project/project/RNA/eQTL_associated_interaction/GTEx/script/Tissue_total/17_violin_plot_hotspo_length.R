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



p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_total/")
read_f<-function(tissue=NULL){
    tissue1 <-gsub("Whole_Blood","whole_blood",tissue)
    file_name <-paste0("../../output/",tissue,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/",tissue1,"_segment_hotspot_cutoff_0.176.bed.gz")
    org <- read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
    print(tissue)
    return(org)
}
tissues <-c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta","Brain_Anterior_cingulate_cortex_BA24","Brain_Caudate_basal_ganglia","Brain_Cerebellum","Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Spinal_cord_cervical_c-1","Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes","Colon_Sigmoid","Colon_Transverse","Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis","Heart_Atrial_Appendage","Heart_Left_Ventricle","Kidney_Cortex","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg","Small_Intestine_Terminal_Ileum","Spleen","Stomach","Uterus","Prostate","Brain_Cerebellar_Hemisphere","Testis","Brain_Nucleus_accumbens_basal_ganglia","Minor_Salivary_Gland","Cells_Cultured_fibroblasts","Pituitary","Vagina","Thyroid","Artery_Tibial","Artery_Coronary","Brain_Hypothalamus","Nerve_Tibial","Brain_Putamen_basal_ganglia","Brain_Amygdala","Breast_Mammary_Tissue","Liver","Lung","Ovary","Pancreas","Whole_Blood")
tmp <-lapply(tissues,read_f)

rs <-do.call(rbind,tmp)

colnames(rs) <- c("chr","start","end")
rs$length_of_hotspot = rs$end-rs$start
unique_rs <-unique(rs)%>%as.data.frame()

pdf("17_violin_plot_hotspo_length.pdf",width=3.5, height=3.5)
p<-ggplot(unique_rs,aes(x=1, y=log10(length_of_hotspot)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1)+ theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("Log10(length of hotspot)")

print(p)
dev.off()

count_result <-data.frame(mean=mean(rs$length_of_hotspot),min=min(rs$length_of_hotspot),max=max(rs$length_of_hotspot),sd=sd(rs$length_of_hotspot))

write.table(count_result ,"17_share_hotspot_union_gene.txt",row.names = F, col.names = T,quote =F,sep="\t")