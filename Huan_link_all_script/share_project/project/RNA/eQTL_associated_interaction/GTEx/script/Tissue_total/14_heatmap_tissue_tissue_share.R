library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tidyverse)
library(pheatmap)
library(reshape2)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------
setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_total/")

org <- read.table("10_tissue_pair_overlap_index.txt.gz",header = T,sep = "\t") %>% as.data.frame()
org2 <-org[,c(1,3,2)]
colnames(org2)[2] ="tissue1"
colnames(org2)[3] ="tissue2"
org_final <-unique(bind_rows(org,org2))
org_final =org_final[,c(2,3,1)]
m <- dcast(data=org_final,tissue1~tissue2)
rownames(m) <-m[,1]
n <-m[,-1]
save(n,file="tissue_tissue_int.Rdata")
names_tissue <-colnames(n)
# names_tissue <- gsub("-", "_", names_tissue)

# annotation_row = data.frame(tissue_name = factor(colnames(n)))
# rownames(annotation_row) <- annotation_row$tissue_name
# annotation_col = data.frame(tissue_name = factor(colnames(n)))
# rownames(annotation_col) <- annotation_col$tissue_name
# annotation_row$tissue_name <- gsub("-", "_", annotation_row$tissue_name)
# annotation_col$tissue_name <- gsub("-", "_", annotation_col$tissue_name)

ann_colors = list(
  tissue = c(Adipose_Subcutaneous ="#f98404", Adipose_Visceral_Omentum = "#f9b208", Adrenal_Gland = "#61b15a", Artery_Aorta="#ee9595", Artery_Coronary="#f4a9a8",  Artery_Tibial="#cf0000", Brain_Amygdala= "#f7ea00",Brain_Anterior_cingulate_cortex_BA24="#f7ea00", Brain_Caudate_basal_ganglia="#f7ea00", Brain_Cerebellar_Hemisphere="#f7ea00", Brain_Cerebellum="#f7ea00", Brain_Cortex="#f7ea00", Brain_Frontal_Cortex_BA9="#f7ea00", Brain_Hippocampus="#f7ea00", Brain_Hypothalamus="#f7ea00", Brain_Nucleus_accumbens_basal_ganglia="#f7ea00", Brain_Putamen_basal_ganglia="#f7ea00", Brain_Spinal_cord_cervical_c_1 ="#f7ea00", Brain_Substantia_nigra="#f7ea00", Breast_Mammary_Tissue="#28abb9", Cells_Cultured_fibroblasts="#a6d6d6", Cells_EBV_transformed_lymphocytes="#8f4068", Colon_Sigmoid="#ffab73", Colon_Transverse="#ffb26b", Esophagus_Gastroesophageal_Junction="#c19277", Esophagus_Mucosa="#532e1c", Esophagus_Muscularis="#f1ae89", Heart_Atrial_Appendage="#93329e", Heart_Left_Ventricle="#440a67", Kidney_Cortex="#ccffbd", Liver="#beca5c", Lung="#c0e218", Minor_Salivary_Gland="#14b1ab", Muscle_Skeletal="#d789d7", Nerve_Tibial="#f0a500", Ovary="#ffc1fa", Pancreas="#87431d", Pituitary="#96bb7c", Prostate="#e8eae6", Skin_Not_Sun_Exposed_Suprapubic="#23049d", Skin_Sun_Exposed_Lower_leg="#845ec2", Small_Intestine_Terminal_Ileum="#1cc5dc",Spleen="#52734d", Stomach="#eac8af", Testis="#9088d4", Thyroid="#064420", Uterus="#f09ae9", Vagina="#fa26a0", Whole_Blood="#de8971")   
)

# pdf("sssss.pdf")
# pheatmap(n)


# pheatmap(n, annotation_col = annotation_col, annotation_row = annotation_row,show_rownames = F,annotation_legend = F,cluster_cols = F,cluster_rows = F)
# pheatmap(n, annotation_col = annotation_col, annotation_row = annotation_row)
#     # annotation_colors = ann_colors,
#     # cluster_cols = F,
#     # cluster_rows = F,
#     # # annotation_legend = T,
#     # show_rownames = F
#     )


annotation_row = data.frame(tissue = factor(colnames(n)))
rownames(annotation_row) <- annotation_row$tissue
annotation_col = data.frame(tissue = factor(colnames(n)))
rownames(annotation_col) <- annotation_col$tissue
annotation_row$tissue <- gsub("-", "_", annotation_row$tissue)
annotation_col$tissue <- gsub("-", "_", annotation_col$tissue)
# pheatmap(n,cluster_col = FALSE, cluster_rows = FALSE,border=FALSE, annotation_col = annotation_col,
#          annotation_row = annotation_row,annotation_colors = ann_colors,
#          cellwidth = 10, cellheight = 6, main = "NAR")
# #------------
# pdf("sssss.pdf",width = 10,height=10)

pheatmap(n,cluster_col = FALSE, cluster_rows = FALSE,border=FALSE, annotation_col = annotation_col,
         annotation_row = annotation_row,annotation_colors = ann_colors,show_colnames = F,fontsize_row=5.5,
         annotation_legend = F,annotation_names_row = F, annotation_names_col = F,
         cellwidth = 6, cellheight = 6, 
         color = colorRampPalette(c("#3edbf0","white",  "#ff4646"))(50),
         
         file='14_heatmap_tissue_tissue_share_hotspot_before_color.pdf')

pheatmap(n,cluster_col = FALSE, cluster_rows = FALSE,border=FALSE, annotation_col = annotation_col,
         annotation_row = annotation_row,annotation_colors = ann_colors,show_colnames = F,fontsize_row=5.5,
         annotation_legend = F,annotation_names_row = F, annotation_names_col = F,
         cellwidth = 6, cellheight = 6, 
         color = colorRampPalette(c("#a2b6df","#0c3483"))(50),
         
         file='14_heatmap_tissue_tissue_share_hotspot_before.pdf')



# pheatmap(n,cluster_col = T, cluster_rows = T,border=FALSE, annotation_col = annotation_col,
#          annotation_row = annotation_row,annotation_colors = ann_colors,show_colnames = F,fontsize_row=5.5,
#          annotation_legend = F,annotation_names_row = F, annotation_names_col = F,
#          cellwidth = 6, cellheight = 6, 
#          color = colorRampPalette(c("#a2b6df","#0c3483"))(50),
         
#          file='14_heatmap_tissue_tissue_share_hotspot.pdf')


pheatmap(n,cluster_col = T, cluster_rows = T,border=FALSE, annotation_col = annotation_col,
         annotation_row = annotation_row,annotation_colors = ann_colors,show_colnames = F,fontsize_row=5.5,
         annotation_legend = F,annotation_names_row = F, annotation_names_col = F,
         cellwidth = 6, cellheight = 6, 
         color = colorRampPalette(c("#cceeff","#33bbff"))(50),
         
         file='14_heatmap_tissue_tissue_share_hotspot_cluster.pdf')

pheatmap(n,cluster_col = F, cluster_rows = F,border=FALSE, annotation_col = annotation_col,
         annotation_row = annotation_row,annotation_colors = ann_colors,show_colnames = F,fontsize_row=5.5,
         annotation_legend = F,annotation_names_row = F, annotation_names_col = F,
         cellwidth = 6, cellheight = 6, 
         color = colorRampPalette(c("#cceeff","#33bbff"))(50),
         
         file='14_heatmap_tissue_tissue_share_hotspot.pdf')
        #  color = colorRampPalette(c("#3edbf0","white",  "#ff4646"))(50),

        #  color = colorRampPalette(c("#3edbf0", "white", "#ffdf6b","firebrick3"))(50),
# dev.off()


#-------------





# annotation_colors = ann_colors




# unique_tissue <-unique(org$tissue2)

# aa <-org%>%dplyr::select(-Index)

# A  <- matrix(NA,nc = length(unique_tissue), nr = length(unique_tissue))

# rownames(A) <- unique_tissue; colnames(A) <- unique_tissue



# org <- unique(sm_ja_index)
# for(i in unique(sm_ja_index$tissue2)){
#     if(i %in% unique(sm_ja_index$tissue1)){
#         print(1)
#     }else{
#         print(i)
#     }
# }


# all_data<-inner_join(all_hotspot,tissue_specific_hotspot,by="Tissue")
# diff <-all_data%>%dplyr::select(-c(2,3))
# colnames(diff)[2] <-"Number"

# tissue_specific_hotspot$Class <- "Specific"
# diff$Class <- "All"
# final <-bind_rows(tissue_specific_hotspot,diff)
# final$Tissue<-factor(final$Tissue,levels= sort(diff$Tissue,decreasing = T))
# pdf("tissue_specific_and_all_tissue.pdf")
# p1 <-ggplot(final,aes(x=Tissue,y=Number,fill=Class))+geom_bar(stat='identity',width=0.5 )+scale_fill_manual(values=c("#f29191","#1eae98"))+p_theme+xlab("Tissue") +ylab("Number of hotspots") +coord_flip() 
# print(p1)   
# dev.off()





# library(pheatmap)
# # a<-read.table("COO.txt",sep="\t",header=T)
# # a<-a[order(a$Predicted.classification),]
# # table(a$Predicted.classification)
# # annotation_col = data.frame(Cell_Type=(c('ABC','GCB')[a$Predicted.classification]))
# # rownames(annotation_col) = paste("", 1:441, sep = "")
# # pheatmap(b,treeheight_row=100,treeheight_col=20,cluster_cols=TRUE,color=colorRampPalette(c("green","black","red"))(1000),
# #  border_color=NA,fontsize=10,fontsize_row=8,fontsize_col=16,file='111.jpg',annotation_col = annotation_col, scale = "row")
# annotation_col = data.frame(CellType = factor(c(rep("exp_D3",3), rep("exp_H", 3)))
# rownames(annotation_col) = paste("Test", 1:10, sep = "")
# head(annotation_col)
# heatmap <- pheatmap(a,scale = "row",cluster_col = FALSE, border=FALSE, treeheight_row = 30,cellwidth = 30, cellheight = 12, main = "Glycolysis PMID:24140020")