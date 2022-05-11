library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(readr)
library(mygene)
library(R.utils)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/prediction/Edgeprediction/script/chr1/output/")

org<-read.table("chr1_huan_refine_temporal_validation_gs_nra.csv",header = T,sep = ",") %>% as.data.frame()

pred <-org[!is.na(org$auc),]

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

pdf("04_chr1_auc_distribution.pdf",width=3.5, height=3.5)
p<-ggplot(pred,aes(x=1, y=auc))+
    # geom_violin(fill="#a3d2ca",width=0.65)+
    geom_boxplot(fill = "#5eaaa8",outlier.color=NA)+ 
    p_theme+
    theme(legend.position ="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
    xlab("")+ylab("AUC")

print(p)
dev.off()






set.seed(123)
valid <- org[sample(nrow(org),round(nrow(org)*0.1)),]
# train <-setdiff(org,valid)
write.table(valid,"./test/021_chr1_hotsopt_valid.csv",quote = F,sep = ",",row.names = F)
# write.table(train,"./train/021_hotsopt_train.csv",quote = F,sep = ",",row.names = F)
# gzip("./train/021_hotsopt_train.csv")
gzip("./test/021_chr1_hotsopt_valid.csv")