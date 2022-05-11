library(tidyverse)
library(stringr)
library(GenomicFeatures)
library(readr)
library(mygene)
library(R.utils)
setwd("/home/huanhuan/project/eQTL_Catalogue/script/prediction/Edgeprediction/script/chr1/output/")

org<-read.table("02_chr1_hotspot.csv.gz",header = T,sep = ",") %>% as.data.frame()

set.seed(123)
valid <- org[sample(nrow(org),round(nrow(org)*0.1)),]
# train <-setdiff(org,valid)
write.table(valid,"./test/021_chr1_hotsopt_valid.csv",quote = F,sep = ",",row.names = F)
# write.table(train,"./train/021_hotsopt_train.csv",quote = F,sep = ",",row.names = F)
# gzip("./train/021_hotsopt_train.csv")
gzip("./test/021_chr1_hotsopt_valid.csv")