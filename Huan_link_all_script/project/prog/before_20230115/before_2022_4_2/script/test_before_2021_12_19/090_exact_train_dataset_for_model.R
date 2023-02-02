library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library("survival")
library("survminer")

setwd("/home/huanhuan/project/prog/output/")
load("07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A.Rdata")
pod_total0=which(dat$new_pod_total==0)
set.seed(112231) #/TOP1
# set.seed(1124)
# set.seed(178)
# set.seed(188)
# set.seed(2)
test_number0 <-sample(x=pod_total0, round(length(pod_total0)*1/3),replace = F)
pod_total1=which(dat$new_pod_total!=0)
set.seed(112231)
# set.seed(1124)
# set.seed(178)
# set.seed(188)
# set.seed(2)
test_number1 <-sample(x=pod_total1, round(length(pod_total1)*1/3),replace = F)

test_set_number = c(test_number0,test_number1)
train_set_number =setdiff(1:nrow(dat),test_set_number)

test=dat[test_set_number,]
train=dat[train_set_number,]
save(test_set_number,train_set_number,test,train,file="09_test_train_dataset.Rdata")
write.table(train,"09_train_dataset.txt",quote=F,sep="\t",row.names=FALSE)
