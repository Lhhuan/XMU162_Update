library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
library(cowplot)
library(ggpubr)
library(Hmisc)
library("survival")
library("survminer")

setwd("/home/huanhuan/project/prog/script/3a_0/output/")
load("07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A_0.Rdata")
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

#------------------------------

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

fit <- survfit(Surv(pfs_month_new, pro_status) ~pod_total_merge12_34, data=test)
# pdf("./figure/08_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C", "#38aa34","#B71C1C"),
                legend.labs=c("Low","Intermediate","High"), #标签
                pval = TRUE,
                legend.title="POD",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/09_pod_merge12_34_pfs_ori_test.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~pod_total_merge12_34, data=test)
# pdf("./figure/08_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                palette=c("#21618C", "#38aa34","#B71C1C"),
                legend.labs=c("Low","Intermediate","High"), #标签
                pval = TRUE,
                legend.title="POD",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/09_pod_merge12_34_os_ori_test.png",p1$plot,dpi=300,height=5.2,width=5)


#----------------
#------------------------------------------------------------------------------------------
# setwd("/home/huanhuan/project/prog/script/3a_0/output/")
# load("07_add_model_index_age_raw_pfs_os_filter_grade_mergrPod_3A_0.Rdata")
# pod_total0=which(dat$pod_total_merge12_34==0)
# set.seed(1121) #/TOP1
# test_number0 <-sample(x=pod_total0, round(length(pod_total0)*1/3),replace = F)

# pod_total1=which(dat$pod_total_merge12_34=="1_2")
# set.seed(1121)
# test_number1 <-sample(x=pod_total1, round(length(pod_total1)*1/3),replace = F)

# pod_total3=which(dat$pod_total_merge12_34=="3_4")
# set.seed(1121)
# test_number3 <-sample(x=pod_total3, round(length(pod_total3)*1/3),replace = F)

# test_set_number = c(test_number0,test_number1,test_number3)
# train_set_number =setdiff(1:nrow(dat),test_set_number)

# test=dat[test_set_number,]
# train=dat[train_set_number,]
# save(test_set_number,train_set_number,test,train,file="09_test_train_dataset.Rdata")
# write.table(train,"09_train_dataset.txt",quote=F,sep="\t",row.names=FALSE)