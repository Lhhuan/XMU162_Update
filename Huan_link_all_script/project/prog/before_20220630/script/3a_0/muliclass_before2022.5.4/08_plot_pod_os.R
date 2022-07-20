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

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}

dat$pod_total_r <-  dat$pod_total
dat$pod_total_r[grep("1|2",dat$pod_total_r)] <- "1_2"
fit <- survfit(Surv(pfs_month_new, pro_status) ~pod_total, data=dat)
# pdf("./figure/10_3_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                # palette=c("#21618C","#B71C1C"),
                legend.labs=c("0","1","2","3","4"), #标签
                pval = TRUE,
                legend.title="POD",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/08_pod01234_pfs_ori.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~pod_total, data=dat)
# pdf("./figure/08_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                # palette=c("#21618C","#B71C1C"),
                legend.labs=c("0","1","2","3","4"), #标签
                pval = TRUE,
                legend.title="POD",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/08_pod01234_os_ori.png",p1$plot,dpi=300,height=5.2,width=5)
#------------合并1,2

fit <- survfit(Surv(pfs_month_new, pro_status) ~pod_total_r, data=dat)
# pdf("./figure/08_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                # palette=c("#21618C","#B71C1C"),
                legend.labs=c("0","1_2","3","4"), #标签
                pval = TRUE,
                legend.title="POD",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/08_pod_merge12_pfs_ori.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~pod_total_r, data=dat)
# pdf("./figure/08_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                # palette=c("#21618C","#B71C1C"),
                legend.labs=c("0","1_2","3","4"), #标签
                pval = TRUE,
                legend.title="POD",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/08_pod_merge12_os_ori.png",p1$plot,dpi=300,height=5.2,width=5)



dat$pod_total_merge34 <-  dat$pod_total
dat$pod_total_merge34[grep("3|4",dat$pod_total_merge34)] <- "3_4"
# dat$pod_total_merge34[grep("1|2",dat$pod_total_merge34)] <- "1_2"


fit <- survfit(Surv(pfs_month_new, pro_status) ~pod_total_merge34, data=dat)
# pdf("./figure/08_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                # palette=c("#21618C","#B71C1C"),
                legend.labs=c("0","1","2","3_4"), #标签
                pval = TRUE,
                legend.title="POD",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/08_pod_merge34_pfs_ori.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~pod_total_merge34, data=dat)
# pdf("./figure/08_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                # palette=c("#21618C","#B71C1C"),
                legend.labs=c("0","1","2","3_4"), #标签
                pval = TRUE,
                legend.title="POD",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/08_pod_merge34_os_ori.png",p1$plot,dpi=300,height=5.2,width=5)


dat$pod_total_merge1234 <-  dat$pod_total
dat$pod_total_merge1234[grep("3|4",dat$pod_total_merge1234)] <- "3_4"
dat$pod_total_merge1234[grep("1|2",dat$pod_total_merge1234)] <- "1_2"
dat$pod1234 <- dat$pod_total_merge1234

fit <- survfit(Surv(pfs_month_new, pro_status) ~pod1234, data=dat)
# pdf("./figure/08_pod24_pfs_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                # palette=c("#21618C","#B71C1C"),
                # legend.labs=c("0","1_2","3_4"), #标签
                pval = TRUE,
                legend.title="POD",
                title="PFS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/08_pod_merge1234_pfs_ori.png",p1$plot,dpi=300,height=5.2,width=5)

fit <- survfit(Surv(os_month_new, dead) ~pod1234, data=dat)
# pdf("./figure/08_pod24_Os_cutoff1_adjust3.pdf",height=5.2,width=5)
p1 <- ggsurvplot(fit,
                # palette=c("#21618C","#B71C1C"),
                # legend.labs=c("0","1_2","3_4"), #标签
                pval = TRUE,
                legend.title="POD",
                title="OS",
                xlab = " Time (Months)",
                xlim=c(0,120),break.time.by=20,
                ggtheme = custom_theme()
                )
ggsave("./figure/08_pod_merge1234_os_ori.png",p1$plot,dpi=300,height=5.2,width=5)


