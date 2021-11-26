library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------

# xQTL = "eQTL"
# cutoff= 7.3

tissues = c("Lung","Whole_Blood")
# intervals=c(6,7,8,9,12,15,18)
j=18

library(Hmisc)

cutoffs <-seq(0.05,0.65,0.05)
cutoffs <-c(0.01,cutoffs)
nums<-c(1:10)


for(tissue in tissues){
    setwd(paste0("/share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/",tissue,"/figure/"))
    for(number in nums){
        print(number)
        figure_list <-list ()
        i=1
        for(cutoff in cutoffs){
            file_name<-paste0("/share/data0/QTLbase/huan/GTEx/random_select_exclude_eQTL/",tissue,"/hotspot/",number,"/interval_18/",tissue,"_segment_hotspot_cutoff_",cutoff,".bed.gz")
            a<-file.info(file_name)$size #the file size
            if (a > 20){ #the null compressed files is 20
                org<-read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
                colnames(org)[1] <-"CHR"
                colnames(org)[2] <-"start"
                colnames(org)[3] <-"end"
                org$hotspot_length <-org$end - org$start - 1

                title_name <-paste0("Cutoff ",cutoff)
                figure_list[[i]] <-ggplot(org, aes(x =hotspot_length)) +geom_histogram(position="identity")+xlab("Length of segment") + ylab("Count")+p_theme+ggtitle(title_name)+
                theme(plot.title = element_text(hjust = 0.5))
                i = i+1
            }
        }
        top_title <-paste0("Length of Hotspot in Random select ",number," in ",tissue)
        # pdf(paste0("histogram_plot_cutoff_7.3_cis_eQTL_in_",tissue,".pdf"),compress = TRUE) 
        pdf(paste0("histogram_plot_length_of_hotspot_cutoff_7.3_cis_eQTL_in_",tissue,"_Random_select_",number,".pdf"))
        p1<-marrangeGrob(figure_list,nrow=2,ncol=3,top = top_title)  
        print(p1)
        dev.off()
    }
}
