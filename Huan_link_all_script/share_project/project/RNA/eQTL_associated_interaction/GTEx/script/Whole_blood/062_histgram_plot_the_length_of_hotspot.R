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
j=18

library(Hmisc)

cutoffs <-seq(0.05,0.65,0.05)
cutoffs <-c(0.01,cutoffs)
# nums<-c(1:10)
tissue <-"Whole_Blood"




figure_list <-list ()
i=1
for(cutoff in cutoffs){
    file_name<-paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_",cutoff,".bed.gz")
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

str_tissue <-str_replace(tissue,"_"," ")
top_title <-paste0("Length of Hotspot ",str_tissue)


setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/figure/",tissue))
pdf(paste0("histogram_plot_length_of_hotspot_cutoff_7.3_cis_eQTL_in_",tissue,".pdf"))
p1<-marrangeGrob(figure_list,nrow=3,ncol=3,top = top_title)  
print(p1)
dev.off()
print("AAA")


