library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

cutoffs <-seq(0.1,0.3,0.01)

# nums<-c(1:10)
tissue <-"Whole_Blood"


rs<-data.frame()

figure_list <-list ()
i=1
for(cutoff in cutoffs){
    # file_name<-paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/",tissue ,"_segment_hotspot_cutoff_",cutoff,".bed.gz")
    file_name<-paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_",cutoff,".bed.gz")
    a<-file.info(file_name)$size #the file size
    if (a > 20){ #the null compressed files is 20
        org<-read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
        colnames(org)[1] <-"CHR"
        colnames(org)[2] <-"start"
        colnames(org)[3] <-"end"
        org$hotspot_length <-org$end - org$start
        org_filter <-filter(org,hotspot_length >5)
        bbb <-table(org$hotspot_length)%>%data.frame
        colnames(bbb)[1]<-"Length"
        bbb$cutoff <-cutoff
        rs<-bind_rows(rs,bbb)

    }
}

output_file <-paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/",tissue,"_length_of_hotspot.txt")
write.table(rs,output_file,row.names = F, col.names = T,quote =F,sep="\t")


