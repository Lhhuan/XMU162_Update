library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(conflicted)
library(parallel)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------
library(Hmisc)

factors<-c("HISTONE_MARK_AND_VARIANT","Human_CHROMATIN_Accessibility","Human_FACTOR")

# for(factor in factors){
ProcessBedGz <-function(factor =NULL){
        setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/normal_cell/",factor))


        org<-read.table("merge_pos_info_narrow_peak.bed.gz",header = F,sep = "\t") %>% as.data.frame()
        org$length <-org$V3-org$V2
        title_name <- "Normal cell line"
        pdf(paste0("Distribution_of_length_of_",factor,".pdf")) 
        p1<-ggplot(org, aes(x =length)) +geom_histogram(position="identity")+
                xlab("Length") + ylab("Count")+p_theme+ggtitle(title_name)+
                theme(plot.title = element_text(hjust = 0.5))

        print(p1)
        dev.off()




        colnames(org)[1] <-"chr"
        colnames(org)[2] <-"start"
        colnames(org)[3] <-"end"
        output1<-paste0("06_",factor,"_length.txt")
        write.table(org,output1,row.names = F, col.names = T,quote =F,sep="\t")

        mean_length<- round(mean(org$length))
        median.length <- median(org$length)
        rs <-data.frame(mean_length =mean_length,median_length = median.length)
        output2<-paste0("06_",factor,"_length_statistics.txt")
        write.table(rs,output2,row.names = F, col.names = T,quote =F,sep="\t")
        print(factor)
}



mclapply(factors, ProcessBedGz, mc.cores = 5)
