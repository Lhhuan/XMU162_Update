library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(parallel)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------

# cutoffs= c(7, 24.25)

xQTL = "eQTL"
cutoffs= 7.3

intervals=c(6,7,8,9,12,15,18)

types= c("cis_1MB","cis_10MB","trans_1MB","trans_10MB")
groups = c("hotspot","non_hotspot")

library(Hmisc)
# tops=c(1500,2000)
top = 2500
bottom = 6
#-------------------------------------------
#-----coord_cartesian(xlim = c(-2, 2))
ProcessBedGz <- function(j = NULL){
    for(group in groups){
        figure_list <-list ()
        Density_figure_list<-list()
        i=1
        setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/")
        for (type in types){ 
            org<-read.table(paste0(group,"/","interval_",j,"_cutoff_7.3_",type,"_eQTL_segment_",group,"_length_of_segment.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
            colnames(org)[4]<-"count"
            new_type <-str_replace(type,"_"," ")
            new_type <-capitalize(new_type) #Capital letter
            org1 <-filter(org,count >=bottom)
            org2 <-filter(org1,count <=top)
            max_limit <- max(org2$count)
            min_limit <-min(org2$count)
            
            # title_name<-paste(new_type,": interval ",j,sep = "")
            title_name<- new_type                
                # print(p1)
            # Density_figure_list[[i]] <- ggplot(org2, aes(x =count)) + geom_density()+
                # xlab("Length of segment") +ylab("Density")+p_theme+
                # ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5))
                
            figure_list[[i]] <-ggplot(org2, aes(x =count)) +geom_histogram(position="identity")+
                xlab("Length of segment") + ylab("Count")+p_theme+ggtitle(title_name)+
                # coord_cartesian(xlim = c(min_limit , max_limit))+
                theme(plot.title = element_text(hjust = 0.5))
            i=i+1
        }

        setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/figure/ALL_eQTL/cis_trans/")
        top_title <-paste("Interval",j,sep=" ")

        # # pdf(paste0("density_plot_cutoff_7.3_interval_",j,"_cis_trans.pdf"),height = 7.5,width = 7) 
        # # pdf(paste0(group,"/Density_plot_length_of_segment_",group,"_cutoff_7.3_interval_",j,".pdf"),height = 7,width = 9) 
        # pdf(paste0(group,"/Density_plot_length_of_segment_",group,"_cutoff_7.3_interval_",j,"_bottom_",bottom,"_top_",top,".pdf")) 
        # p1<-marrangeGrob(Density_figure_list,nrow=2,ncol=2,top = top_title)  
        # print(p1)
        # dev.off()
        #------------

        # pdf(paste0(group,"/Histogram_plot_length_of_segment_",group,"_cutoff_7.3_interval_",j,".pdf"),height = 7,width = 9) 
        pdf(paste0(group,"/Histogram_plot_length_of_segment_",group,"_cutoff_7.3_interval_",j,"_bottom_",bottom,"_top_",top,".pdf")) 
        p1<-marrangeGrob(figure_list,nrow=2,ncol=2,top = top_title)  
        print(p1)
        dev.off()
        # print(j)
        print(group)

    }
}


mclapply(intervals, ProcessBedGz, mc.cores = 10)


#--------------
# for(j in intervals){
# ProcessBedGz <- function(j = NULL){
#     for(group in groups){
#         figure_list <-list ()
#         Density_figure_list<-list()
#         i=1
#         setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/output/ALL_eQTL/cis_trans/")
#         for (type in types){ 
#             org<-read.table(paste0(group,"/","interval_",j,"_cutoff_7.3_",type,"_eQTL_segment_",group,"_length_of_segment.bed.gz"),header = F,sep = "\t") %>% as.data.frame()
#             colnames(org)[4]<-"count"
#             new_type <-str_replace(type,"_","")
#             new_type <-capitalize(new_type) #Capital letter
#             # title_name<-paste(new_type,": interval ",j,sep = "")
#             title_name<- new_type

#             Density_figure_list[[i]] <- ggplot(org, aes(x =count)) + geom_density()+scale_y_log10()+
#                 xlab("length") +ylab("Density")+p_theme+
#                 ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5))

#             figure_list[[i]] <-ggplot(org, aes(x =count)) +geom_histogram(position="identity")+scale_y_log10()+
#                 xlab("length") + ylab("Log(Count)")+p_theme+ggtitle(title_name)+
#                 #coord_cartesian(xlim = c(0, 1))+
#                 theme(plot.title = element_text(hjust = 0.5))

#             i=i+1
#         }

#         setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLdb/figure/ALL_eQTL/cis_trans/")
#         top_title <-paste("Interval",j,sep=" ")

#         pdf(paste0(group,"/Density_plot_length_of_segment_",group,"_cutoff_7.3_interval_",j,".pdf"),height = 7.5,width = 7) 
#         p1<-marrangeGrob(Density_figure_list,nrow=2,ncol=2,top = top_title)  
#         print(p1)
#         # marrangeGrob(Density_figure_list,nrow=2,ncol=2,top = top_title)
#         dev.off()

#         pdf(paste0(group,"/Histogram_plot_length_of_segment_",group,"_cutoff_7.3_interval_",j,".pdf"),height = 7.5,width = 7) 
#         p1<-marrangeGrob(figure_list,nrow=2,ncol=2,top = top_title)  
#         print(p1)
#         # marrangeGrob(figure_list,nrow=2,ncol=2,top = top_title) 
#         # marrangeGrob(figure_list,nrow=2,ncol=2) 
#         dev.off()

#         # pdf(paste0("density_plot_cutoff_7.3_interval_",j,"_cis_trans.pdf"),height = 7.5,width = 7) 

#         print(j)
#         print(group)

#     }
# }


