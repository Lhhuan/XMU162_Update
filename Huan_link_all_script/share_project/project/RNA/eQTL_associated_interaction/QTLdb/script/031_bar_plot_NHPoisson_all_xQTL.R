library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

QTLs = c("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL","cerQTL","lncRNAQTL","metaQTL","miQTL","riboQTL")
# xQTL = "QTL"


p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
#-------------------------------------------------------

# cutoffs= c(7, 24.25)
cutoffs= c(7, 7.3)
for(xQTL in QTLs ){
    current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",xQTL,"/", sep = "")
    setwd(current_directory)
    for (cutoff in cutoffs){
        input_file <-paste("../../output/ALL_",xQTL,"/count_number_by_emplambda_in_different_interval_",cutoff,"_all_",xQTL,".txt", sep = "")
        org<-read.table(input_file,header = T,sep = "\t") %>% as.data.frame()

        interval<-seq(from=1000, to=5000, by=1000)

        figure_list <-list ()
        i=1
        for (j in interval){
            org1<-filter(org, interval==j)
            title_name<-paste("Bar plot of emplambda interval",j, "bp",sep = " ")
            
            figure_list[[i]] <-ggplot(org1, aes(x = emplambda, y =number ))+ 
                geom_bar(stat = "identity")+
                xlab("Emplambda")+
                p.theme+
                ylab("Number of snp") +
                ggtitle(title_name)+
                theme(plot.title = element_text(hjust = 0.5))+
                theme(axis.title.x = element_text(size = 15), axis.text.y = element_text(size = 15),axis.title.y = element_text(size = 15), axis.text.x = element_text(size = 15))+
                theme(plot.title = element_text(size = 15))
                i=i+1
            # print(p1)
        }
        figure_name <-paste("bar_plot_cutoff_",cutoff,"_snp_number_in_different_interval_all_",xQTL,".pdf", sep = "")
        top_title<-paste(xQTL,", cutoff = ",cutoff, sep = "")
        pdf(figure_name,height = 15,width = 30) 
        p1<-marrangeGrob(figure_list,nrow=2,ncol=3,top = top_title)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
        print(p1)
        dev.off()
    }
}

    # interval1<-seq(from=2000, to=8000, by=2000)
    # interval2<-seq(from=10000, to=100000, by=10000)
    # interval=c(interval1,interval2) #两个interval合并

# org<-read.table("../output/count_number_by_emplambda_in_different_interval_cutoff_9.txt",header = T,sep = "\t") %>% as.data.frame()

# interval<-seq(from=10000, to=60000, by=10000)

# figure_list <-list ()
# i=1
# for (j in interval){
#     org1<-filter(org, interval==j)
#     title_name<-paste("Bar plot of emplambda interval",j, "bp",sep = " ")
    
#     figure_list[[i]] <-ggplot(org1, aes(x = emplambda, y =number ))+ 
#         geom_bar(stat = "identity")+
#         xlab("Emplambda")+
#         p.theme+
#         ylab("Number of snp") +
#         ggtitle(title_name)+
#         theme(plot.title = element_text(hjust = 0.5))
#         i=i+1
#     # print(p1)
# }

# pdf("cutoff9_snp_number_in_different_interval.pdf",height = 20,width = 28) 
# marrangeGrob(figure_list,nrow=3,ncol=2,top = "cutoff = 9" )  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
# # marrangeGrob(figure_list[1:13],nrow=4,ncol=4)  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
# dev.off()