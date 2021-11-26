library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")


p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

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
#--------------------------------------------------------cutoff=7


org<-read.table("../output/count_number_by_emplambda_in_different_interval_cutoff_7_all_eQTL.txt",header = T,sep = "\t") %>% as.data.frame()


interval1<-seq(from=2000, to=8000, by=2000)
interval2<-seq(from=10000, to=60000, by=10000)
interval=c(interval1,interval2) #两个interval合并

# interval<-seq(from=10000, to=60000, by=10000)

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
        theme(plot.title = element_text(hjust = 0.5))
        i=i+1
    # print(p1)
}

pdf("bar_plot_cutoff7_snp_number_in_different_interval_all_eQTL.pdf",height = 40,width = 56) 
marrangeGrob(figure_list,nrow=4,ncol=3,top = "cutoff = 7" )  #指定list的具体位置，也可以直接用整个list组装，比如：figure_list
dev.off()



# #-----------------------------
# for (j in interval){
#     org1<-filter(org, interval==j)
#     title_name<-paste("Bar plot of emplambda int",j, "bp",sep = " ")
#     p1<-ggplot(org1, aes(x = emplambda, y =number ))+ 
#         geom_bar(stat = "identity")+
#         xlab("Emplambda")+
#         p.theme+
#         ylab("Number of snp") +
#         ggtitle(title_name)+
#         theme(plot.title = element_text(hjust = 0.5))
#     print(p1)
# }




# #-------------------------



#     figure_name <-paste("Bar_plot_of_emplambda_int",j,".png", sep = "")
#     title_name<-paste("Bar plot of emplambda int",j, "bp",sep = " ")
#     png(figure_name,height = 800,width = 1000) 
#     p1<-ggplot(date, aes(x = pos, y =emplambda ))+ 
#         geom_bar(stat = "identity")+
#         xlab("Emplambda")+
#         p.theme+
#         ylab("pos")+
#         ggtitle(title_name)+
#         theme(plot.title = element_text(hjust = 0.5))
#     print(p1)
#     dev.off()
# }


# #------------------------------------------------------