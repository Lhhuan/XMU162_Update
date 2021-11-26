library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)


QTLs = c("eQTL","caQTL","edQTL","hQTL","mQTL","pQTL","reQTL","sQTL","QTL")
n<-length(QTLs)

# xQTL = "eQTL"


# cutoff =24.25
cutoff =7.3
intervals =1000
j=intervals 

fun=function(i){
    xQTL =QTLs[i]
    print(i)
    print(xQTL)
    current_directory <- paste("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_",xQTL,"/", sep = "")
    setwd(current_directory)
    input_file <-paste("../../output/ALL_", xQTL, "/NHPoisson_emplambda_interval_",j,"_cutoff_",cutoff,"_all_",xQTL,"_cis_trans_distance.txt.gz", sep = "")
    org<-fread(input_file,header = T,sep = "\t") %>% as.data.frame()

    # figure_name <-paste("Point_plot_emplambda_distance_cutoff_",cutoff,"_interval_",j,"_all_",xQTL,".pdf", sep = "")
    # pdf(figure_name,height = 5,width = 10)

    figure_name <-paste("Point_plot_emplambda_distance_cutoff_",cutoff,"_interval_",j,"_all_",xQTL,".png", sep = "")
    png(figure_name,height = 1000,width = 3000)

    # p1<-ggplot(org, aes(x=SYMBOL, y=number, group=1)) + geom_line(linetype="solid", color ="#4a9ff5") +geom_point(size=0.6, color ="#4a9ff5" )
    # p1<-ggplot(org, aes(x = reorder(SYMBOL,number), y=number, group=1)) +geom_point(size=0.6, color ="black" ) 
        # p1<-p1 +coord_flip() #x和y交换坐标轴
    p1<-ggplot(org, aes(x = reorder(distance,emplambda), y=emplambda, group=1)) +geom_point(size=0.1, color ="black" ) 

    p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.title.y = element_text(size = 25),axis.title.x = element_text(size = 25),axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=25),axis.text.x = element_text(size=25,angle=90)) 
    p2<-p2+xlab("Distance of snp and trait") + ylab("Emplambda")
    p2
    dev.off()



    #-------------------------------------------------density plot

    p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

    figure_name_density <-paste("density_plot_distance_cutoff_",cutoff,"_interval_",j,"_all_",xQTL,".pdf", sep = "")
    pdf(figure_name_density,height = 15,width = 30)

    title_name<-paste("Density of distance between snp and trait")
    p<-ggplot (org,aes(x=distance))+ 
        geom_density() + 
        xlab("distance")+
        p.theme+
        ylab("Density")+
        ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
    p
    dev.off()
}

library(parallel)
a=mclapply(1:n,function(x){fun(x)},mc.cores=9)

