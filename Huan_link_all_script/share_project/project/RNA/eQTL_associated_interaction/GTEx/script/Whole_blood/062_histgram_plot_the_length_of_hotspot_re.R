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

cutoffs <-0.176
# cutoffs <-c(0.01,cutoffs)
# nums<-c(1:10)
tissue <-"Whole_Blood"





for(cutoff in cutoffs){
    file_name<-paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/whole_blood_segment_hotspot_cutoff_",cutoff,".bed.gz")
    a<-file.info(file_name)$size #the file size
    if (a > 20){ #the null compressed files is 20
        org<-read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
        colnames(org)[1] <-"CHR"
        colnames(org)[2] <-"start"
        colnames(org)[3] <-"end"
        org$hotspot_length <-org$end - org$start
        #----------
        aa <-table(org$hotspot_length)%>%data.frame()

        #--------
        breaks1 = seq(0,200000,10000)
        title_name <-paste0("Cutoff ",cutoff)
        p <-ggplot(aa, aes(x =Var1,y=Freq)) +geom_bar(stat="identity")+xlab("Length of segment") + ylab("Count")+p_theme+ggtitle(title_name)+
        #scale_x_continuous(breaks=breaks1)+
        theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(angle=90)))

        setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/figure/",tissue))
        pdf(paste0("bar_plot_length_of_hotspot_cutoff_0.176_7.3_cis_eQTL_in_",tissue,".pdf"))
        print(p)
        dev.off()





        title_name <-paste0("Cutoff ",cutoff)
        p <-ggplot(org, aes(x =hotspot_length)) +geom_bar(position="identity")+xlab("Length of segment") + ylab("Count")+p_theme+ggtitle(title_name)+
        theme(plot.title = element_text(hjust = 0.5))
        setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/figure/",tissue))
        pdf(paste0("bar_plot_length_of_hotspot_0.176_cutoff_7.3_cis_eQTL_in_",tissue,".pdf"))
        print(p)
        dev.off()
        org_filter <-filter(org,hotspot_length<=5000)
        nrow(org_filter)
        p2 <-ggplot(org_filter, aes(x =hotspot_length)) +geom_bar(position="identity")+xlab("Length of segment") + ylab("Count")+p_theme+ggtitle(title_name)+theme(plot.title = element_text(hjust = 0.5))
                
        pdf(paste0("bar_plot_length_of_hotspot_cutoff_0.176_7.3_less_than5000_cis_eQTL_in_",tissue,".pdf"))
        print(p2)
        dev.off()


    }
}

str_tissue <-str_replace(tissue,"_"," ")
top_title <-paste0("Length of Hotspot ",str_tissue)


setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/figure/",tissue))
pdf(paste0("histogram_plot_length_of_hotspot_0.176_cutoff_7.3_cis_eQTL_in_",tissue,"_re.pdf"))
p1<-marrangeGrob(figure_list,nrow=3,ncol=3,top = top_title)  
print(p1)
dev.off()
print("AAA")


