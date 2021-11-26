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

cutoffs <-seq(0.1,0.3,0.01)

# nums<-c(1:10)
tissue <-"Lung"




figure_list <-list ()
i=1
for(cutoff in cutoffs){
    file_name<-paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/hotspot_cis_eQTL/interval_18/",tissue ,"_segment_hotspot_cutoff_",cutoff,".bed.gz")
    a<-file.info(file_name)$size #the file size
    if (a > 20){ #the null compressed files is 20
        org<-read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
        colnames(org)[1] <-"CHR"
        colnames(org)[2] <-"start"
        colnames(org)[3] <-"end"
        org$hotspot_length <-org$end - org$start
    #     org_filter <-filter(org,hotspot_length >5 &&hotspot_length <5000)
    #     bbb <-table(org$hotspot_length)%>%data.frame
        title_name <-paste0("Cutoff ",cutoff)

    # #------------
    # org$hotspot_length[which(org$hotspot_length >5000)]=5001
    # org$hotspot_length[which(org$hotspot_length <6)]=5
    # bbb <-table(org$hotspot_length)%>%data.frame

        # figure_list[[i]] <-ggplot(org_filter, aes(x =hotspot_length)) +geom_histogram(position="identity")+xlab("Length of segment") + ylab("Count")+p_theme+ggtitle(title_name)+
        # figure_list[[i]] <-ggplot() +geom_bar(bbb, aes(x =Var1,y=log10(Freq)),position="identity")+xlab("Length of segment") + ylab("Log(Count)")+p_theme+ggtitle(title_name)+
        # figure_list[[i]] <-ggplot(bbb) +geom_bar( mapping=aes(x =Var1,y=log10(Freq)),position="identity")+xlab("Length of segment") + ylab("Log(Count)")+p_theme+ggtitle(title_name)+
        #  ggplot(bbb, aes(x =Var1,y=log10(Freq))) +geom_bar(stat="identity")+xlab("Length of segment") + ylab("Log(Count)")+p_theme+ggtitle(title_name)+
         figure_list[[i]] <-ggplot(org, mapping=aes(x =hotspot_length)) +geom_bar()+xlab("Length of segment") + ylab("Count")+p_theme+ggtitle(title_name)+
         theme(plot.title = element_text(hjust = 0.5))+
        scale_x_continuous(limits=c(5, 5000),breaks = c(1000,2000,3000,4000,5000))
        # ggplot(org,aes(y=hotspot_length))+geom_boxplot()
        i = i+1
    }
}

# output_file <-paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/",tissue,"/Cis_eQTL/",tissue,"_length_of_hotspot.txt")
# write.table(rs,output_file1,row.names = F, col.names = T,quote =F,sep="\t")


top_title <-paste0("Length of Hotspot ",tissue)


setwd(paste0("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/figure/",tissue))
pdf(paste0("histogram_plot_length_of_hotspot_cutoff_7.3_cis_eQTL_in_",tissue,"_re2.pdf"))
p1<-marrangeGrob(figure_list,nrow=3,ncol=3,top = top_title)  
print(p1)
dev.off()
print("AAA")


