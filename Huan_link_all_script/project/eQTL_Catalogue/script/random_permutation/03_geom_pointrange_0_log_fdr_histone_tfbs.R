library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(data.table)
library(Seurat)

setwd("/home/huanhuan/project/eQTL_Catalogue/output/enrichment/figure/0/")
org<-read.table("compare_0_jaacard_index_fisher_test_histone_tfbs_two_side.txt",header = T,sep = "\t") %>% as.data.frame()
org$marker<-str_replace(org$marker,"CHROMATIN_Accessibility","CA")
org <-org%>%arrange(OR) #order to plot
org$FDR <-p.adjust(org$p_value,method="fdr")
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))
credplot.gg_main <- function(d){
  # org$marker<-str_replace(org$marker,"_"," ")
p1<-ggplot(data=org, aes(x=marker,y=log(OR), ymin=log(conf_int_bottom) , ymax=log(conf_int_up)))
  p1<-p1+geom_pointrange(size = 0.001)+ scale_y_continuous(limits= c(-0.3, 1.35), breaks= seq(0,1.2,0.4))
  p1<-p1+coord_flip()# +  # flip coordinates (puts labels on y axis)
  p1<-p1+ylab("Odds ratio (log scale)")+xlab("Markers") 
  p1 <-p1+p_theme+theme(axis.text.x = element_text(color = "black"),axis.text.y = element_text(color = "black"))
  p1 <-p1+geom_text(aes(x=marker,y=1.3,label = significant)) 
}

# pdf("fisher_test_0_jaacard_index_forest_plot.pdf",height = 3.2,width = 4)
pdf("fisher_test_0_jaacard_index_forest_plot_log_histone_tfbs.pdf",height = 3.5,width = 3.8) #
org$marker <- factor(org$marker, levels=org$marker) #marker ranking
p<-credplot.gg_main(org)
p
dev.off()

#------------------------------------

#----------
p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), 
                                                panel.grid.minor = element_blank(), 
                                                panel.background = element_rect(color="black",size=1.2), 
                                                axis.title.x = element_text(size =12),
                                                axis.text.x = element_text(color="black",size=9),
                                                axis.text.y = element_text(color="black",size = 12),
                                                axis.ticks.y=element_blank())

p1<-ggplot(data=org, aes(x=log(OR),y=marker))+
  geom_vline(xintercept=0,colour="#12497f",linetype=2,size=0.6)+
  geom_errorbarh(aes(xmax =log(conf_int_up), xmin = log(conf_int_bottom)), height = 0.05)+
  geom_point(size=1,fill="#c00024",colour="#c00024",)+ #,shape=22
  theme_bw()+
  p_theme +
  labs(x="log(Odds ratio)",y="")+
  scale_x_continuous(limits= c(-0.3, 1.3), breaks= seq(0,1.2,0.4))

p2 <-ggplot(data=org, aes(x=-log10(FDR),y=marker))+
    geom_bar(stat = "identity", width = 0.3,fill="#12497f") +
    geom_vline(xintercept=-log10(0.05),colour="#c00024",linetype=2,size=0.6)+
    theme_bw()+
    p_theme+
    labs(x="-log10(FDR)",y="")+
    scale_x_continuous(limits= c(0, 72.2), breaks= seq(0,75,15))+
    theme(axis.text.y=element_blank(),
    axis.title.y = element_blank()) 
    # +margin(2, 2, 2, 2, "cm")

library(patchwork)
p3<- (p1+theme(plot.margin=unit(c(0,0,0,0),"cm")))+(p2+theme(plot.margin=unit(c(0,0,0,0),"cm"))) + plot_layout(nrow = 1, width = c(100,30))

# p3 <-p1+plot_spacer()+ p2 +plot_layout(nrow = 1, width = c(100,0,40))
ggsave("fisher_test_0_jaacard_index_forest_plot_log_histone_tfbs_forest.png",p3,dpi=300,width=6.5,height=5)
