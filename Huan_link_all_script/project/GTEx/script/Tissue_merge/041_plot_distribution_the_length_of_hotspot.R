library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)

p_theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

#---------------------------------------------------------
library(Hmisc)

org<-read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/Tissue_merge/Cis_eQTL/hotspot_cis_eQTL/interval_18/Tissue_merge_segment_hotspot_cutoff_0.176.bed.gz",header = F,sep = "\t") %>% as.data.frame()
colnames(org) <-c("CHR","start","end")
org$hotspot_length <-org$end - org$start
org$Log10_len <-log10(org$hotspot_length)
setwd("/home/huanhuan/project/GTEx/script/Tissue_merge/figure/")
pdf("041_boxplot_distribution_original_hotspot_length_log10.pdf",width=3.5, height=3.5)

p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1,outlier.color=NA)+ 
    theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y=element_text(size=5,color="black"))+xlab("")+ylab("log10(length of hotspot)")+
    coord_cartesian(ylim = c(0, 6.5))+scale_y_continuous(breaks=seq(0,6.5,0.25))

print(p)
dev.off()
#---------------------------
org$class <-NA
org$class[which(org$hotspot_length ==1)] <- "1"
org$class[which(org$hotspot_length <=5 & org$hotspot_length >1)] <- "2-5"
org$class[which(org$hotspot_length <=10 & org$hotspot_length >5)] <- "6-10"
org$class[which(org$hotspot_length <=50 & org$hotspot_length >10)] <- "11-50"
org$class[which(org$hotspot_length <=100 & org$hotspot_length >51)] <- "51-100"
org$class[which(org$hotspot_length <=300 & org$hotspot_length >101)] <- "101-300"
org$class[which(org$hotspot_length <=1000 & org$hotspot_length >301)] <- "301-1000"
org$class[which(org$hotspot_length <=10000 & org$hotspot_length >1001)] <- "1001-10000"
org$class[which(org$hotspot_length >10000)] <- ">10000"

org$class<-factor(org$class,levels=c("1","2-5","6-10","11-50","51-100","101-300","301-1000","1001-10000",">10000"))
l_c <- table(org$class)%>%data.frame()
colnames(l_c) <-c("Length","Frequency")

m_10000 <-filter(org,hotspot_length >10000)
m_10000 <-m_10000[,1:4]
write.table(m_10000,"length_more10000.bed",row.names = F, col.names = F,quote =F,sep="\t")

pdf("./041_barplot_hotspot_distribution.pdf",width = 5,height = 5)
p1 <-ggplot(data = l_c, mapping = aes(x =Length , y = Frequency)) + geom_bar(stat = 'identity', fill = "#5eaaa8", width=0.6)+
  p_theme+scale_y_continuous(expand=c(0,0))+
  theme(axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black",angle=60,hjust=1))
print(p1)
dev.off()

#coord_flip()+

# , width=0.5

#-------------------

quantile(org$hotspot_length)
o_25_75 <-filter(org,hotspot_length>24 &hotspot_length <314)%>%select("CHR","start","end")
write.table(o_25_75,"/home/huanhuan/project/GTEx/output/Tissue_merge/hotspot_quantile/org_hotspot_quantile25_75.bed",row.names = F, col.names = F,quote =F,sep="\t")





a <-as.data.frame(table(org$Log10_len))
a<-a[order(a$Freq,a$Var1),]
a$Var1 <- as.numeric(as.character(a$Var1))
a$Freq <- as.numeric(as.character(a$Freq))
llen_c <-a[1,1]
# b <-filter(a,Var1 <1)
# b<-b[order(b$Var1),]

# 3.0948203803548
# hist(org$Log10_len)
# dev.off()

orgfl <- filter(org,Log10_len <=llen_c)

pdf("041_boxplot_distribution_filter_L_hotspot_length_outlier.pdf",width=3.5, height=3.5)

p<-ggplot(orgfl,aes(x=1, y=hotspot_length))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1)+ 
    theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("length of hotspot")

print(p)
dev.off()

pdf("041_boxplot_distribution_filter_L_hotspot_length_log10_outlier.pdf",width=3.5, height=3.5)
p<-ggplot(orgfl,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1)+ 
    theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("log10(length of hotspot)")
print(p)
dev.off()



lo_25_75 <-filter(orgfl,hotspot_length>=23 &hotspot_length <=302)




hist(orgfl$hotspot_length)
dev.off()

c <- table(orgfl$hotspot_length)%>%as.data.frame()
c<-apply(c,2,as.character)
c<-apply(c,2,as.numeric)
# a$Var1 <- as.numeric(as.character(a$Var1))
# a$Freq <- as.numeric(as.character(a$Freq))
c<-c[order(c$Freq,c$Var1),]





pdf("041_boxplot_distribution_original_hotspot_length_log10_outlier.pdf",width=3.5, height=3.5)

p<-ggplot(org,aes(x=1, y=log10(hotspot_length)))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1)+ 
    theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("log10(length of hotspot)")

print(p)
dev.off()
#----------------
pdf("041_boxplot_distribution_original_hotspot_length_outlier.pdf",width=3.5, height=3.5)

p<-ggplot(org,aes(x=1, y=hotspot_length))+geom_violin(fill="#a3d2ca",width=0.65)+geom_boxplot(fill = "#5eaaa8",width=0.1)+ 
    theme(legend.position ="none")+p_theme+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+xlab("")+ylab("length of hotspot")

print(p)
dev.off()


stat <-data.frame(mean = mean(org$hotspot_length),min=min(org$hotspot_length),max=max(org$hotspot_length),sd=sd(org$hotspot_length))%>%t()%>%data.frame()
colnames(stat) <-"length"
stat <-add_column(stat, statistic=rownames(stat), .before = "length")
# stat <-mutate(stat, statistic=rownames(stat), .before = 1)
# # write.table(stat,"")
write.table(stat,"041_distribution_of_hotspot_length.txt",row.names = F, col.names = T,quote =F,sep="\t")

setwd("/home/huanhuan/project/GTEx/output/Tissue_merge/")
len1 <-filter(org,hotspot_length <2)
write.table(len1,"/home/huanhuan/project/GTEx/output/Tissue_merge/hotspot_len1.txt",row.names=F,col.names=F,quote=F,sep="\t")


#------------------------------------------------
org$chr= as.numeric(str_replace(org$CHR,"chr",""))

a <- filter(org,hotspot_length >=6 & Log10_len <3.0948203803548)

a_chr1_6 <-filter(a,chr<= 6)
a_chr1_6<-a_chr1_6 %>% select(CHR,start,end)
write.table(a_chr1_6,"/home/huanhuan/project/GTEx/output/Tissue_merge/hotspot_chr1_6/041_filter_length_6_1243_hotspot_chr1_6.bed",row.names=F,col.names=F,quote=F,sep="\t")
a_chr1_11 <-filter(a,chr<= 11)
a_chr1_11<-a_chr1_11 %>% select(CHR,start,end)
write.table(a_chr1_11,"/home/huanhuan/project/GTEx/output/Tissue_merge/hotspot_chr1_11/041_filter_length_6_1243_hotspot_chr1_11.bed",row.names=F,col.names=F,quote=F,sep="\t")
a<-a %>% select(CHR,start,end)
write.table(a,"/home/huanhuan/project/GTEx/output/Tissue_merge/041_filter_length_6_1243_hotspot.bed",row.names=F,col.names=F,quote=F,sep="\t")
