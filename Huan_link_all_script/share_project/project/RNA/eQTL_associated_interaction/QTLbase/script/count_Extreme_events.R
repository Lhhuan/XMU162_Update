library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

pwd_dir = "/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/ALL_QTL"
setwd(pwd_dir)
org<-read.table("../../output/merge_QTL_all_QTLtype_pop.txt.gz",header = T,sep = "\t") %>% as.data.frame()

org1<-org[,c("SNP_chr","SNP_pos","Pvalue","QTL_type")]
org1$SNP_ID<-paste(org1$SNP_chr,org1$SNP_pos,sep="_")

#-------------------------------------------------density

p.theme<-theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 10),
                                                axis.title.x = element_text(size = 10),
                                                axis.line = element_line(colour = "black"))

pdf("p_value_density.pdf",height = 5,width = 5)
p <-ggplot (org,aes(x=-log10(Pvalue)))+ 
    geom_density(alpha=.2) + 
    xlab("Emplambda")+
    p.theme+
    ylab("Density")
print(p)
dev.off()



#----------------------------------count extreme events
unique_snp_id <-org1$SNP_ID %>%unique()
org_p <- org1%>%dplyr::select(SNP_ID,Pvalue)%>%unique()
org_pg <- group_by(org_p, SNP_ID)
org_min_P<-summarise(org_pg,min_p = min(Pvalue))%>%as.data.frame()
all_org_n <-nrow(org_min_P)
#11002442

org_extreme<- filter(org_min_P,-log10(min_p)>7.3)
org_extreme_n <-nrow(org_extreme)
#5774258
#------------------

#
org_all<-read.table("../../data/all.txt.gz",header = T,sep = "\t") %>% as.data.frame()
max(org_all$Pvalue) 
#0.04998182
min(org_all$Pvalue)
#4.940656e-324
#------------------------venn
library(VennDiagram)

org_QTL<-org1[,c("SNP_ID","QTL_type")]
org_QTL_unique <-unique(org_QTL)
# org_QTL_unique <- org1%>%dplyr::select(SNP_ID,QTL_type)%>%unique()
unique_QTLs <- unique(as.character(org_QTL_unique$QTL_type))

i=1
name_list=list()
file_list=list()
# for(qtl in unique_QTLs ){
for(j in c(1:length(unique_QTLs))){
    file_list[i]=filter(org_QTL_unique,QTL_type ==unique_QTLs[j])%>%dplyr::select(SNP_ID)
    name_list[i]=unique_QTLs[j]
    i=i+1}
# my_col =c('#f6d365','#fda085','#fbc2eb','#a6c1ee','#84fab0','#8fd3f4','#d4fc79','#96e6a1','#ff8177',
#       '#b12a5b','#43e97b','#30cfd0','#d299c2')
my_col =c('#f6d365','#fda085','#fbc2eb','#a6c1ee','#84fab0')
# my_col =c('#f6d365','#fda085','#fbc2eb','#a6c1ee','#84fab0','#8fd3f4','#d4fc79')
# T<-venn.diagram(list(caQTL=file_list[1],cerQTL=file_list[2],eQTL=file_list[3],edQTL=file_list[4],hQTL=file_list[5],lncRNAQTL=file_list[6],mQTL=file_list[7],metaQTL=file_list[8],miQTL=file_list[9],
# pQTL=file_list[10],reQTL=file_list[11],riboQTL=file_list[12],sQTL =file_list[13]),filename=NULL,col=my_col,fill=my_col)

T<-venn.diagram(list(caQTL=file_list[1],eQTL=file_list[3],hQTL=file_list[5],mQTL=file_list[7],sQTL =file_list[13]),filename=NULL,col=my_col,fill=my_col)
grid.draw(T)
