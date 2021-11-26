library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(circlize)


circos.clear()

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")
org<-fread("../output/NHPoisson_emplambda_interval_1000_cutoff_7.3_all_merge_eQTL_caQTL_edQTL_hQTL_mQTL_pQTL_reQTL_sQTL.txt.gz",header = T,sep = "\t") %>% as.data.frame()
#-------------------------------------

org_eQTL <-filter(org,xQTL == "eQTL")
org_eQTL$chr<- as.factor(org_eQTL$chr)
org_eQTL_a <- org_eQTL[,1:3]
colnames(org_eQTL_a)[2] <- 'start'
org_eQTL_a$end <-org_eQTL_a$start-1
org_eQTL_a=org_eQTL_a[,c(3,2,4,1)]
# org_eQTL_a$chr<-paste("chr",org_eQTL_a$chr,"-",org_eQTL_a$start,"-",,org_eQTL_a$start)
# input_file <-paste("../../../output/ALL_", xQTL, "/NHPoisson_emplambda_interval_",j,"cutoff_",cutoff,"_all_",xQTL,".txt", sep = "")
#---------------------------------------------------
f1 = function() {
  circos.par(gap.after = 2, start.degree = 90)
  circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22), 
                                plotType = c("ideogram", "labels"), ideogram.height = 0.03)
}

circos.trackPoints(org_eQTL$chr, org_eQTL$t, org_eQTL$emplambda, col = col, pch = 16, cex = 0.5)


#---------------------------------------
rs1 <- data.frame()
xQTLs<-unique(org$xQTL)
for(QTL in xQTLs ){
  ALL_qtl <-filter(org,xQTL==QTL)
  for(i in c(1:22)){
    all_ALL_qtl <-filter(ALL_qtl,chr==i)
    tmp<-head(all_ALL_qtl,50)
    rs1 <- bind_rows(rs1,tmp)}}

write.table(rs1,"../output/test_demo.txt",row.names = F, col.names = T,quote =F,sep="\t")


#----------------------------------------------


library(CMplot)
library(reshape2)
library(stringr)
library(tidyverse)
data <- read.table("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/output/test_demo.txt",head=T)
# data <- read.table("../output/test_demo.txt",head=T)
data$SNP <- paste0(data$chr,":",data$t)
data <- data[,c('SNP','emplambda','xQTL')]
data2 <- spread(data,xQTL,emplambda)
head(data2)

Position <- str_extract(data2$SNP,":(.*)")
Position <- gsub(":","",Position)
Chromosome <- str_extract(data2$SNP,"(.*):")
Chromosome <- gsub(":","",Chromosome)

result <- data.frame(data2$SNP,Chromosome,Position,data2[,2:ncol(data2)])
head(result)
CMplot(result,type="p",plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=0.6,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
       memo="",dpi=300,cex= 0.08,file.output=T,verbose=TRUE,width=10,height=10)
#------------------------------




#-----------------------------------
org$chr<- as.factor(org$chr)
circos.initialize(factors = org$chr, x = org$t )
circos.trackPlotRegion(factors = org$chr, y=org$emplambda, panel.fun = function(x, y) {
  circos.axis(labels.cex=0.5, labels.font=1, lwd=0.8)
})
#-------------------
circos.track(factors = org$chr, y = org$emplambda,
panel.fun = function(x, y) {
circos.text(x = get.cell.meta.data("xcenter"),
y = get.cell.meta.data("cell.ylim")[2] + uy(5, "mm"),
labels = get.cell.meta.data("sector.index"))
circos.axis(labels.cex = 0.5)
})


#--------------------------------------------------------


pdf("circlize_emplambda_interval_1000_cutoff_7.3.pdf",height = 8,width = 8)

org_initialize <- unique(org[,1:3])
circos.initialize(factors = org_initialize$chr, x = org_initialize$t )

circos.trackPlotRegion(factors = org_initialize$chr, y=org_initialize$emplambda, panel.fun = function(x, y) {
  circos.axis(labels.cex=0.5, labels.font=1, lwd=0.8)
})


#--------------------------------------------
bbb<-fread("../output/test_a_aaa.txt",header = T,sep = "\t") %>% as.data.frame()
colnames(bbb)[2] <- 'start'
bbb$end<-bbb$start
bbb=bbb[,c(3,2,5,1,4)]
#--------------------------------------
bbb<-fread("../output/test_a.txt",header = T,sep = "\t") %>% as.data.frame()
colnames(bbb)[2] <- 'start'
bbb$end<-bbb$start
bbb=bbb[,c(3,2,5,1,4)]

#--------------------------------------------
colnames(org)[2] <- 'start'
org$end<-org$start
org=org[,c(3,2,5,1,4)]
circos.clear()
circos.initializeWithIdeogram(org, sort.chr = FALSE)
#---------------------------------------------------

org_eQTL <-filter(org,xQTL == "eQTL")



circos.par("track.height" = 0.1)
circos.initialize(factors = org$chr,x = org$t)
org_eQTL <-filter(org,xQTL == "eQTL")


circos.track(factors = org_eQTL$chr, y = org_eQTL$emplambda,
panel.fun = function(x, y) {
circos.text(x = get.cell.meta.data("xcenter"),
y = get.cell.meta.data("cell.ylim")[2] + uy(5, "mm"),
labels = get.cell.meta.data("sector.index"))
circos.axis(labels.cex = 0.5)
})

circos.track(factors = org$chr, y = org$emplambda,
panel.fun = function(x, y) {
circos.text(x = get.cell.meta.data("xcenter"),
y = get.cell.meta.data("cell.ylim")[2] + uy(5, "mm"),
labels = get.cell.meta.data("sector.index"))
circos.axis(labels.cex = 0.5)
})



col <-

rep(c("skyblue", "red"), 4)

circos.trackPoints(factors = df$factors,

x = df$x, y = df$y,

col = col, pch = 16, cex = 0.5)

circos.text(x = -1, y = 0.5, labels = "text",

sector.index = "a", track.index = 1)




pdf("histogram_count_number_by_emplambda_interval_1000_cutoff_7.3.pdf",height = 6,width = 6)
p1<-ggplot(org, aes(x =emplambda, y=number,fill =xQTL)) +geom_histogram(position="identity",alpha =0.4) 
print(p1)
dev.off()

for(i in c(1:22)){
    png("histogram_count_number_by_emplambda_interval_1000_cutoff_7.3.pdf",height = 1000,width = 2000)
    p1<-ggplot(org2, aes(x =t, y=emplambda,colour =xQTL)) +geom_point(size=0.1) 

    p2<-p1+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.title.y = element_text(size = 10),axis.title.x = element_text(size = 10),axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10),axis.text.x = element_text(size=10)) 
    p2<-p2+xlab("Position") + ylab("Emplambda") 
    p2<-p2+ggtitle(title_name)
    p2<-p2+theme(legend.text =element_text(size =20),legend.title =element_text(size =20),plot.title=element_text(size =20))
    print(p2)
    dev.off()
    print(i)
}
