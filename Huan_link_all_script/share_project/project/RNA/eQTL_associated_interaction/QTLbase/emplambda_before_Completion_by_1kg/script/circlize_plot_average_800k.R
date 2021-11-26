library(gcookbook)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
library(stringr)
library(circlize)
library(CMplot)
library(reshape2)
library(tidyverse)


setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/QTLbase/figure/")
org<-fread("../output/all_NHPoisson_emplambda_interval_1000_cutoff_7.3_win_800kb_bed.gz",header = T,sep = "\t") %>% as.data.frame()
colnames(org)[1] <- "chr"
colnames(org)[2] <- "start1"
colnames(org)[3] <- "end1"
colnames(org)[5] <- "start2"
colnames(org)[6] <- "end2"
colnames(org)[7] <- "xQTL"
colnames(org)[8] <- "emplambda"
org_need<-org[,c('chr','start1','end1','start2','end2','xQTL','emplambda')]
org_need$region<-paste(org_need$start1,org_need$end1,sep="_")
org_a<-org_need%>%group_by(xQTL,chr,region)%>%summarise(mean(emplambda))%>%as.data.frame()
colnames(org_a)[4] <- "mean_emplambda"
data<- spread(org_a,xQTL,mean_emplambda)
nrow(data)
#3398

org_all <-filter(org,xQTL =="QTL")
org_xqtl <-filter(org,xQTL !="QTL")
org_all$xQTL <-paste("ALL ",org_all$xQTL, sep = "")
data<-rbind(org_xqtl,org_all)
# data<-org
data$SNP <- paste0(data$chr,":",data$t)
data <- data[,c('SNP','emplambda','xQTL')]
data2 <- spread(data,xQTL,emplambda)
head(data2)

Position <- str_extract(data2$SNP,":(.*)")
Position <- gsub(":","",Position)
Chromosome <- str_extract(data2$SNP,"(.*):")
Chromosome <- gsub(":","",Chromosome)

result <- data.frame(data2$SNP,Chromosome,Position,data2[,2:ncol(data2)])


#--------------------------

final_file<-result[,c('Chromosome','Position','caQTL','hQTL','sQTL','mQTL','eQTL','ALL.QTL')]
colnames(final_file)[1] <- "chr"
colnames(final_file)[2] <- "end"
final_file$end<-as.integer(final_file$end)
final_file$start <-final_file$end -1
final_file1<- final_file[,c("chr","start","end",'caQTL','hQTL','sQTL','mQTL','eQTL','ALL.QTL')]
final_file1[is.na(final_file1)]=0
final_file1$chr <-as.character(paste("chr",final_file1$chr,sep=""))
unique_chr <-unique(final_file1$chr)
# unique_chr_n<-length(unique_chr)
rs<-data.frame()

# file_list=list()
qtl_types<-c('caQTL','hQTL','sQTL','mQTL','eQTL','ALL.QTL')


for(j in c(1:length(unique_chr))){
    file_chr=filter(final_file1,chr == unique_chr[j])
    max_pos = max(file_chr$end)
    min_pos =min(file_chr$end)   
    for(qtl in qtl_types ){
      file_list=filter(final_file1,chr == unique_chr[j])%>%dplyr::select(qtl,end)
      end_list = c(file_list$end)
      for()
      i=i+1
    }}


file_chr <- group_by(final_file1, chr)


library(circlize)
circos.clear()
# circos.initializeWithIdeogram()
pdf("circos_heatmap_plot.pdf")
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
col_fun = colorRamp2(c(-1, 0, 1), c("green", "black", "red"))
circos.genomicHeatmap(final_file1, col = col_fun, side = "inside", border = "white")

dev.off()
#----------------------------------





head(result)
 
result1<-result[,c('data2.SNP','Chromosome','Position','caQTL','hQTL','sQTL','mQTL','eQTL','ALL.QTL')]
CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
       memo="ylim_0.6_cex_0.005",dpi=300,cex= 0.005,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.6)

head(result1)

#---------------------------------------------------------------- #green #渐变色
threshold_c <-seq(from=0.1, to=0.6, by=0.1 )
mycol = colorRampPalette(c("#65FDF0", "#1D6FA3"),bias=1)(n=61) #渐变色盘
barplot(rep(1,61),col=mycol,border=NA,main="bias=1")
mycol<-mycol[61:1]
signal_c<-mycol[1:60]
col_c <-mycol[61]
CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",col=col_c,
       memo="ylim_0.6_cex_0.005_col_green",dpi=300,cex= 0.005,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.6,
       threshold=threshold_c,
        threshold.col=NULL,
        signal.col=signal_c,signal.cex=0.005)
#------------------------------------------------------------blue
mycol = colorRampPalette(c("#B6DAFD", "#4facfe"),bias=1)(n=61)
mycol<-mycol[61:1]
signal_c<-mycol[1:60]
col_c <-mycol[61]
CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",col=col_c,
        memo="ylim_0.6_cex_0.005_col_blue_2",dpi=300,cex= 0.005,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.6,
        threshold=threshold_c,
         threshold.col=NULL,
        signal.col=signal_c,signal.cex=0.005
       )






# result1<-result[,c('data2.SNP','Chromosome','Position','reQTL','pQTL','edQTL','caQTL','hQTL','sQTL','mQTL','eQTL','ALL.QTL')]

# CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
#        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
#        memo="ylim_0.6",dpi=300,cex= 0.01,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.6)


# CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
#        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
#        memo="ylim_0.75",dpi=300,cex= 0.01,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.75)



#------------------------------------------
# CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
#        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
#        memo="ylim_0.6_cex_0.005",dpi=300,cex= 0.005,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.6)


# CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
#        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
#        memo="ylim_0.55",dpi=300,cex= 0.01,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.55)

# CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
#        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
#        memo="ylim_0.55_cex_0.005",dpi=300,cex= 0.005,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.55)

# CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
#        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
#        memo="ylim_0.6",dpi=300,cex= 0.01,file.output=T,verbose=TRUE,LOG10=FALSE,ylim=0.6)

# CMplot(result1,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=1,cir.legend=TRUE,
#        outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
#        memo="",dpi=300,cex= 0.01,file.output=T,verbose=TRUE,LOG10=FALSE)
#--------------------------------------------------------------------------





















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
CMplot(result,plot.type="c",chr.labels=paste("Chr",c(1:22),sep=""),r=0.6,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
       memo="",dpi=300,cex= 0.08,file.output=T,verbose=TRUE)
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
