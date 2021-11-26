# library(circlize)


# set.seed(999)
# n = 1000
# df = data.frame(factors = sample(letters[1:8], n, replace = TRUE),
#     x = rnorm(n), y = runif(n))

# circos.par("track.height" = 0.1)

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
# result1<-result[,c('data2.SNP','Chromosome','Position','reQTL','pQTL','edQTL','caQTL','hQTL','sQTL','mQTL','eQTL')]

final_file<-result[,2:ncol(result)]
colnames(final_file)[1] <- "chr"
colnames(final_file)[2] <- "end"
final_file$end<-as.integer(final_file$end)
final_file$start <-final_file$end -1
final_file$start<-as.integer(final_file$start)
final_file1<- final_file[,c(1,11,2,3,4,5,6,7,8,9,10)]
# final_file1<- final_file[,c(1,11,2,10)]
# colnames(final_file1)[4] <- "value1"
final_file1[is.na(final_file1)]=0
final_file1$chr <-as.character(paste("chr",final_file1$chr,sep=""))

# aaa<-sample(1:nrow(final_file1),150)
aaa<-sample(1:nrow(final_file1),1500)
bbb<-final_file1[aaa,]
library(circlize)
circos.clear()
# circos.initializeWithIdeogram()
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
col_fun = colorRamp2(c(0, 0.2, 0.4), c("green","red", "blue"))
circos.genomicHeatmap(bbb, col = col_fun, side = "inside", border = "white")
dev.off()
head(final_file1)



