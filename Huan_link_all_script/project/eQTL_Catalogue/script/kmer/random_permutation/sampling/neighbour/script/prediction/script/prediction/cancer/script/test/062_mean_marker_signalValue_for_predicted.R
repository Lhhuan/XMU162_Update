library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)
library(tibble)
library(R.utils)
library(data.table)


setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/cancer/output/")
catf <- read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/ca_tf_TCGA_sample.txt",header = T,sep = "\t") %>% as.data.frame()
catf <- unique(catf[,c(2,1)])
catf$marker <-gsub("Human_CHROMATIN_Accessibility","CHROMATIN_Accessibility",catf$marker)
catf$marker <-gsub("Human_FACTOR","TFBS",catf$marker)

histone <- read.table("/share/Projects/huanhuan/project/RNA/eQTL_associated_interaction/annotation/annotation_data/cistromeDB/cancer_cell/TCGA_mark.txt",header = F,sep = "\t") %>% as.data.frame()
colnames(histone)<- c("TCGA","marker")

dat <- bind_rows(catf,histone)
dat$marker_path <- paste0(dat$TCGA,"/",dat$marker)

# markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","CTCF","TFBS","CHROMATIN_Accessibility")
# markers = dat$marker_path


predcited_signal<-function(marker= NULL){
    file_name <-paste0(marker,"_warm_hotspot_region_win4518_large_than6.bed.gz")
    if(file_test("-f",file_name)){
        org <- read.table(file_name,header = F,sep = "\t") %>% as.data.frame()
        colnames(org)[c(1:3,7)] <- c("Chr","start","end","signalvalue")
        dat <- org%>%group_by(Chr,start,end)%>%summarize(mean_signalvalue =mean(signalvalue),median_signalvalue =median(signalvalue),max_signalvalue = max(signalvalue))%>%as.data.frame()
        output_file <-paste0(marker,"_max_mean_median_warm_hotspot_region_win4518_large_than6.bed")
        write.table(dat,output_file,row.names = F, col.names = T,quote =F,sep="\t")
        print(marker)
        system(paste("gzip",output_file,sep=" "))
        return(marker)
    }
}

tmp <- lapply(dat$marker_path,predcited_signal)
rs <- do.call(rbind,tmp)
write.table(rs,"062_cancer_marker_list.txt",row.names = F, col.names = F,quote =F,sep="\t")

#=============================================================================================computer_matrix

# markers = c("H3K4me1","H3K4me3","H3K9ac","H3K9me3","H3K27ac","H3K27me3","H3K36me3","TFBS","CHROMATIN_Accessibility")

# unique_cancer <- unique(dat$TCGA)


blood <- fread("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/Blood/output/07_predicted_warm_hotspot_region_dataset.txt.gz",header = T,sep = "\t") %>% as.data.frame()
kmer_gc <- fread("/home/huanhuan/project/eQTL_Catalogue/script/kmer/random_permutation/sampling/neighbour/script/prediction/script/prediction/warm_hotspot/07_predicted_warm_hotspot_region_gc_kmerdataset.txt.gz",header = T,sep = "\t") %>% as.data.frame()
# kmer_gc[,402:410] <-NA 
# colnames(kmer_gc)[402:410] <- colnames(blood)[402:410] 


for (cancer in unique(dat$TCGA)){
    print("start")
    fdat <- kmer_gc 
    for(i in 402:410){
        # cat(i, "\n")
        marker <-colnames(blood)[i]
        file_name <- paste0(cancer,"/",marker,"_max_mean_median_warm_hotspot_region_win4518_large_than6.bed.gz")
        if(file_test("-f",file_name)){
            org <- read.table(file_name,header = T,sep = "\t") %>% as.data.frame()
            # print("here \n")
            fdat <- left_join(fdat,org[,1:4],by=c("Chr","start","end"))
            colnames(fdat)[i] <- marker
            
        }
        else{
            fdat$marker <- NA 
            colnames(fdat)[i] <- marker
        }
    }
    fdat[is.na(fdat)] <-0
    out_file <- paste0(cancer,"/warm_hotspot_region_win4518_large_than6_gc_kmer_marker.bed")
    write.table(fdat,out_file,row.names = F, col.names = T,quote =F,sep="\t")
    gzip(out_file)
    print(cancer)
    print("finish")
}







