library(dplyr)


org<-read.tabel("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/output/random_segmention/01_count_number_and_length_of_hotspot_chr_in_Whole_Blood.txt.gz",sep="\t",header=T)

aa <-org%>%group_by(Cutoff,Length)%>%summarise(sum(Number))

write.table(aa,"/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/random_sample/aa.txt",col=F,row= F)