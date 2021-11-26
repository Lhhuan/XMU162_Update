ProcessBedGz <- function(file = NULL){
  name <- str_replace(file,"QTLbase_all_eQTL_","")  #replace
  name <- str_replace(name,".bed.gz","")  #replace
  name <- str_replace(name,"_narrow_peak","")  #replace
  print(name)
  if (name != "hic_loops_1" && name != "hic_loops_2" && name !="phastCons100way" ){
    print("klk")
    org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
    org2$name <-paste(org2$V6,org2$V7,org2$V8,sep=":")
    colnames(org2)[9]<-name
    org2$key <-paste(org2$V1,org2$V2,sep=":")
    org2<-org2[,-c(1:8)] #remove top 8 column
    #tmp <-left_join(tmp,org2,by = "key")
  } else if(name == "hic_loops_1" || name == "hic_loops_2" ){ #hi-C
    org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
    org2$name <-paste(org2$V6,org2$V7,org2$V8,org2$V12,sep=":")
    colnames(org2)[13]<-name
    org2$key <-paste(org2$V1,org2$V2,sep=":")
    org2<-org2[,-c(1:12)] #remove top 12 column
    #tmp <-left_join(tmp,org2,by = "key")
  }else{
    org2<-read.table(file,header = F,sep = "\t") %>% as.data.frame()
    colnames(org2)[9]<-"conservation"
    org2$key <-paste(org2$V1,org2$V2,sep=":")
    org2<-org2[,-c(1:8)]
    #tmp <-left_join(tmp,org2,by = "key")
  }
  gc()
  return(org2)
}

files2 <- files
org2 <- mclapply(files2, ProcessBedGz, mc.cores = 20)

for(i in 4:length(org2)){
  tmp <-left_join(tmp,org2[[i]],by = "key")
}

gc()

microbenchmark::microbenchmark(read.table("QTLbase_all_eQTL_CTCF.bed.gz",header = F,sep = "\t"),
                               read.delim("QTLbase_all_eQTL_CTCF.bed.gz",header = F,sep = "\t"),
                               read.delim2("QTLbase_all_eQTL_CTCF.bed.gz",header = F,sep = "\t"),
                               fread("QTLbase_all_eQTL_CTCF.bed.gz",header = F,sep = "\t"),
                               read.csv("QTLbase_all_eQTL_CTCF.bed.gz",header = F,sep = "\t"),
                               times = 1)

a <- read.delim2("QTLbase_all_eQTL_CTCF.bed.gz",header = F,sep = "\t")
b <- fread("QTLbase_all_eQTL_CTCF.bed.gz",header = F,sep = "\t")
c <- read.table("QTLbase_all_eQTL_CTCF.bed.gz",header = F,sep = "\t")


split -l 300000 QTLbase_all_eQTL_CTCF.bed.gz -d -a 3 QTLbase_all_eQTL_CTCF_ &&ls| grep QTLbase_all_eQTL_CTCF_ | xargs -n1 -i mv {} {}.bed.gz



# 
# ###
#
zcat QTLbase_all_eQTL_circRNA.bed.gz | split -l 300000 --filter='gzip > $FILE.bed.gz' &&ls | grep x | xargs -n1 -i mv {} QTLbase_all_eQTL_circRNA_{}

notationName <- list.files()
finalOrg <- lapply(notationName, function(i){
  files <- list.files(path = paste0("./",i), pattern = paste0(i,"_"))
  org <- mclapply(paste0(i,"/",files), 
                  function(file){
                    print(file)
                    ProcessBedGz(file)
                  }, mc.cores = 15)
  org
})
names(finalOrg) <- notationName

# fenpi
a <- tmp
b <- tmp

##1
for(i in 1:length(finalOrg[[1]])){ a <- left_join(a, finalOrg[[1]][[i]], by = "key") }
##2
for(i in 1:length(finalOrg[[1]])){ colnames(finalOrg[[1]][[i]])[1] <- "CTCF/QTLbase_all_eQTL_CTCF" }
    
x <- do.call(rbind,finalOrg[[1]])
b <- left_join(b, x, by = "key")

microbenchmark::microbenchmark(for(i in 1:length(finalOrg[[1]])){ a <- left_join(a, finalOrg[[1]][[i]], by = "key") },
                               b <- left_join(b, x, by = "key"),
                               times = 1)

b <- lineprof(for(i in 1:length(finalOrg[[1]])){ a <- left_join(a, finalOrg[[1]][[i]], by = "key") })


a[!duplicated(a$V5),]