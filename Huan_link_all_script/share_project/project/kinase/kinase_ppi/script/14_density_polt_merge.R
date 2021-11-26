library(ggplot2)
library(dplyr)
library(readr)
setwd("/home/huanhuan/project/kinase/kinase_ppi/script/")
di<-"~/project/kinase/kinase_ppi"

ran<-read.table(file.path(di, "output/all_random_rwr_result.txt"),header = T,sep = "\t") %>% as.data.frame()
ki<-read.table(file.path(di, "output/all_kinase_rwr_result.txt"),header = T,sep = "\t") %>% as.data.frame()
# random<-dplyr::rename(ran, Affinity=X0.7024283533) #重命名
names(ran)<-c("Affinity")
random<-subset(ran,Affinity<0.6) #把RWR的第一个去掉，即rwr的start 
# kinase<-dplyr::rename(ki,Affinity=X0.7023344719)
names(ki)<-c("Affinity") #修改列名
kinase<-subset(ki,Affinity<0.6) #把RWR的第一个去掉，即rwr的start 
random$class <-c("Random")
kinase$class <-c("Kinase")
rs<-bind_rows(random,kinase)
setwd("/home/huanhuan/project/kinase/kinase_ppi/figure/")

pdf("14_all_random_kinase_Affinity_density.pdf",height = 3.5,width = 5)
p<-ggplot(rs,aes(x=Affinity,colour= class))+geom_density(alpha=.2) + xlab("Affinity") #alpha设置透明度
p1<-p+theme(panel.grid =element_blank())+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                panel.background = element_blank(), axis.title.y = element_text(size = 13),
                                                axis.title.x = element_text(size = 13),
                                                axis.line = element_line(colour = "black")) #去背景
p1<-p1+ylab("Density")
p1<-p1+labs(colour = "Class") #修改图例名字
p1
dev.off()

pdf("14_all_random_kinase_Affinity_density_x_log.pdf",height = 3.5,width = 5)
p2<-p1+scale_x_log10()+xlab("log10(Affinity)") #对x轴进行log变换
p2
dev.off()

random<-ran[,1]
kinase<-ki[,1]
ks.test(random,kinase)
ks.test(random,kinase,alternative ="greater") #random 随机小于kinase
