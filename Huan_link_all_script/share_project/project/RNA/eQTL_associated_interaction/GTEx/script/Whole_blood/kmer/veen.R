library(VennDiagram)
A = 1:150
B = c(121:170,300:320)
C = c(20:40,141:200)
Length_A<-length(A)
Length_B<-length(B)
Length_C<-length(C)

setwd("/home/huanhuan/project/RNA/eQTL_associated_interaction/GTEx/script/Whole_blood/kmer/")

Length_AB<-length(intersect(A,B))

Length_BC<-length(intersect(B,C))

Length_AC<-length(intersect(A,C))

Length_ABC<-length(intersect(intersect(A,B),C))

T<-venn.diagram(list(A=A,B=B),filename=NULL,lwd=1,lty=2, ,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'),rotation.degree=90)
grid.draw(T)

T<-venn.diagram(list(A=A,B=B,C=C),filename=NULL,lwd=1,lty=2,col=c('red','green','blue') ,fill=c('red','green','blue'),cat.col=c('red','green','blue'),reverse=TRUE)

grid.draw(T)