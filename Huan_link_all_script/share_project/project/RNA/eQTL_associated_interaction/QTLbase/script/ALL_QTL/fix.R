
CHR1<-filter(rs, chr==1)
CHR1 <- CHR1[order(CHR1$t), ]


gz1 <- gzfile("sorted_chr1.txt.gz", "w")
#write.csv(df1, gz1)
write.table(CHR1,gz1,row.names = F, col.names = T,quote =F,sep="\t")
close(gz1)



CHR2<-filter(rs, chr==2)
CHR2 <- CHR2[order(CHR2$t), ]
# write.table(CHR2,"sorted_chr2.txt",row.names = F, col.names = T,quote =F,sep="\t")


gz1 <- gzfile("sorted_chr2.txt.gz", "w")
#write.csv(df1, gz1)
write.table(CHR2,gz1,row.names = F, col.names = T,quote =F,sep="\t")
close(gz1)


library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

# R >= 3.1.1 is needed

if (!require("NHPoisson")) {
    install.packages("NHPoisson")
    library("NHPoisson")
}

# setting a random seed for reproducibility
set.seed(1138)

POTevents.fun <- function (T, thres, date = NULL)
{
    if (is.null(date))
        date <- c(1:length(T))
    date <- as.matrix(date)
    if ((is.null(date) == FALSE) & (dim(date)[1] != length(T)))
        stop("T and date must have the same number of observations")
    exc <- (T > thres)
    inrachtx <- (diff(c(0, exc)) == 1)
    Pi <- c(1:length(T))[inrachtx == 1]
    numerachtx <- cumsum(inrachtx)[exc == 1]
    intentx <- (T - thres)[exc == 1]
    Im <- tapply(intentx, INDEX = numerachtx, FUN = mean)
    Ix <- tapply(intentx, INDEX = numerachtx, FUN = max)
    L <- tapply(intentx, INDEX = numerachtx, FUN = length)
    Px <- Pi + tapply(intentx, INDEX = numerachtx, FUN = which.max) - 1
    inddat <- 1 - exc
    inddat[Px] <- 1
    datePi <- date[Pi, ]
    datePx <- date[Px, ]
    cat("Number of events: ", length(Im), fill = TRUE)
    cat("Number of excesses over threshold", thres, ":", sum(exc),
        fill = TRUE)
    return(list(Pi = Pi, datePi = datePi, Px = Px, datePx = datePx,
        Im = Im, Ix = Ix, L = L, inddat = inddat, T = T, thres = thres,
        date = date))
}


a<-data.frame()
a<-data.frame(1:3000,0.05)
a[1500,2] <- 5e-9
a[1498,2] <- 5e-9


cutoff =7.3
j =5

colnames(a)[1] <-"SNP_pos"
colnames(a)[2] <-"min_p"
BarEv <- POTevents.fun(T = -log10(a$min_p), thres = cutoff, date = a$SNP_pos)
tB <-a$SNP_pos
emplambdaB <- emplambda.fun(posE = BarEv$Px,inddat=BarEv$inddat, t = tB, lint = j, tit = "title")


LAL <-emplambdaB$emplambda


indice[1:5]

sum(indice[iini:ifin])/sum(inddat[iini:ifin]