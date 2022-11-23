library(ggplot2)
library(Rcpp)
library(readxl)
library(dplyr)
library(stringr)
library(ggpubr)
library(gridExtra)
library(ggpval)
library(Seurat)
library(reshape2)
library(R.utils)
library(mclust)
library(kmer)

.kdist <- function(x, from, to, seqlengths, k) {
    .Call('_kmer_kdist', PACKAGE = 'kmer', x, from, to, seqlengths, k)
}

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/figure/")
load("sig_kmer_0_100.Rdata")
seqalongx=seq_along(1:nrow(Sorg))
k=6
kcounts <-as.matrix(Sorg)
seqlengths =apply(kcounts, 1, sum) + k - 1


d <- .kdist(kcounts, from = seqalongx - 1, to = seqalongx - 1,
                        seqlengths = seqlengths, k = k)

save(d,file ="08_kdistance.Rdata")