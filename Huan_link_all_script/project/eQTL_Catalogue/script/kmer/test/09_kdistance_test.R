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

setwd("/home/huanhuan/project/eQTL_Catalogue/script/kmer/figure/")
load("sig_kmer_0_100.Rdata")

.kdist <- function(x, from, to, seqlengths, k) {
    .Call('_kmer_kdist', PACKAGE = 'kmer', x, from, to, seqlengths, k)
}

# Check if object is DNA
.isDNA <- function(x){
  if(inherits(x, "DNAbin")){
    return(TRUE)
  }else if(inherits(x, "AAbin")){
    return(FALSE)
  }else if(mode(x) == "character"){
    return(FALSE)
  }else if(mode(x) == "raw"){
    return(all(x %in% as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48,
                               224, 176, 208, 112, 240, 4, 2))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x, use.names = FALSE) %in%
                   as.raw(c(136, 72, 40, 24, 192, 160, 144, 96, 80, 48,
                                         224, 176, 208, 112, 240, 4, 2))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

# Check if object is amino acid sequence
.isAA <- function(x){
  if(inherits(x, "AAbin")){
    return(TRUE)
  }else if(inherits(x, "DNAbin")){
    return(FALSE)
  }else if(mode(x) == "character"){
    return(FALSE)
  }else if(mode(x) == "raw"){
    return(all(x %in% as.raw(c(65:90, 42, 45, 63))))
  }else if(mode(x) == "list"){
    if(length(x) > 0){
      return(all(unlist(x, use.names = FALSE) %in% as.raw(c(65:90, 42, 45, 63))))
    }else{
      return(FALSE)
    }
  }else{
    return(FALSE)
  }
}

.alphadetect <- function(sequences, residues = NULL, gap = "-",
                         endchar = "?"){
  if(identical(toupper(residues), "RNA")){
    residues <- c("A", "C", "G", "U")
  }else if(.isDNA(sequences) | identical(toupper(residues), "DNA")){
    residues <- c("A", "C", "G", "T")
  }else if(.isAA(sequences) | identical(residues, "AA") |
           identical (toupper(residues), "AMINO")){
    residues <- LETTERS[-c(2, 10, 15, 21, 24, 26)]
  }
  else if(is.null(residues)){
    residues <- sort(unique(as.vector(unlist(sequences, use.names = FALSE))))
    if(!is.null(gap)) residues <- residues[residues != gap]
    if(!is.null(endchar)) residues <- residues[residues != endchar]
  }else{
    if(!is.null(gap)) residues <- residues[residues != gap]
    if(!is.null(endchar)) residues <- residues[residues != endchar]
  }
  if(!(length(residues) > 0)){# & mode(residues) == "character")){
    stop("invalid residues argument")
  }
  return(residues)
}


.unalign <- function(x, gap = "-"){
  DNA <- .isDNA(x)
  AA <- .isAA(x)
  gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
  if(is.list(x)){
    if(length(x) == 1){
      tmpname <- names(x)
      x <- x[[1]]
      if(is.null(dim(x))){
        x <- matrix(x, nrow = 1)
        rownames(x) <- tmpname
      }
    }
  }
  res <- vector(mode = "list", length = nrow(x))
  for(i in 1:nrow(x)) res[[i]] <- x[i, x[i, ] != gap, drop = TRUE]
  if(AA){
    res <- lapply(res, unclass)
    class(res) <- "AAbin"
  }else if(DNA){
    res <- lapply(res, unclass)
    class(res) <- "DNAbin"
  }
  names(res) <- rownames(x)
  return(res)
}


library(ape)
data(woodmouse)
x <-woodmouse
k = 5
residues = NULL
gap = "-"
named = TRUE
compress = TRUE
encode = FALSE
DNA <- .isDNA(x)
if(DNA) class(x) <- "DNAbin" else if(AA) class(x) <- "AAbin"
residues <- .alphadetect(x, residues = residues, gap = gap)
gap <- if(AA) as.raw(45) else if(DNA) as.raw(4) else gap
if(is.matrix(x)) x <- .unalign(x, gap = gap)
nseq <- length(x)
if(is.null(names(x))) names(x) <- paste0("S", 1:nseq)
seqalongx <- seq_along(x)
kcounts <- kcount(x, k = k, residues = residues, gap = gap, compress = compress)
seqlengths <- apply(kcounts, 1, sum) + k - 1



seqalongx <- seq_along(woodmouse)
kcounts <- kcount(woodmouse, k = 3)
k=3
seqlengths <- apply(kcounts, 1, sum) + k - 1

d <- kdist(kcounts, from = seqalongx - 1, to = seqalongx - 1,
                        seqlengths = seqlengths, k = k)


dat <-as.matrix(Sorg)
seqlengths <- apply(kcounts, 1, sum) + k - 1