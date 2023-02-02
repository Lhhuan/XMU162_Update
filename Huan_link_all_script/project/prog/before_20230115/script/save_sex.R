setwd("D:\\Huan_R\\prog\\data")
library(ggplot2)
library(dplyr)
library(stringr)
options(stringsAsFactors = FALSE)
all <-read.csv("all_data.2022-3-21-new.csv",na.strings = "")

colnames(all)[3]  <-"Sex"
sex <-all[,c("No","Sex")]
save(sex,file="sex.Rdata")
