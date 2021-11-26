# R >= 3.1.1 is needed

if (!require("NHPoisson")) {
    install.packages("NHPoisson")
    library("NHPoisson")
}


library(ggplot2)
library(dplyr)
library(Rcpp)
library(readxl)
library(stringr)
library(gcookbook)
library(gridExtra)
library(ggpubr)

# setting a random seed for reproducibility
set.seed(1138)

# Data preparation
# Variables
# dia: Postion of the day in the year, from 121 (1st of May) to 253 (30th of September).
# mes: Month of the year, from 5 to 9.
# ano: Year, from 1951 to 2004.
# diames: Position of the day in the month, from 1 to 30 or 31.
# Tx: Daily maximum temperature.
# Tn: Daily minimum temperature.
# Txm31: Local maximum temperature signal. Lowess of Tx with a centered window of 31 days.
# Txm15: Local maximum temperature signal. Lowess of Tx with a centered window of 15 days.
# Tnm31: Local minimum temperature signal. Lowess of Tn with a centered window of 31 days.
# Tnm15: Local minimum temperature signal. Lowess of Tn with a centered window of 15 days.
# TTx: Long term maximum temperature signal. Lowess of Tx with a centered 40% window.
# TTn: Long term minimum temperature signal. Lowess of Tn with a centered 40% window.

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



data("BarTxTn", package = "NHPoisson")
dateB <- cbind(BarTxTn$ano, BarTxTn$mes, BarTxTn$diames)
BarEv <- POTevents.fun(T = BarTxTn$Tx, thres = 318, date = dateB)

# This function calculates the characteristics of the extreme events of a series(xi) defined using a peak
# over threshold (POT) method with an extreme threshold. The initial and the maximum intensity
# positions, the mean excess, the maximum excess and the length of each event are calculated.


Px <- BarEv$Px
inddat <- BarEv$inddat

tB <- BarTxTn$ano + rep(c(0:152)/153, 55)
emplambdaB <- emplambda.fun(posE = Px, inddat = inddat, t = tB, lint = 153, tit = "Barcelona")

tmp2<-cbind(emplambda=emplambdaB$emplambda, t=emplambdaB$t) %>% as.data.frame()
write.table(tmp2,"emplambdaB.txt",row.names = F, col.names = T,quote =F,sep="\t")

mod1Bind <- fitPP.fun(covariates = NULL, posE = BarEv$Px, inddat = BarEv$inddat, 
  tit = "BAR  Intercept ", start = list(b0 = 1))
covB <- cbind(cos(2 * pi * BarTxTn$dia/365), sin(2 * pi * BarTxTn$dia/365), BarTxTn$Txm31, 
  BarTxTn$Tnm31, BarTxTn$TTx, BarTxTn$TTn)
dimnames(covB) <- list(NULL, c("Cos", "Sin", "Txm31", "Tnm31", "TTx", "TTn"))
aux <- stepAICmle.fun(ImlePP = mod1Bind, covariatesAdd = covB, startAdd = c(1, -1, 
  0, 0, 0, 0), direction = "both")


modB.1 <- fitPP.fun(covariates = NULL, posE = Px, inddat = inddat, tit = "BARCELONA Tx; Intercept", 
  start = list(b0 = 1), dplot = FALSE, modCI = FALSE)