# R >= 3.1.1 is needed

if (!require("NHPoisson")) {
    install.packages("NHPoisson")
    library("NHPoisson")
}

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

# Fitting the model
  #------------ fitPP.fun
  #This function fits by maximum likelihood a NHPP where the intensity λ(t) is formulated as a function of covariates. It also calculates and plots approximate confidence intervals for λ(t).
mod1Bind <- fitPP.fun(covariates = NULL, posE = BarEv$Px, inddat = BarEv$inddat, 
  tit = "BAR  Intercept ", start = list(b0 = 1))
covB <- cbind(cos(2 * pi * BarTxTn$dia/365), sin(2 * pi * BarTxTn$dia/365), BarTxTn$Txm31, 
  BarTxTn$Tnm31, BarTxTn$TTx, BarTxTn$TTn)
dimnames(covB) <- list(NULL, c("Cos", "Sin", "Txm31", "Tnm31", "TTx", "TTn"))
aux <- stepAICmle.fun(ImlePP = mod1Bind, covariatesAdd = covB, startAdd = c(1, -1, 
  0, 0, 0, 0), direction = "both")


modB.1 <- fitPP.fun(covariates = NULL, posE = Px, inddat = inddat, tit = "BARCELONA Tx; Intercept", 
  start = list(b0 = 1), dplot = FALSE, modCI = FALSE)


covB <- cbind(cos(2 * pi * BarTxTn$dia/365), sin(2 * pi * BarTxTn$dia/365))
modB.2 <- fitPP.fun(covariates = covB, posE = Px, inddat = inddat, tit = "BARCELONA Tx; Cos, Sin", 
  start = list(b0 = -100, b1 = 1, b2 = 1), modSim = TRUE, dplot = FALSE, modCI = FALSE)
aux <- testlik.fun(ModG = modB.2, ModR = modB.1) #testlik.fun: Likelihood ratio test to compare two nested models

covB <- cbind(cos(2 * pi * BarTxTn$dia/365), sin(2 * pi * BarTxTn$dia/365), cos(4 * 
  pi * BarTxTn$dia/365), sin(4 * pi * BarTxTn$dia/365))
modB.3 <- fitPP.fun(covariates = covB, posE = Px, inddat = inddat, tit = "BARCELONA Tx; Cos, Sin, Cos2, Sin2", 
  start = list(b0 = -100, b1 = 1, b2 = 1, b3 = 1, b4 = 1), modSim = TRUE, dplot = FALSE, 
  modCI = FALSE)
aux <- testlik.fun(ModG = modB.3, ModR = modB.2)


covB.final <- cbind(cos(2 * pi * BarTxTn$dia/365), sin(2 * pi * BarTxTn$dia/365), 
  BarTxTn$Txm31)
dimnames(covB.final) <- list(NULL, c("Cos", "Sin", "Txm31"))
modB.final <- fitPP.fun(covariates = covB.final, posE = Px, inddat = inddat, tim = tB, 
  tit = "BARCELONA Tx; Cos, Sin, Txm31", start = list(b0 = -100, b1 = 1, b2 = 1, 
    b3 = 0))

aux <- LRTpv.fun(modB.final)

summary(modB.final)

confint(modB.final)

confintAsin.fun(modB.final)

# Model validation

  #graphrate.fun: Plot Of The Fitted And Empirical Occurrence Rate
posEHB <- transfH.fun(modB.final)$posEH
resB <- unifres.fun(posEHB)
graphresU.fun(unires = resB$unires, posE = modB.final@posE, Xvariables = cbind(covB.final, 
  BarTxTn$dia), namXv = c("cos", "sin", "Txm31", "summer day index"), tit = "BARCELONA; cos, sin, Txm31", 
  addlow = FALSE)

ResDB <- CalcResD.fun(mlePP = modB.final, lint = 153)

qqnorm(ResDB$RawRes)

graphrate.fun(ResDB, tit = "BARCELONA; cos, sin, Txm31")

covBtemp <- cbind(BarTxTn$Txm31, BarTxTn$TTx)
aux <- graphResCov.fun(mlePP = modB.final, nint = 50, X = covBtemp, namX = c("Txm31", 
  "TTx"), tit = "BARCELONA; cos, sin, Txm31", typeRes = "Raw", indgraph = TRUE)

graphres.fun(objres = ResDB, typeRes = "Pearson", addlow = TRUE, indgraph = TRUE, 
  tit = "BARCELONA; cos, sin, Txm31")
aux <- resQQplot.fun(nsim = 100, objres = ResDB, covariates = covB.final,
  tit = "BARCELONA; cos, sin, Txm31", fixed.seed = 123)


aux <- globalval.fun(mlePP = modB.final, lint = 153, covariates = covB.final, typeI = "Disjoint", 
  typeResLV = "Raw", nintLP = 50, tit = "BARCELONA; cos, sin, Txm31", Xvar = covBtemp, 
  namXvar = c("Txm31", "TTx"), resqqplot = TRUE, fixed.seed = 123)


# simulation based inference

aux <- GenEnv.fun(lambda = modB.final@lambdafit[8263:8415], fun.name = "length", 
  nsim = 1000, fixed.seed = 123)


aux <- GenEnv.fun(lambda = modB.final@lambdafit[8263:8415], fun.name = "posk.fun", 
  fun.args = 1, nsim = 1000, fixed.seed = 123)

L <- NULL
for (i in c(1:1000)) L[i] <- length(simNHP.fun(lambda = modB.final@lambdafit[8263:8415],
  fixed.seed = (123 + i))$posNH)
P0 <- sum(L == 0)/1000
P0 
