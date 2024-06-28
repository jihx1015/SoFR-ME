#--- Sample code
#==================================================================
#---- Source all the R codes in r (from SoFRBIV_Main.R)
#==================================================================

#source("~/Library/CloudStorage/OneDrive-IndianaUniversity/StatMedSubmissionFIle/SoFRBIV_Main.R")
source("~/Library/CloudStorage/OneDrive-IndianaUniversity/StatMedSubmissionFIle/RCodeSIMSub.R")
#==================================================================
#--- We will simulate Y, Z, W, and M with n = 300
#==================================================================

set.seed(0)
Data = DataSimFunc()


#==================================================================
#--- Plot the Data
#==================================================================
par(mfrow=c(2,2))
library(MASS)
hist(Data$Y)

Y = Data$Y
W = Data$Data$W
M = Data$Data$M
Z = Data$Z
a = Data$a

# Nsim = 100; sig02.del = 1.0; tauval = 1; 
# alpx = 1; Gx = 3; Ge = 3; Gx0 = TRUE; BasFun = "Bspline"; PsiAccept = 0.05; lag = 1

#==================================================================
#--- Plot W
#==================================================================
matplot(Data$a, t(Data$Data$W), type = "l", col="lightgrey", lwd=1)
lines(Data$a, colMeans(Data$Data$W), col = "hotpink") ## Plot the mean

#==================================================================
#--- Plot M
#==================================================================
matplot(Data$a, t(Data$Data$M), type = "l", col="lightgrey", lwd=1)
lines(Data$a, colMeans(Data$Data$M), col = "hotpink") ## Plot the mean


#--- Sample estimation

SumRes4.pr.752ax <- BayesFIVMult(Nsim = 10000,
                                 Y = c(Data$Y), W = Data$Data$W, M = Data$Data$M, Z = Data$Z , PsiAccept = 1,
                                 sig2.propdelt = 0.75, sig2.priordelt = 1.5, sig02.del = 1)#sig2.propdelt = 3.75, sig2.priordelt = 5.0

mean(SumRes4.pr.752ax$AcceptDelt)



SumRes4.pr.752x <- BayesFIVSing(Nsim = 10000,
                                Y = c(Data$Y), W = Data$Data$W, M = Data$Data$M, Z = Data$Z , PsiAccept = 1,
                                sig2.propdelt = 0.75, sig2.priordelt = 0.8, sig02.del = 1.2)

mean(SumRes4.pr.752x$AcceptDelt)

library(ggplot2)
library(tidyverse)
PostSum <- PosteriorSum(SumRes4.pr.752ax, T0 = 50, dfm = 7, Data = Data,lwd = 2)

PostSum <- PosteriorSum(SumRes4.pr.752x, T0 = 50, dfm = 7, Burn = NULL, probs = 0.95,Data = Data,lwd = 2)






SumRes4.pr.7520x <- BayesFIVSingV2(Nsim = 2000,
                                   Y = c(Data$Y), W = Data$Data$W, M = Data$Data$M, Z = Data$Z , PsiAccept = 1, Ge=2,
                                   sig2.propdelt = 1.5, sig2.priordelt = 2.0)

mean(SumRes4.pr.7520x$AcceptDelt)


SumRes4.pr.752x <- BayesFIVSing(Nsim = 2000,
                                Y = c(Data$Y), W = Data$Data$W, M = Data$Data$M, Z = Data$Z , PsiAccept = 1,
                                sig2.propdelt = 1.75, sig2.priordelt = 2.0, sig02.del = 0.75)#sig2.propdelt = 3.75, sig2.priordelt = 5.0

mean(SumRes4.pr.752x$AcceptDelt)

SumRes4.pr.752ax <- BayesFIVMult(Nsim = 10000,
                                Y = c(Data$Y), W = Data$Data$W, M = Data$Data$M, Z = Data$Z , PsiAccept = 1,
                                sig2.propdelt = 1.75, sig2.priordelt = 2.0, sig02.del = 0.75)#sig2.propdelt = 3.75, sig2.priordelt = 5.0

mean(SumRes4.pr.752ax$AcceptDelt)


#--- Look at posterior estimates

PostSum <- PosteriorSum(SumRes4.pr.752x, T0 = 50, dfm = 10)


#--- PLot the posterior estimate of beta(t)
par(mfrow=c(1,1))
matplot(a, t(PostSum$GamSm), type= "l", col="lightgrey")
lines(a, colMeans(PostSum$GamSm), col="red")
lines(a, Data$Betat, col="green")

#--- Plot posterior estimates of delta(t)
#--- Plot the delta(t) 
matplot(a, t(PostSum$DeltatSm), type= "l", col="lightgrey")
lines(a, colMeans(PostSum$DeltatSm), col="red",lwd = 4)
lines(a, Data$Delt, col="green",lwd=4)
