library(survival)
library(flexsurv)

#Simulate survival data
simulSurv <- function(N, lambda, rho, beta, rateC, sigma, DGM.Y, DGM.X, a, b, k.0, k.1, rateC.0, rateC.1)
{
  #main binary covariate denoting group membership
  x1 <- sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5)) #control group: "0"; #treatment: "1"
  
  #additional covariates
  num.coef = 5 #number of additional *continuous* regression coefficients
  X <- cbind(x1, matrix(rep(runif(N*num.coef, 0, 1)),ncol=num.coef))
  
  #5-level *categorical* covariates
  x7 <- as.factor(sample(x=c(0, 1, 2, 3, 4), size=N, replace=TRUE, prob=c(0.05, 0.1, 0.2, 0.2, 0.45)))
  x7.1 <- ifelse(x7=="1", 1, 0)
  x7.2 <- ifelse(x7=="2", 1, 0)
  x7.3 <- ifelse(x7=="3", 1, 0)
  x7.4 <- ifelse(x7=="4", 1, 0)
  
  x8 <- as.factor(sample(x=c(0, 1, 2, 3, 4), size=N, replace=TRUE, prob=c(0.05, 0.1, 0.2, 0.2, 0.45)))
  x8.1 <- ifelse(x8=="1", 1, 0)
  x8.2 <- ifelse(x8=="2", 1, 0)
  x8.3 <- ifelse(x8=="3", 1, 0)
  x8.4 <- ifelse(x8=="4", 1, 0)
  
  x9 <- as.factor(sample(x=c(0, 1, 2, 3, 4), size=N, replace=TRUE, prob=c(0.05, 0.1, 0.2, 0.2, 0.45)))
  x9.1 <- ifelse(x9=="1", 1, 0)
  x9.2 <- ifelse(x9=="2", 1, 0)
  x9.3 <- ifelse(x9=="3", 1, 0)
  x9.4 <- ifelse(x9=="4", 1, 0)
  
  x10 <- as.factor(sample(x=c(0, 1, 2, 3, 4), size=N, replace=TRUE, prob=c(0.05, 0.1, 0.2, 0.2, 0.45)))
  x10.1 <- ifelse(x10=="1", 1, 0)
  x10.2 <- ifelse(x10=="2", 1, 0)
  x10.3 <- ifelse(x10=="3", 1, 0)
  x10.4 <- ifelse(x10=="4", 1, 0)
  
  X <- cbind(X, x7.1, x7.2, x7.3, x7.4, x8.1, x8.2, x8.3, x8.4, 
             x9.1, x9.2, x9.3, x9.4, x10.1, x10.2, x10.3, x10.4)
  
  
  #Weibull, AFT or GG latent event times
  IndY.1 <- ifelse(DGM.Y=="1", 1, 0)
  IndY.2 <- ifelse(DGM.Y=="2", 1, 0)
  IndY.3 <- ifelse(DGM.Y=="3", 1, 0)
  
  #Linear, Poly, Sine, or Piecewise X's
  IndX.1 <- ifelse(DGM.X=="1", 1, 0)
  IndX.2 <- ifelse(DGM.X=="2", 1, 0)
  IndX.3 <- ifelse(DGM.X=="3", 1, 0)
  IndX.4 <- ifelse(DGM.X=="4", 1, 0)
  
  Linear <- X%*%beta[1:22]
    
  QuadInt <- X[,1]*beta[1] + X[,2]^2*beta[2] + X[,3]*beta[3] + 
    X[,4]^2*beta[4] + X[,3]*X[,4]*beta[5] + X[,4]^2*X[,5]*beta[6] +
    X[,6]^3*beta[7] + X[,7:22]%*%beta[8:23]
  
  Sine <- X[,1]*beta[1] + sin(X[,2])*beta[2] + X[,3]*beta[3] + 
    X[,4]*beta[4] + X[,5]*beta[5] + X[,6]*beta[6] +
    X[,7:22]%*%beta[7:22]
    
  Piecewise <- X[,1]*beta[1] + (1*X[,2]<0.3)*beta[2] + 
    (1*X[,2]>=0.3 & 1*X[,2]<0.8)*beta[3] + (1*X[,2]>=0.8)*beta[4] +
    (1*X[,3]<0.3)*beta[5] + (1*X[,3]>=0.3 & 1*X[,3]<0.8)*beta[6] +
    (1*X[,3]>=0.8)*beta[7] + X[,4]*beta[8] + X[,5]*beta[9] + X[,6]*beta[10] +
    X[,7:22]%*%beta[11:26]
  
  LP <- ifelse(DGM.X=="1", list(Linear), 
          ifelse(DGM.X=="2", list(QuadInt), 
            ifelse(DGM.X=="3", list(Sine), list(Piecewise))))
  LP <- LP[[1]]
  
  v <- runif(N)
  
  Weib <-  (- log(v) / (lambda * exp(LP)))^(1 / rho) 
  AFT <- exp(LP + rnorm(N, 0, sigma))
  GG <-  rep(NA, N)
  group0 <- which(X[,1]==0)
  group1 <- which(X[,1]==1)
  #C <- IndY.3*( rep(NA, N) )
  GG[group0] <- IndY.3*( rgengamma.orig(length(group0), shape=b, scale=a*exp(X[group0,]%*%beta[1:ncol(X)]), k=k.0) )
  #Generalized Gamma (Weibull) for group 1
  GG[group1] <- IndY.3*( rgengamma.orig(length(group1), shape=b, scale=a*exp(X[group1,]%*%beta[1:ncol(X)]), k=k.1) )
  
  Tlat <- ifelse(DGM.Y=="1", list(Weib), 
            ifelse(DGM.Y=="2", list(AFT), list(GG)))
  Tlat <- Tlat[[1]]
  
  #generate censoring times with Exponential distributed data
  C.WA <- rexp(n=N, rate=rateC)
  
  C.GG <- rep(NA, N)
  C.GG[group0] <- IndY.3*( rexp(n=length(group0), rate=rateC.0) )
  C.GG[group1] <- IndY.3*( rexp(n=length(group1), rate=rateC.1) )
  
  C <- ifelse(DGM.Y=="3", list(C.GG), list(C.WA))
  C <- C[[1]]
  
  # follow-up times and event indicators
  time <- pmin(Tlat, C)
  status <- as.numeric(Tlat <= C)
  
  # data set
  dat <- data.frame(id=1:N,
                    time=time,
                    Tlat=Tlat,
                    status=status,
                    X=X)
}

#Settings for data simulation
set.seed(1600)
reps=10
flc.no=48 #number of FLCs
dat.list <- vector(mode = "list", length = reps*flc.no)

lambda=0.8 #Weibull scale
rho=1.8 #Weibull shape
sigma=0.52 #AFT SD
k.0=1/16 #GG param
k.1=1 #GG param
a=1.5 #GG param
b=2 #GG param
beta1.1=c(-5, 4, -2.5, -4, -5, 0, -3, 3.5, -1.5, 3.5, -2.5, -3, 4.5, -2, 0, -5.5, 0, -2.5, 0, 0, 0, 0) #Weibull, Linear
beta1.2=c(-5, 4, -2.5, -3, 2, -4, 2, -3, 3.5, -1.5, 3.5, -2.5, -3, 4.5, -2, 0, -5.5, 0, -2.5, 0, 0, 0, 0) #Weibull, Poly
beta1.3=c(-5, 4, -2.5, -4, -5, 0, -3, 3.5, -1.5, 3.5, -2.5, -3, 4.5, -2, 0, -5.5, 0, -2.5, 0, 0, 0, 0) #Weibull, Sine
beta1.4=1.6*c(-5, 4, 2, -1.5, -2.5, -1, 2,  -4, -5, 0, -3, 3.5, -1.5, 3.5, -2.5, -3, 4.5, -2, 0, -5.5, 0, -2.5, 0, 0, 0, 0) #Weibull, Piecewise
beta2.1=c(-0.5, 2, -1.5, 1.75, -1.5, 1.4, 1.0, 0.5, -0.75, 0.5, -1.5, -2, 2.5, -1, 1, 0.5, 0, -1.5, 0, 0, 0, 0) #AFT, Linear
beta2.2=c(-0.5, 2, -1.5, 1.5, -1 -1.5, 2, 1.4, 1.0, 0.5, -0.75, 0.5, -1.5, -2, 2.5, -1, 1, 0.5, 0, -1.5, 0, 0, 0, 0) #AFT, Poly
beta2.3=c(-0.5, 2, -1.5, 1.75, -1.5, 1.4, 1.0, 0.5, -0.75, 0.5, -1.5, -2, 2.5, -1, 1, 0.5, 0, -1.5, 0, 0, 0, 0) #AFT, Sine
beta2.4=c(-0.5, 2, 1, 0.5, -1.5, -0.5, 0.5, 1.75, -1.5, 1.4, 1.0, 0.5, -0.75, 0.5, -1.5, -2, 2.5, -1, 1, 0.5, 0, -1.5, 0, 0, 0, 0) #AFT, Piecewise
beta3.1=2.3*c(-2, 0.5, 0.2, 0.7, 0.5, 0, 0.7, 0.8, 0.5, 0.3, 0.4, 0.4, 0.5, 0.8, 0, 0.5, 0, 0.5, 0, 0, 0, 0) #GG, Linear
beta3.2=2.1*c(-2, 0.5, 0.2, 0.6, 0.4, 0.3, 0.2, 0.7, 0.8, 0.5, 0.3, 0.4, 0.4, 0.5, 0.8, 0, 0.5, 0, 0.5, 0, 0, 0, 0) #GG, Poly
beta3.3=2.4*c(-2, 0.5, 0.2, 0.7, 0.5, 0, 0.7, 0.8, 0.5, 0.3, 0.4, 0.4, 0.5, 0.8, 0, 0.5, 0, 0.5, 0, 0, 0, 0) #GG, Sine
beta3.4=4.3*c(-2, 1.7, 0.5, 0.2, 0.9, 0.1, -0.1, 0.7, 0.5, 0, 0.7, 0.8, 0.5, 0.3, 0.4, 0.4, 0.5, 0.8, 0, 0.5, 0, 0.5, 0, 0, 0, 0) #GG, Piecewise

#Loop through rows of this dataframe
flc <- read.csv("flc.csv", na.strings=c("","NA"))

for(row in 1:nrow(flc)){
  for(k in 1:reps){
  dat.list[[row*reps-10+k]] <- simulSurv(N=flc$N[row], rateC=flc$rateC[row], beta=get(flc$beta[row]), DGM.X=flc$DGM.X[row], DGM.Y=flc$DGM.Y[row], #THIS CHANGES - 1
                             lambda=lambda, rho=rho, sigma=sigma, #lambda/rho only for Weibull; sigma only for AFT
                             a=a, b=b, k.0=k.0, k.1=k.1, #only for GG
                             rateC.0=flc$rateC.0[row], rateC.1=flc$rateC.1[row]) #only for GG 
  }
}

