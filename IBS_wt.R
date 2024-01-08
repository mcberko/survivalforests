
## YUAN 2008 AIBS Function 
## -- adaptation for Matt
## Vinnie Liu
## November 29, 2021

## Source ----
library(dplyr)
library(survival)
library(randomForestSRC)


## Empirical Curve Fxn ----
##  parameters needed: 
## - death time (dt.i) and death index (di.i) for patient i from the test set
## - survival info for all patients from training set (trngY)
## - Whether mod/pred is RSF or not (T/F)
### ---> returns: empirical curve for patient i 

"emp.fxn" <- function(dt.i, di.i, trngY,  RSF){
  trng.death <- data.frame(death.times = trngY[ , 1], death.indx = trngY[ , 2])
  trng.dt <- trng.death[trng.death$death.indx == TRUE, ]
  t <- sort(unique(trng.dt$death.times))
 
   if(RSF == "TRUE"){
    if(di.i == TRUE){ ## FOR CENSORING
      I <- matrix(0, nrow = 2, ncol = length(t))
    }else{ I <- matrix(NA, nrow = 2, ncol = length(t))}
    
    I[1,] <- t
    
    if(dt.i < min(t)){
      I[2,] <- 0   }
    if(min(t) <= dt.i && dt.i <= max(t)){
      res <- max(which(t <= dt.i))
      I[2,(1:res)] <- 1   }
    if(dt.i > max(t)){
      I[2,] <- 1   }
  }
  
  if(RSF == "FALSE"){ ## TO TRUNCATE WEIB/COX FXNS
    trng.dtF <- trng.death
    trunc <- which(unique(trng.dtF$death.times) <= max(t))
    tmp <- trng.dtF$death.times[trunc]
    t2 <- sort(unique(trng.dtF$death.times)[trunc])
    
    I <-matrix(0, nrow = 2, ncol = length(t2))
    
    I[1,] <- t2
    
    if(dt.i < min(t2)){
      I[2,] <- 0   }
    if(min(t2) <= dt.i && dt.i <= max(t2)){
      res <- max(which(t2 <= dt.i))
      I[2,(1:res)] <- 1   }
    if(dt.i > max(t2)){
      I[2,] <- 1   }
  }
  return(I[2,])
}

## >>  Examples 
# patient 1
#tmp1 <- emp.fxn(dt.i= testY[1,1], di.i = testY[1,2], trngY, RSF = T)
## patient 84
#tmp2 <- emp.fxn(dt.i= testY[200,1], di.i = testY[200,2], trngY, RSF = T)

## EIBS Fxn ----
##  parameters needed: 
## - predicted survival curve (pred.i) and empirical curve (emp.i) for patient i 
##      from the test set
## - survival info for all patients from training set (trngY)
## - Whether mod/pred is RSF or not (T/F)
### --- > returns: IBS for patient i

"eibs.fxn.wt" <- function(pred.i, emp.i, trngY, RSF){
  trng.death <- data.frame(death.times = trngY[ , 1], death.indx = trngY[ , 2])
  trng.dt <- trng.death[trng.death$death.indx == TRUE, ]
  t <- sort(unique(trng.dt$death.times)) 
 
   if(RSF == "TRUE"){
    bs <- matrix(0, nrow = 2, ncol = length(t))
    bs[1,] <- t
    bs[2,] <- (pred.i - emp.i)^2
    
    reim <- c()
    reim[1] <- bs[2,1] * t[1]
    for(i in 2:length(t)){
      reim[i] <- bs[2,i] * (t[i] - t[i-1])/max(t)
    }
   } 
  
  if(RSF == "FALSE"){ ## TRUNCATE FOR WEIB/COX
    trng.dtF <- trng.death
    trunc <- which(unique(trng.dtF$death.times) <= max(t))
    tmp <- trng.dtF$death.times[trunc]
    t2 <- sort(unique(trng.dtF$death.times)[trunc])
    trunc.scurve <- pred.i[(1:length(t2))]
    
    bs <- matrix(0, nrow = 2, ncol = length(t2))
    bs[1,] <- t2
    bs[2,] <- (trunc.scurve - emp.i)^2
    
    reim <- c()
    reim[1] <- bs[2,1] * t2[1]
    for(i in 2:length(t)){
      reim[i] <- bs[2,i] * (t2[i] - t2[i-1])
    }
  } 
  return(sum(reim, na.rm=TRUE))
}

## >>  Examples 
#emp1 <- emp.fxn(dt.i= testY[24, 1], di.i = testY[24, 2], trngY, RSF = TRUE)
#eibs1 <- eibs.fxn(rf.pred$survival[24, ], emp1, trngY, RSF = TRUE)


## AIBS FXN ----
## parameters needed: 
## - predicted survival matrix (survmat)
## - whether the survmat is RSF or not (T/F)
## - survival info for all patients from training set (trngY) and test set (testY)
### ---> returns: AIBS for all patients in test set (testY)

"df.aibs.fxn.wt" <- function(survmat, trngY, testY, RSF){
  indiv.eibs <- c()
  for(i in 1:nrow(testY)){
    emp <- emp.fxn(dt.i = testY[i, 1], di.i = testY[i, 2], trngY, RSF = RSF)
    eibs.i <- eibs.fxn.wt(pred.i = survmat[i, ], emp.i = emp, trngY, RSF = RSF)
    indiv.eibs[i] <- eibs.i
  }
  aibs <- mean(indiv.eibs)
  return(aibs)
}

