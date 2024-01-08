library(randomForestSRC)
library(party)
library(pec)
library(rotsf)
library(aorsf)
library(tictoc)
library(rockchalk) #for the noise function

#----------------------------------------------------#
#                  Noise function                    #
#----------------------------------------------------#
create.noise <- function(n){
  mu <- rep(0, 30)
  stddev <- rep(1, 30)
  corMat <- lazyCor(0.5, 30) #first number here is the correlation coefficient
  covMat <- stddev %*% t(stddev) * corMat
  dat1 <- pnorm(mvrnorm(n = n, mu = mu, Sigma = covMat, empirical = FALSE))
  
  mu2 <- rep(0, 20)
  stddev2 <- rep(1, 20)
  corMat2 <- lazyCor(0, 20)
  covMat2 <- stddev2 %*% t(stddev2) * corMat2
  dat2 <- pnorm(mvrnorm(n = n, mu = mu2, Sigma = covMat2, empirical = FALSE))
  
  noise.dat <- as.data.frame(cbind(dat1, dat2))
}
#-----------------------------------------------------#
#                  RUN SIMULATION                     #
#-----------------------------------------------------#

#First set parameters, initialize empty vectors, and read in various accompanying R scripts

flc.no <- 48 #48 is 1/3 of all FLCs. I ran the code below separately for the FLCs with 1) No ExtraVars; 2) Noisy variables; and 3) Confounders
reps <- 10
Num <- rep(seq(1:flc.no),each=reps)

tau = c(0.025,0.1,0.5,0.8,0.9,0.975) #quantiles to estimate

C.index.1 <- cbind(matrix(NA, reps*flc.no, length(tau)), Num)
C.index.6 <- C.index.5 <- C.index.4 <- C.index.3 <- C.index.2 <- C.index.1
q_loss.1 <- cbind(matrix(NA, reps*flc.no, length(tau)), Num)
q_loss.6 <- q_loss.5 <- q_loss.4 <- q_loss.3 <- q_loss.2 <- q_loss.1
IBS.1 <- cbind(matrix(NA, reps*flc.no, 1), Num)
IBS.6 <- IBS.5 <- IBS.4 <- IBS.3 <- IBS.2 <- IBS.1 
Yquant0.025.1 <- cbind(matrix(NA, reps*flc.no, 200), Num)
Yquant0.025.6 <- Yquant0.025.5 <- Yquant0.025.4 <- Yquant0.025.3 <- Yquant0.025.2 <- Yquant0.025.1 
Yquant0.10.1 <- cbind(matrix(NA, reps*flc.no, 200), Num)
Yquant0.10.6 <- Yquant0.10.5 <- Yquant0.10.4 <- Yquant0.10.3 <- Yquant0.10.2 <- Yquant0.10.1
Yquant0.50.1 <- cbind(matrix(NA, reps*flc.no, 200), Num)
Yquant0.50.6 <- Yquant0.50.5 <- Yquant0.50.4 <- Yquant0.50.3 <- Yquant0.50.2 <- Yquant0.50.1
Yquant0.80.1 <- cbind(matrix(NA, reps*flc.no, 200), Num)
Yquant0.80.6 <- Yquant0.80.5 <- Yquant0.80.4 <- Yquant0.80.3 <- Yquant0.80.2 <- Yquant0.80.1
Yquant0.90.1 <- cbind(matrix(NA, reps*flc.no, 200), Num)
Yquant0.90.6 <- Yquant0.90.5 <- Yquant0.90.4 <- Yquant0.90.3 <- Yquant0.90.2 <- Yquant0.90.1
Yquant0.975.1 <- cbind(matrix(NA, reps*flc.no, 200), Num)
Yquant0.975.6 <- Yquant0.975.5 <- Yquant0.975.4 <- Yquant0.975.3 <- Yquant0.975.2 <- Yquant0.975.1

#Create vectors for the actual max quantiles (AMQ) available for each observation, for all 6 methods
AMQ.1 <- cbind(matrix(NA, reps*flc.no, 200), Num)
AMQ.6 <- AMQ.5 <- AMQ.4 <- AMQ.3 <- AMQ.2 <- AMQ.1
source("C_index_function.r")
source("IBS_wt.r")
source("l1survival_func.r")
source("RISTfunctions.r")

source("help_functions.r")
source("help_functions_l1.r")
source("help_functions_rotsf.r")
source("help_functions_orsf.r")
source("help_functions_RIST.r")

#RIST parameters
K = 7     # number of covariate considered per spilt, usually sqrt(p) or p/3
nmin = 6  # minimum number of observed data in each node, default is 6.
M = 50    # number of trees in each fold, default is 50
L = 1 #i.e., RIST_1 -- number of folds, 1-5 are recommended.

maxtime <- rep(NA, reps*flc.no)
for(l in 1:(reps*flc.no)){
  maxtime[l] <- max(dat.list[[l]]$time)
  dat.list[[l]]$status = (dat.list[[l]]$status==1)
}
tao = 1.2*max(maxtime) # length of study. Must be larger than all survival times.
#As long as it is larger than all survival times, the value of it will not make any difference.


tic()
for(k in 1:(reps*flc.no)){
  n=length(dat.list[[k]]$time)
  n_train=length(dat.list[[k]]$time)-200
  train <- dat.list[[k]][1:n_train,-c(1,3)]
  test.error <- dat.list[[k]][(n_train+1):n,]
  test <- test.error[,-c(1,3)]
  
  ##-----add in noise - comment out when not relevant-----#
  #noise.dat <- create.noise(n)
  #train <- cbind(train, noise.dat[1:n_train,])
  #test <- cbind(test, noise.dat[(n_train+1):n,])
  #-------------------------------------------------------#
  
  ##-----add in confounders - comment out when not relevant-----#
  #length1 <- length(train)
  #sig.sq <- 0.25
  #conf1 <- rnorm(n, 0, sig.sq) + c(train$X.V2, test$X.V2) #continuous var 1
  #conf2 <- rnorm(n, 0, sig.sq) + c(train$X.V4, test$X.V4) #continuous var 2
  #conf3 <- rnorm(n, 0, sig.sq) + c(train$X.V5, test$X.V5) #continuous var 3
  #train <- cbind(train, conf1[1:n_train], conf2[1:n_train], conf3[1:n_train])
  #test <- cbind(test, conf1[(n_train+1):n], conf2[(n_train+1):n], conf3[(n_train+1):n])
  #colnames(train)[(length1+1):length(train)] <- c('conf1','conf2','conf3')
  #colnames(test)[(length1+1):length(test)] <- c('conf1','conf2','conf3')
  #-------------------------------------------------------#
  
  #Method 1: RSF-log-rank
  
  v.obj <- rfsrc(Surv(time, status) ~ ., data=train, ntree = 100, mtry=7, splitrule="logrank")
  rsf.pred <- predict(v.obj, test)
  
  #Pull out Y's from K-M curve in order to compute QL
  Yquant0.025.1[k,1:200] <- find_quantile(surv = rsf.pred, max_value = max(train$time), tau = tau[1])[[1]]
  Yquant0.10.1[k,1:200] <- find_quantile(surv = rsf.pred, max_value = max(train$time), tau = tau[2])[[1]]
  Yquant0.50.1[k,1:200] <- find_quantile(surv = rsf.pred, max_value = max(train$time), tau = tau[3])[[1]]
  Yquant0.80.1[k,1:200] <- find_quantile(surv = rsf.pred, max_value = max(train$time), tau = tau[4])[[1]]
  Yquant0.90.1[k,1:200] <- find_quantile(surv = rsf.pred, max_value = max(train$time), tau = tau[5])[[1]]
  Yquant0.975.1[k,1:200] <- find_quantile(surv = rsf.pred, max_value = max(train$time), tau = tau[6])[[1]]
  AMQ.1[k,1:200] <- find_quantile(surv = rsf.pred, max_value = max(train$time), tau = tau[1])[[2]]
  q_loss.1[k,1] <- mean((test.error$Tlat - Yquant0.025.1[k,1:200])*(tau[1] - 1*(test.error$Tlat - Yquant0.025.1[k,1:200] < 0)), na.rm = TRUE) #2.5th quantile
  q_loss.1[k,2] <- mean((test.error$Tlat - Yquant0.10.1[k,1:200])*(tau[2] - 1*(test.error$Tlat - Yquant0.10.1[k,1:200] < 0)), na.rm = TRUE) #10th quantile
  q_loss.1[k,3] <- mean((test.error$Tlat - Yquant0.50.1[k,1:200])*(tau[3] - 1*(test.error$Tlat - Yquant0.50.1[k,1:200] < 0)), na.rm = TRUE) #50th quantile
  q_loss.1[k,4] <- mean((test.error$Tlat - Yquant0.80.1[k,1:200])*(tau[4] - 1*(test.error$Tlat - Yquant0.80.1[k,1:200] < 0)), na.rm = TRUE) #80th quantile
  q_loss.1[k,5] <- mean((test.error$Tlat - Yquant0.90.1[k,1:200])*(tau[5] - 1*(test.error$Tlat - Yquant0.90.1[k,1:200] < 0)), na.rm = TRUE) #90th quantile
  q_loss.1[k,6] <- mean((test.error$Tlat - Yquant0.975.1[k,1:200])*(tau[6] - 1*(test.error$Tlat - Yquant0.975.1[k,1:200] < 0)), na.rm = TRUE) #97.5th quantile
  
  #C-index on 2.5th,10th,50th (median),90th,97.5th quantiles
  C.index.1[k,1] <- c.index(test$time, test$status, Yquant0.025.1[k,1:200])
  C.index.1[k,2] <- c.index(test$time, test$status, Yquant0.10.1[k,1:200])
  C.index.1[k,3] <- c.index(test$time, test$status, Yquant0.50.1[k,1:200])
  C.index.1[k,4] <- c.index(test$time, test$status, Yquant0.80.1[k,1:200])
  C.index.1[k,5] <- c.index(test$time, test$status, Yquant0.90.1[k,1:200])
  C.index.1[k,6] <- c.index(test$time, test$status, Yquant0.975.1[k,1:200])
  
  trngY <- cbind(train$time, train$status)
  testY <- cbind(test$time, test$status)
  
  IBS.1[k,1] <- df.aibs.fxn.wt(rsf.pred$survival, trngY, testY, RSF = TRUE)
  
  #Method 2: L1
  v.obj.l1 <- rfsrc.l1(Surv(time, status) ~ ., train, test, ntree = 100)
  
  surv.prob <- matrix(NA, 200, 1500)
  surv.time <- matrix(NA, 200, 1500)
  surv.cens <- matrix(NA, 200, 1500)
  
  for (i in 1:nrow(test)){
    shatl1=survfit(Surv(time, status) ~ 1,type="kaplan-meier",data=train[v.obj.l1[[i]][1,],],weights=v.obj.l1[[i]][2,])
    surv.censor <- shatl1$n.event
    surv.censor[surv.censor >= 1] <- 1 
    shatl1=cbind(shatl1$time,shatl1$surv)
    surv.cens[i,1:nrow(shatl1)] <- surv.censor
    surv.prob[i,1:nrow(shatl1)] <- shatl1[,2]
    surv.time[i,1:nrow(shatl1)] <- shatl1[,1]
  }
  
  #Pull out Y's from K-M curve in order to compute QL
  Yquant0.025.2[k,1:200] <- find_quantile.l1(surv1 = surv.prob, surv2 = surv.time, tau = tau[1])[[1]]
  Yquant0.10.2[k,1:200] <- find_quantile.l1(surv1 = surv.prob, surv2 = surv.time, tau = tau[2])[[1]]
  Yquant0.50.2[k,1:200] <- find_quantile.l1(surv1 = surv.prob, surv2 = surv.time, tau = tau[3])[[1]]
  Yquant0.80.2[k,1:200] <- find_quantile.l1(surv1 = surv.prob, surv2 = surv.time, tau = tau[4])[[1]]
  Yquant0.90.2[k,1:200] <- find_quantile.l1(surv1 = surv.prob, surv2 = surv.time, tau = tau[5])[[1]]
  Yquant0.975.2[k,1:200] <- find_quantile.l1(surv1 = surv.prob, surv2 = surv.time, tau = tau[6])[[1]]
  AMQ.2[k,1:200] <- find_quantile.l1(surv1 = surv.prob, surv2 = surv.time, tau = tau[1])[[2]]
  q_loss.2[k,1] <- mean((test.error$Tlat - Yquant0.025.2[k,1:200])*(tau[1] - 1*(test.error$Tlat - Yquant0.025.2[k,1:200] < 0)), na.rm = TRUE) #2.5th quantile
  q_loss.2[k,2] <- mean((test.error$Tlat - Yquant0.10.2[k,1:200])*(tau[2] - 1*(test.error$Tlat - Yquant0.10.2[k,1:200] < 0)), na.rm = TRUE) #10th quantile
  q_loss.2[k,3] <- mean((test.error$Tlat - Yquant0.50.2[k,1:200])*(tau[3] - 1*(test.error$Tlat - Yquant0.50.2[k,1:200] < 0)), na.rm = TRUE) #50th quantile
  q_loss.2[k,4] <- mean((test.error$Tlat - Yquant0.80.2[k,1:200])*(tau[4] - 1*(test.error$Tlat - Yquant0.80.2[k,1:200] < 0)), na.rm = TRUE) #80th quantile
  q_loss.2[k,5] <- mean((test.error$Tlat - Yquant0.90.2[k,1:200])*(tau[5] - 1*(test.error$Tlat - Yquant0.90.2[k,1:200] < 0)), na.rm = TRUE) #90th quantile
  q_loss.2[k,6] <- mean((test.error$Tlat - Yquant0.975.2[k,1:200])*(tau[6] - 1*(test.error$Tlat - Yquant0.975.2[k,1:200] < 0)), na.rm = TRUE) #97.5th quantile
  
  #C-index on 2.5th,10th,50th (median),90th,97.5th quantiles
  C.index.2[k,1] <- c.index(test$time, test$status, Yquant0.025.2[k,1:200])
  C.index.2[k,2] <- c.index(test$time, test$status, Yquant0.10.2[k,1:200])
  C.index.2[k,3] <- c.index(test$time, test$status, Yquant0.50.2[k,1:200])
  C.index.2[k,4] <- c.index(test$time, test$status, Yquant0.80.2[k,1:200])
  C.index.2[k,5] <- c.index(test$time, test$status, Yquant0.90.2[k,1:200])
  C.index.2[k,6] <- c.index(test$time, test$status, Yquant0.975.2[k,1:200])
  
  trngY <- cbind(train$time, train$status)
  testY <- cbind(test$time, test$status)
  
  #Compute alternative IBS for L1 -- we ultimately discarded this, hence why it's commented out
  #IBS.2[k,1] <- df.aibs.fxn.wt(rsf.pred$survival, trngY, testY, RSF = TRUE)
  
  
  #Method 3: CIF
  v.obj <- pecCforest(Surv(time, status) ~ ., data=train, controls = cforest_unbiased())
  
  #Now use 'cforest' object to predict at uncensored survival times from training data
  train.death <- data.frame(death.times = train[ , 1], death.indx = train[ , 2])
  train.dt <- train.death[train.death$death.indx == TRUE, ]
  unique.times <- unique(train.dt$death.times)
  
  rsf.pred <- predictSurvProb(v.obj,newdata=test[,3:length(test)],times=unique.times)
  
  #Pull out quantiles / compute QL
  Yquant0.025.3[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[1])[[1]]
  Yquant0.10.3[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[2])[[1]]
  Yquant0.50.3[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[3])[[1]]
  Yquant0.80.3[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[4])[[1]]
  Yquant0.90.3[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[5])[[1]]
  Yquant0.975.3[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[6])[[1]]
  AMQ.3[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[1])[[2]]
  q_loss.3[k,1] <- mean((test.error$Tlat - Yquant0.025.3[k,1:200])*(tau[1] - 1*(test.error$Tlat - Yquant0.025.3[k,1:200] < 0)), na.rm = TRUE) #2.5th quantile
  q_loss.3[k,2] <- mean((test.error$Tlat - Yquant0.10.3[k,1:200])*(tau[2] - 1*(test.error$Tlat - Yquant0.10.3[k,1:200] < 0)), na.rm = TRUE) #10th quantile
  q_loss.3[k,3] <- mean((test.error$Tlat - Yquant0.50.3[k,1:200])*(tau[3] - 1*(test.error$Tlat - Yquant0.50.3[k,1:200] < 0)), na.rm = TRUE) #50th quantile
  q_loss.3[k,4] <- mean((test.error$Tlat - Yquant0.80.3[k,1:200])*(tau[4] - 1*(test.error$Tlat - Yquant0.80.3[k,1:200] < 0)), na.rm = TRUE) #80th quantile
  q_loss.3[k,5] <- mean((test.error$Tlat - Yquant0.90.3[k,1:200])*(tau[5] - 1*(test.error$Tlat - Yquant0.90.3[k,1:200] < 0)), na.rm = TRUE) #90th quantile
  q_loss.3[k,6] <- mean((test.error$Tlat - Yquant0.975.3[k,1:200])*(tau[6] - 1*(test.error$Tlat - Yquant0.975.3[k,1:200] < 0)), na.rm = TRUE) #97.5th quantile
  
  #C-index on 2.5th,10th,50th (median),90th,97.5th quantiles
  C.index.3[k,1] <- c.index(test$time, test$status, Yquant0.025.3[k,1:200])
  C.index.3[k,2] <- c.index(test$time, test$status, Yquant0.10.3[k,1:200])
  C.index.3[k,3] <- c.index(test$time, test$status, Yquant0.50.3[k,1:200])
  C.index.3[k,4] <- c.index(test$time, test$status, Yquant0.80.3[k,1:200])
  C.index.3[k,5] <- c.index(test$time, test$status, Yquant0.90.3[k,1:200])
  C.index.3[k,6] <- c.index(test$time, test$status, Yquant0.975.3[k,1:200])
  
  trngY <- cbind(train$time, train$status)
  testY <- cbind(test$time, test$status)
  
  IBS.3[k,1] <- df.aibs.fxn.wt(rsf.pred, trngY, testY, RSF = TRUE)

  
  #Method 4: RotSF
  v.obj <- rotsfpca(Surv(time, status) ~ ., data=train, trlength = 500, m = 2,
                    na.action = na.omit, vari_status = FALSE)
  
  #Now use 'rotsf' object to predict at uncensored survival times from training data
  train.death <- data.frame(death.times = train[ , 1], death.indx = train[ , 2])
  train.dt <- train.death[train.death$death.indx == TRUE, ]
  unique.times <- unique(train.dt$death.times)
  rsf.pred <- rotsfpca.surv_predict(v.obj, test[,3:length(test)], uniquetimes = unique.times, trlength = 500)
  
  #Pull out quantiles / compute QL
  Yquant0.025.4[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[1])[[1]]
  Yquant0.10.4[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[2])[[1]]
  Yquant0.50.4[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[3])[[1]]
  Yquant0.80.4[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[4])[[1]]
  Yquant0.90.4[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[5])[[1]]
  Yquant0.975.4[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[6])[[1]]
  AMQ.4[k,1:200] <- find_quantile.cif_rotsf(surv = rsf.pred, max_value = max(train$time), tau = tau[1])[[2]]
  q_loss.4[k,1] <- mean((test.error$Tlat - Yquant0.025.4[k,1:200])*(tau[1] - 1*(test.error$Tlat - Yquant0.025.4[k,1:200] < 0)), na.rm = TRUE) #2.5th quantile
  q_loss.4[k,2] <- mean((test.error$Tlat - Yquant0.10.4[k,1:200])*(tau[2] - 1*(test.error$Tlat - Yquant0.10.4[k,1:200] < 0)), na.rm = TRUE) #10th quantile
  q_loss.4[k,3] <- mean((test.error$Tlat - Yquant0.50.4[k,1:200])*(tau[3] - 1*(test.error$Tlat - Yquant0.50.4[k,1:200] < 0)), na.rm = TRUE) #50th quantile
  q_loss.4[k,4] <- mean((test.error$Tlat - Yquant0.80.4[k,1:200])*(tau[4] - 1*(test.error$Tlat - Yquant0.80.4[k,1:200] < 0)), na.rm = TRUE) #80th quantile
  q_loss.4[k,5] <- mean((test.error$Tlat - Yquant0.90.4[k,1:200])*(tau[5] - 1*(test.error$Tlat - Yquant0.90.4[k,1:200] < 0)), na.rm = TRUE) #90th quantile
  q_loss.4[k,6] <- mean((test.error$Tlat - Yquant0.975.4[k,1:200])*(tau[6] - 1*(test.error$Tlat - Yquant0.975.4[k,1:200] < 0)), na.rm = TRUE) #97.5th quantile
  
  #C-index on 2.5th,10th,50th (median),90th,97.5th quantiles
  C.index.4[k,1] <- c.index(test$time, test$status, Yquant0.025.4[k,1:200])
  C.index.4[k,2] <- c.index(test$time, test$status, Yquant0.10.4[k,1:200])
  C.index.4[k,3] <- c.index(test$time, test$status, Yquant0.50.4[k,1:200])
  C.index.4[k,4] <- c.index(test$time, test$status, Yquant0.80.4[k,1:200])
  C.index.4[k,5] <- c.index(test$time, test$status, Yquant0.90.4[k,1:200])
  C.index.4[k,6] <- c.index(test$time, test$status, Yquant0.975.4[k,1:200])
  
  trngY <- cbind(train$time, train$status)
  testY <- cbind(test$time, test$status)
  
  IBS.4[k,1] <- df.aibs.fxn.wt(rsf.pred, trngY, testY, RSF = TRUE)
  
  
  #Method 5: ORSF
  orsf.mod=orsf(Surv(time, status) ~ ., data=train,n_tree=100)
  times = sort(train$time[train$status==1])
  rsf.pred <- predict(orsf.mod, new_data=test, pred_type = 'surv', pred_horizon=times)
  
  #Pull out quantiles / compute QL
  Yquant0.025.5[k,1:200] <- find_quantile.orsf(surv = rsf.pred, max_value = max(train$time), tau = tau[1])[[1]]
  Yquant0.10.5[k,1:200] <- find_quantile.orsf(surv = rsf.pred, max_value = max(train$time), tau = tau[2])[[1]]
  Yquant0.50.5[k,1:200] <- find_quantile.orsf(surv = rsf.pred, max_value = max(train$time), tau = tau[3])[[1]]
  Yquant0.80.5[k,1:200] <- find_quantile.orsf(surv = rsf.pred, max_value = max(train$time), tau = tau[4])[[1]]
  Yquant0.90.5[k,1:200] <- find_quantile.orsf(surv = rsf.pred, max_value = max(train$time), tau = tau[5])[[1]]
  Yquant0.975.5[k,1:200] <- find_quantile.orsf(surv = rsf.pred, max_value = max(train$time), tau = tau[6])[[1]]
  AMQ.5[k,1:200] <- find_quantile.orsf(surv = rsf.pred, max_value = max(train$time), tau = tau[1])[[2]]
  q_loss.5[k,1] <- mean((test.error$Tlat - Yquant0.025.5[k,1:200])*(tau[1] - 1*(test.error$Tlat - Yquant0.025.5[k,1:200] < 0)), na.rm = TRUE) #2.5th quantile
  q_loss.5[k,2] <- mean((test.error$Tlat - Yquant0.10.5[k,1:200])*(tau[2] - 1*(test.error$Tlat - Yquant0.10.5[k,1:200] < 0)), na.rm = TRUE) #10th quantile
  q_loss.5[k,3] <- mean((test.error$Tlat - Yquant0.50.5[k,1:200])*(tau[3] - 1*(test.error$Tlat - Yquant0.50.5[k,1:200] < 0)), na.rm = TRUE) #50th quantile
  q_loss.5[k,4] <- mean((test.error$Tlat - Yquant0.80.5[k,1:200])*(tau[4] - 1*(test.error$Tlat - Yquant0.80.5[k,1:200] < 0)), na.rm = TRUE) #80th quantile
  q_loss.5[k,5] <- mean((test.error$Tlat - Yquant0.90.5[k,1:200])*(tau[5] - 1*(test.error$Tlat - Yquant0.90.5[k,1:200] < 0)), na.rm = TRUE) #90th quantile
  q_loss.5[k,6] <- mean((test.error$Tlat - Yquant0.975.5[k,1:200])*(tau[6] - 1*(test.error$Tlat - Yquant0.975.5[k,1:200] < 0)), na.rm = TRUE) #97.5th quantile
  
  #C-index on 2.5th,10th,50th (median),90th,97.5th quantiles
  C.index.5[k,1] <- c.index(test$time, test$status, Yquant0.025.5[k,1:200])
  C.index.5[k,2] <- c.index(test$time, test$status, Yquant0.10.5[k,1:200])
  C.index.5[k,3] <- c.index(test$time, test$status, Yquant0.50.5[k,1:200])
  C.index.5[k,4] <- c.index(test$time, test$status, Yquant0.80.5[k,1:200])
  C.index.5[k,5] <- c.index(test$time, test$status, Yquant0.90.5[k,1:200])
  C.index.5[k,6] <- c.index(test$time, test$status, Yquant0.975.5[k,1:200])
  
  trngY <- cbind(train$time, train$status)
  testY <- cbind(test$time, test$status)
  
  IBS.5[k,1] <- df.aibs.fxn.wt(rsf.pred, trngY, testY, RSF = TRUE)
  
  
  #Method 6: RIST-1
  train <- dat.list[[k]][1:n_train,]
  train <- train[,c(5:26,4,2)]
  test.error <- dat.list[[k]][(n_train+1):n,]
  test.error <- test.error[,c(1,5:26,4,2,3)]
  test <- dat.list[[k]][(n_train+1):n,]
  test <- test[,c(5:26,4,2)]
  
  R_Muti_ERT_build = Muti_ERT_fit(train, M, K, L, nmin, SupLogRank=1, tao, impute="random")
  R_Muti_ERT_predict= Muti_ERT_Predict(test, R_Muti_ERT_build$Forest_seq[[L]], R_Muti_ERT_build$SurvMat_seq[[L]], R_Muti_ERT_build$time_intrest)
  
  Yquant0.025.6[k,1:200] <- find_quantile.rist(surv = R_Muti_ERT_predict, tau = tau[1])[[1]]
  Yquant0.10.6[k,1:200] <- find_quantile.rist(surv = R_Muti_ERT_predict, tau = tau[2])[[1]]
  Yquant0.50.6[k,1:200] <- find_quantile.rist(surv = R_Muti_ERT_predict, tau = tau[3])[[1]]
  Yquant0.80.6[k,1:200] <- find_quantile.rist(surv = R_Muti_ERT_predict, tau = tau[4])[[1]]
  Yquant0.90.6[k,1:200] <- find_quantile.rist(surv = R_Muti_ERT_predict, tau = tau[5])[[1]]
  Yquant0.975.6[k,1:200] <- find_quantile.rist(surv = R_Muti_ERT_predict, tau = tau[6])[[1]]
  AMQ.6[k,1:200] <- find_quantile.rist(surv = R_Muti_ERT_predict, tau = tau[1])[[2]]
  q_loss.6[k,1] <- mean((test.error$Tlat - Yquant0.025.6[k,1:200])*(tau[1] - 1*(test.error$Tlat - Yquant0.025.6[k,1:200] < 0)), na.rm = TRUE) #2.5th quantile
  q_loss.6[k,2] <- mean((test.error$Tlat - Yquant0.10.6[k,1:200])*(tau[2] - 1*(test.error$Tlat - Yquant0.10.6[k,1:200] < 0)), na.rm = TRUE) #10th quantile
  q_loss.6[k,3] <- mean((test.error$Tlat - Yquant0.50.6[k,1:200])*(tau[3] - 1*(test.error$Tlat - Yquant0.50.6[k,1:200] < 0)), na.rm = TRUE) #50th quantile
  q_loss.6[k,4] <- mean((test.error$Tlat - Yquant0.80.6[k,1:200])*(tau[4] - 1*(test.error$Tlat - Yquant0.80.6[k,1:200] < 0)), na.rm = TRUE) #80th quantile
  q_loss.6[k,5] <- mean((test.error$Tlat - Yquant0.90.6[k,1:200])*(tau[5] - 1*(test.error$Tlat - Yquant0.90.6[k,1:200] < 0)), na.rm = TRUE) #90th quantile
  q_loss.6[k,6] <- mean((test.error$Tlat - Yquant0.975.6[k,1:200])*(tau[6] - 1*(test.error$Tlat - Yquant0.975.6[k,1:200] < 0)), na.rm = TRUE) #97.5th quantile
  
  #C-index on 2.5th,10th,50th (median),90th,97.5th quantiles
  C.index.6[k,1] <- c.index(test$time, test$status, Yquant0.025.6[k,1:200])
  C.index.6[k,2] <- c.index(test$time, test$status, Yquant0.10.6[k,1:200])
  C.index.6[k,3] <- c.index(test$time, test$status, Yquant0.50.6[k,1:200])
  C.index.6[k,4] <- c.index(test$time, test$status, Yquant0.80.6[k,1:200])
  C.index.6[k,5] <- c.index(test$time, test$status, Yquant0.90.6[k,1:200])
  C.index.6[k,6] <- c.index(test$time, test$status, Yquant0.975.6[k,1:200])
  
  trngY <- cbind(train$time, train$status)
  testY <- cbind(test$time, test$status)
  
  IBS.6[k,1] <- df.aibs.fxn.wt(R_Muti_ERT_predict[[2]], trngY, testY, RSF = TRUE)
  
  print(paste0("Finished data set #", k))
  
}
toc()

#Assemble results for analysis
Full_results1 <- cbind(C.index.1, C.index.2, C.index.3,
                       C.index.4, C.index.5, C.index.6,
                       IBS.1, IBS.2, IBS.3, IBS.4, IBS.5, IBS.6,
                       q_loss.1, q_loss.2, q_loss.3,
                       q_loss.4, q_loss.5, q_loss.6)

#Do the same for the other two thirds of the FLCs, then merge in Excel. The result is the .csv file called "Full_results_merged".
