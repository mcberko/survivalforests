
#############  Function to build a L1 forest with survival data

rfsrc.l1=function(formula,datrain,datest,ntrees=100){


### run survival random forests with default or custom (L1 here) splitting rule and extract the BOPs for test data
# formula = model formula
# datrain = training data frame
# datest = test data frame (for which we want predictions)
# ntrees = number of trees in the forest

# output:
# List with the BOPs (Bag of Observations for Prediction). Element i in the list is the BOP for the ith test data, in matrix form.
#     	row 1 = indices of the observation in datrain
#	row 2 = weight of the corresponding observation (number of times the observation appears in the BOP) 


ntrain=dim(datrain)[1]
ntest=dim(datest)[1]

boptest=vector("list",length=ntest)

for (i in 1:ntrees)
{
  ind=sample(1:ntrain,ntrain,replace=TRUE)  	# indices of the bootstrap sample
  dat1=datrain[ind,]      			# bootstrap sample to build the tree

		# The L1 splitting rule is in the custom2 slot
   fit=rfsrc(formula,data=dat1,ntree=1,bootstrap="none",sampsize=n,splitrule="custom2",membership=TRUE) 

   fitmem=rep(NA,ntrain)
   fitmem[ind]=fit$membership			# get the terminal node IDs for the training data

   pptest=predict(fit,newdata=datest,membership=TRUE)
   testmem=pptest$membership			# get the node terminal node IDs for the test data

	for (i  in 1:ntest)
	{
	boptest[[i]]=c(boptest[[i]],which(fitmem==testmem[i]))	# find which training observations are in the same terminal node as a test observation
	}
}

	for(i in 1:ntest)	
	{
	tab=table(boptest[[i]])
	boptest[[i]]=rbind(as.numeric(names(tab)),as.vector(tab))	# convert the BOPs into matrix form to save space
	}

return(boptest)

}


# Example

library(randomForestSRC)
library(survival)

ntrain=100      # training size
ntest=5	# test size
n=ntrain+ntest

# generate data
x1=runif(n)
x2=runif(n)
time=x1+x2+rexp(n)
# 20% censoring
status=rbinom(n,1,.8)

dat=data.frame(x1,x2,time,status)
datrain=dat[1:ntrain,]
datest=dat[(ntrain+1):n,]

###### off-the-shelf forest
fitlogrank=rfsrc(Surv(time, status) ~ x1+x2,data=datrain,ntree=100)
predlogrank=predict(fitlogrank,newdata=datest)
# estimated survival function for the ith test observation
i=1
shatlogrank=cbind(predlogrank$time.interest,predlogrank$survival[i,])


###### L1 forest
fitl1=rfsrc.l1(Surv(time, status) ~ x1+x2,datrain,datest,ntrees=100)
# The function returns the BOPs
fitl1

# With these BOPs, you can compute any summary quantity you want.
# For example, here is how to compute the estimated survival function for the ith test observation
i=1
shatl1=survfit(Surv(time, status) ~ 1,type="kaplan-meier",data=datrain[fitl1[[i]][1,],],weights=fitl1[[i]][2,])
shatl1=cbind(shatl1$time,shatl1$surv)





