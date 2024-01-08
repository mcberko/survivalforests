#This functions takes 3 arguments: observed time, censoring status, and the predicted median (or other quantile)
c.index <- function(Y, status, predmedian){
  delta=c()
  omit=c()
  ind=1
  for(i in 1:(n-n_train-1)) {
    for(j in (i+1):(n-n_train)) {
      
      if(Y[i]<Y[j]) { 
        Ymin=i
        Ymax=j 
      } else {
        Ymin=j
        Ymax=i
      }  
      
      #Compare predicted quantiles with empirical quantiles
      if(!is.na(predmedian[Ymin]) & !is.na(predmedian[Ymax])) {
        
        if(Y[Ymin]<Y[Ymax] & status[Ymin]==1 & predmedian[Ymin]<predmedian[Ymax]) { delta[ind]=1
        } else if(Y[Ymin]<Y[Ymax] & status[Ymin]==1 & predmedian[Ymin]==predmedian[Ymax]) { delta[ind]=0.5
        } else if(Y[Ymin]==Y[Ymax] & status[Ymin]==1 & status[Ymax]==1 & predmedian[Ymin]==predmedian[Ymax]) { delta[ind]=1
        } else if(Y[Ymin]==Y[Ymax] & status[Ymin]==1 & status[Ymax]==1 & predmedian[Ymin]!=predmedian[Ymax]) { delta[ind]=0.5
        } else delta[ind]=0
        
      } else delta[ind]=NA
      
      #Omit non-permissible pairs
      if(status[Ymin]==0){omit[ind]=1
      } else if(Y[Ymin]==Y[Ymax] & status[Ymin]==0 & status[Ymax]==0){omit[ind]=1
      } else omit[ind]=0
      
      ind=ind+1
    }
  }
  total=sum(!is.na(delta))-sum(omit)
  Cindex=sum(delta,na.rm=T)/total
  return(Cindex)
}

