# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#

quantile_loss <- function(u, tau) {
  return (u*(tau - (u < 0)))
}

# help functions
get.event.info <- function(obj, subset = NULL) {
  ## survival case
  if (grepl("surv", obj$family)) {
    if (!is.null(obj$yvar)) {
      if (is.null(subset)) {
        subset <- (1:nrow(cbind(obj$yvar)))
      }
      r.dim <- 2
      time <- obj$yvar[subset, 1]
      cens <- obj$yvar[subset, 2]
      ## censoring must be coded coherently
      if (!all(floor(cens) == abs(cens), na.rm = TRUE)) {
        stop("for survival families censoring variable must be coded as a non-negative integer")
      }
      ## Extract the unique event types.
      event <- na.omit(cens)[na.omit(cens) > 0]
      event.type <- sort(unique(event))
    }
    ##everything else
    else {
      r.dim <- 0
      event <- event.type <- cens <- cens <- time <- NULL
    }
    ## Set grid of time points.
    time.interest <- obj$time.interest
  }
  else {
    ## NULL for other families
    if ((obj$family == "regr+") | (obj$family == "class+")) {
      r.dim <- dim(obj$yvar)[2]
    }
    else {
      r.dim <- 1
    }
    event <- event.type <- cens <- time.interest <- cens <- time <- NULL
  }
  return(list(event = event, event.type = event.type, cens = cens,
              time.interest = time.interest, time = time, r.dim = r.dim))
}

find_quantile <- function(surv, max_value, tau){
  b_size <- dim(surv$survival)[1]
  event.info <- get.event.info(surv)
  bin <- event.info$time.interest
  quantiles <- rep(0, b_size)
  actual.max.quant <- rep(0, b_size)
  for (i in 1:b_size){
    check <- surv$survival[i,] >= 1-tau
    j <- max(sum(check), 1)
    quantile <- bin[j]
    quantiles[i] <- quantile
    actual.max.quant[i] <- 1-surv$survival[i,sum(!is.na(bin))]
  }
  
  return(list(v1=quantiles,v2= actual.max.quant))
}

#Matt:
#surv=rsf.pred
#tau=0.9
