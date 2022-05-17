require(survival)

################################################################################
# Aki Havulinna's code

getStrata <- function(coxfit, newdata=NULL, warn=TRUE) {
   # this function extracts strata for each observation in the coxfit (or model frame) object
   # safer to use "newdata"
   temp <- untangle.specials(terms(coxfit), "strata")
   if(length(temp$vars) < 1) {
      if(warn) warning('no strata in the model')
      return(NULL)
   }
   if(!is.null(newdata)) {
      f <- model.frame(coxfit, data=newdata)
   } else mf <- model.frame(coxfit)
   if (length(temp$vars) == 1) stratum <- mf[[temp$vars]]
   else stratum <- strata(mf[, temp$vars], shortlabel = TRUE)
   stratum
}

findhaz <- function(bh, time0, time1=NULL, stratum=NULL)
{
   if(is.null(bh$strata)) { 
      if(!is.null(stratum)) stop('did not expect strata')
      af <- approxfun(bh$time, bh$hazard, rule=2)
      h <- -af(time0)
      if(!is.null(time1)) h <- -h - af(time1)
   } else {
      if(is.null(stratum)) {
	 stop('basehaz has strata, expected strata as an argument')
      }
      h <- rep(NA,length(time0))
      us <- unique(bh$strata) 
      if(length(setdiff(us, unique(stratum)))) {
	 stop('stratum has categories which are not in baseline')
      }
      for(st in us) {
	 idx <- bh$strata == st
	 idx2 <- stratum == st
	 af <- approxfun(bh$time[idx], bh$hazard[idx], rule=2)
	 h[idx2] <- -af(time0[idx2]) # this is h(end) if time1 is not given
	 if(!is.null(time1)) {
	    h[idx2] <- -h[idx2] - af(time1[idx2])
	 }
      }
   }
   h
}

absrisk <- function(surv, bh, lp=0, stratum=NULL, years=10)
{
   # calculates absolute risk given survival indicators,
   # baseline hazard and linear predictor + possible strata
   n <- nrow(surv)
   type <- attr(surv, 'type')
   if(type == 'counting') {
      haz <- findhaz(bh, surv[,1], surv[,1] + years, stratum=stratum)
   } else if (type == 'right') {
      haz <- findhaz(bh, rep(years, n), stratum=stratum)
   } else {
      stop('Unsupported censoring type: ', type) 
   }
   1 - exp(haz * exp(lp))  
}

################################################################################
# Absolute risk from a Cox model
Coxar <- function(mod, years=10, bh=NULL)
{
   if(is.null(bh)) {
      bh <- basehaz(mod, centered=TRUE)
   }
   stratum <- getStrata(mod, warn=FALSE)
   return(absrisk(surv=mod$y, bh=bh, lp=mod$linear.predictor,
      stratum=stratum, years=years))
}


