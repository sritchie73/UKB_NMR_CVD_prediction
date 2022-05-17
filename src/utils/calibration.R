library(survival)
library(data.table)
library(foreach)
source("src/utils/aki_absrisk.R") # extract predicted absolute risk from Cox models

calibration.fit <- function(formula, data, time, predicted_risk, byage=FALSE) {
  formula <- gsub(" +", " ", formula)
  surv_term <- gsub(" ?~.*", "", formula) 
  follow_col <- gsub("Surv\\(", "", gsub(",.*", "", surv_term))
  event_col <- gsub(".*, ?", "", gsub("\\)", "", surv_term))

  depend_vars <- strsplit(gsub(".*~ ?", "", formula), " ?\\+ ?")[[1]]
  strata_terms <- depend_vars[depend_vars %like% "^strata\\("]
  strata_vars <- gsub("strata\\(", "", gsub("\\)", "", strata_terms))

  # Drop censored individuals
  non_censored <- copy(data)
  setnames(non_censored, c(event_col, follow_col), c("event", "follow"))
  non_censored <- non_censored[event == 1 | follow >= time]
  setnames(non_censored, c("event", "follow"), c(event_col, follow_col))
  data <- non_censored

  # Get predicted absolute risk
  if(missing(predicted_risk)) {
		# Fit cox model
		cph <- coxph(as.formula(formula), data=data, x=TRUE)

		abs_risk <- Coxar(cph, time)

		# Build table of strata terms and absolute risk
		dt <- data[, ..strata_vars]
		dt[, rn := .I]
		dt <- dt[as.integer(rownames(cph$x))]
		dt[, predicted_risk := abs_risk]
  } else {
    dt <- data[, ..strata_vars]
    dt[, rn := .I]
    dt[, predicted_risk := predicted_risk] 
  }

  if (!byage) {
		# Bin risk into deciles
		dt <- dt[order(predicted_risk),.SD,by=strata_vars]
		dt[, risk_decile := floor((seq_len(.N)-1)/(.N/10))+1, by=strata_vars]
  } else {
    dt[data[,.(rn=.I, age)], on = .(rn), age := age]
    dt[, risk_decile := floor(age/5)*5] 
  }

  # Get strata groups
  strata_groups <- data[,.N,by=strata_vars]

  # Fit kaplan meier survival estimates in each risk decile in each strata
  obs_risk <- foreach(sg = strata_groups[,.I], .combine=rbind) %do% {
    foreach(rd = dt[,unique(risk_decile)], .combine=rbind) %do% {
      rd_dat <- data[dt[risk_decile == rd, rn]]
      rd_dat <- rd_dat[strata_groups[sg], on = c(strata_vars)]
			km <- survfit(as.formula(sprintf("%s ~ 1", surv_term)), data=rd_dat)
			kms <- summary(km, time=time)
      cbind(strata_groups[sg], data.table(risk_decile=rd, obs.risk = 1 - kms$surv, obs.L95 = 1 - kms$upper, obs.U95 = 1 - kms$lower))
    }
  }

  # Summarise predicted risk per decile
  pred_risk <- dt[, .(
      pred.risk=mean(predicted_risk), 
      pred.L95=quantile(predicted_risk, prob=0.025), 
      pred.U95=quantile(predicted_risk, prob=0.975)), 
    by=c(strata_vars, "risk_decile")]

  # Merge
  comp <- pred_risk[obs_risk, on = c(strata_vars, "risk_decile")]

  # If by age group, add
  if (byage) {
    setnames(comp, "risk_decile", "age_group")
    comp[, age_group := sprintf("%s-%s", age_group, age_group+4)]
  }
  
  # Return
  return(comp)
} 
