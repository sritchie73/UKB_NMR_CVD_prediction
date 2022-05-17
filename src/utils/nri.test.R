library(data.table)
library(survival)
library(nricens)

source("src/utils/aki_absrisk.R") # extract predicted absolute risk from Cox models

# Tabulate net reclassification index:
# - continuous
# - Categorical: <5%, 5-10%, >10% (NICE 2014 guidelines)
# - Categorical: <5%, 5-7.5%, >7.5% (ACC/AHA 2019 guidelines)
nri.test <- function(base_model, new_model, time_horizon, dat, event_col, follow_col, bootstraps) {
  # Fit cox models for both models to compare
  if (inherits(base_model, "formula")) {
    dat <- as.data.frame(dat)
		base_cph <- coxph(base_model, data=dat, x=TRUE, y=TRUE)
		new_cph <- coxph(new_model, data=dat, x=TRUE, y=TRUE)

		base_cph_risk <- Coxar(base_cph, time_horizon)
		new_cph_risk <- Coxar(new_cph, time_horizon) 

		event_vec <- dat[,event_col]
		time_vec <- dat[,follow_col]
  } else {
    base_cph_risk <- base_model
    new_cph_risk <- new_model

    event_vec <- event_col
    time_vec <- follow_col
  }

  if (length(event_vec) != length(base_cph_risk) && 
      length(event_vec) != length(new_cph_risk)) {
    stop("Cannot work with missing data in either model's variables")
  }

  # Calculate continuous NRI
  contNRI <- nricens(event = event_vec, time = time_vec, p.std = base_cph_risk, p.new = new_cph_risk,
                     updown = "diff", cut = 0, t0 = time_horizon, niter=bootstraps)

  # Calculate NRI for NICE 2014
  niceNRI <- nricens(event = event_vec, time = time_vec, p.std = base_cph_risk, p.new = new_cph_risk,
                    updown = "category", cut = c(0.05, 0.1), t0 = time_horizon, niter=bootstraps)

  # Calculate NRI for AHA/ACC 2019
  ahaaccNRI <- nricens(event = event_vec, time = time_vec, p.std = base_cph_risk, p.new = new_cph_risk,
                      updown = "category", cut = c(0.05, 0.075), t0 = time_horizon, niter=bootstraps)

  list("Continuous NRI"=contNRI, "NICE 2014 Categorical NRI"=niceNRI, "ACC/AHA 2019 Categorical NRI"=ahaaccNRI)
}

