library(rms)
library(data.table)
source("src/utils/aki_absrisk.R")

#' Recalibrate 10-year risk of CVD
#' 
#' Recalibrates each individual's predicted 10-year risk of CVD based on 
#' 10-year incidence rates in CPRD. See supplement of Sun L, *et al*. 
#' Use of polygenic risk scores and other molecular markers to enhance 
#' cardiovascular risk prediction: prospective cohort study and modelling 
#' analysis. *bioRxiv* (2019). doi: 10.1101/744565. for details.
#' 
#' @param cph Object to use to predict each individuals risk. May be 
#'            a vector giving each person's risk of CVD in the 10-year
#'            follow-up period, or may be an object returned by the 
#'            `cph` function in the `rms` package from which 
#'            we will predict each person's 10-year risk. If providing
#'            a `cph` object, the model must be fit with time in study 
#'            as time-scale for a 10-year study follow-up, adjusting for
#'            age as a covariate. The model should also be stratified by
#'            sex. Additional covariates are also permitted 
#'            in the model.
#' @param IID The column from the original data the `cph` model was fit
#'            on giving the sample's identifier.
#' @param age The column from the original data the `cph` model was fit 
#'            on containing each person's age.
#' @param sex The column from the original data the `cph` model was fit
#'            on indicating whether each person is male or female. By
#'            default, assumes 0 is female, and 1 is male (use the `male` 
#'            argument to change this).
#' @param male value in the `sex` column assigned to males.
#' 
#' @return a data.table/data.frame containing three columns: 
#'    (1) each individual's sample identifier (IID).
#'    (2) each individual's predicted 10-year risk based on the input
#'        Cox Proportional Hazards modle (cph).
#'    (3) each individual's 10-year CVD risk re-calibrated to the annual
#'        10-year incidence rates in CPRD.
#'        
#' @section **Adapting this function**
#' To adapt this to other endpoint definitions or follow-up times, 
#' incidence rates in CPRD need to be calculated based on the desired 
#' endpoint, in age-group buckets of size follow-up-time / 2, for the
#' range of baseline ages + follow-up-time / 2 in the dataset of interest.
#' 
#' E.g., say we wanted to recalibrated 6-year risk of hard CAD in UK 
#' Biobank rather than 10-year risk of CVD, we would need to obtain
#' incidence rates of hard CAD in 3 year age-groups (6-year follow-up 
#' time divided by two) for people aged 39 to 77 (Taking age-groups as 
#' 36-38, 39-41, ..., 72-74, because UK Biobank pariticipants are aged 
#' 37-73, then we'd need corresponding incidence rates in CPRD age-groups 
#' 39-41, 42-44, ..., 75-77 to match the desired 6-year follow-up time).
#'        
recalibrate_risk = function(cph, IID, age, sex, male="M") {
  # Construct table from input columns, this will also be used to hold
  # predicted risk, recalibrated risk, and various intermediate columns
  # for doing the recalibration.
  dt = data.table(IID, age, sex=as.character(sex))
  dt[sex != male, sex := "female"]
  dt[sex == male, sex := "male"]
  
  # Age- and sex-specific incidence rates (per 1000 person-years) of CVD 
  # among CPRD participants (n=3,117,544) from the Appendix table of 
  # Sun L, *et al*. Use of polygenic risk scores and other molecular 
  # markers to enhance cardiovascular risk prediction: prospective 
  # cohort study and modelling analysis. *bioRxiv* (2019). 
  # doi: 10.1101/744565.
  CPRD = data.table(
    age_at_risk_start = seq(40, 75, by=5),
    age_at_risk_end = seq(44, 79, by=5),
    male_rate_per_1000py = c(1.578, 2.964, 5.141, 7.561, 10.848, 13.731, 19.271, 25.667),
    female_rate_per_1000py = c(0.726, 1.309, 2.165, 2.969, 4.560, 7.007, 10.840, 17.007)
  )
  
  ##
  ## Step 1. Estimate the 10-year CVD risk in each person from the Cox 
  ## Proportional Hazards model provided by the user.
  ##

  if (inherits(cph, "coxph") && !inherits(cph, "rms")) {
    stop("Error: Cox model must be fit with the rms package not the survival package")
  } else if (inherits(cph, "rms")) {
    if (cph$maxtime != 10) {
      warning("Estimating absolute individual risk over maximum follow-up time of ", cph$maxtime, 
              " years, then recalibrating to 10-year risk")
    }
    # Estimate each individual's risk based on their dependent variables
    prob_surv <- survest(cph, linear.predictors = cph$linear.pred, time = cph$maxtime, conf.int = FALSE)
    individual_risk <- 1 - prob_surv$surv
    
    # Handle the fact that these predictions may not always be made in all
    # individuals (predictions can only be made in complete cases in the 
    # data, missing values in any terms/covariates lead to that person
    # not being present in the model residuals).
    names(individual_risk) = dt$IID[as.numeric(names(individual_risk))] 
    individual_risk = data.table(IID = as.integer(names(individual_risk)), predicted_risk = individual_risk)
    dt[individual_risk, on = .(IID), predicted_risk := i.predicted_risk]
  } else {
    dt = cbind(dt, predicted_risk=cph) # user provided individual risk
  }

  ##
  ## Step 2. Calculate the mean of the predicted 10-year CVD risk in
  ## each five year age group in UK Biobank, stratified by sex.
  ##
  
  # Assign each person to a 5-year age group
  dt[, age_group_start := floor(age/5)*5]
  dt[, age_group_end := floor(age/5)*5 + 4]

  # Determine average risk of developing CVD in the 10-year follow-up
  # for each age group, stratified by sex.
  mean_risk = dt[, .(mean_predicted_risk=mean(predicted_risk, na.rm=TRUE)), 
                  by = .(sex, age_group_start, age_group_end)]
  mean_risk = dcast(mean_risk, age_group_start + age_group_end ~ sex, 
                    value.var="mean_predicted_risk")
  setnames(mean_risk, c("male", "female"), c("male_predicted_risk", "female_predicted_risk"))
  
  ##
  ## Step 3. Calculate the expected 10-year CVD risk for each age-group
  ## based on the incidence rates in CPRD. Expected 10-year risk for each
  ## 5-year age group is based on the incidence rates at the mid-point of
  ## the interval ahead. E.g. for the 40 to 44 year age-group the 
  ## incidence rate for 45 to 49 years was used.
  ##
  
  # Calculate annual CVD incidence in each age-at-risk group from the 
  # 1000 person-year CVD rates:
  CPRD[, male_annual_incidence := male_rate_per_1000py / 1000]
  CPRD[, female_annual_incidence := female_rate_per_1000py / 1000]
  
  # The expected 10-year CVD risk is calculated for each age group
  # based on the mid-point of the next interval ahead, e.g. for the 
  # 40-44 year age-group the expected 10-year risk is calculated based
  # on the annual incidence rates in the 45-49 year olds in CPRD:
  CPRD[, age_group_start := age_at_risk_start - 5]
  CPRD[, age_group_end := age_at_risk_end - 5]

  # The expected risk is calculated assuming exponential survival (i.e.
  # constant hazard) in the 10-years ahead.
  CPRD[, male_expected_risk := 1 - exp(-male_annual_incidence * 10)]
  CPRD[, female_expected_risk := 1 - exp(-female_annual_incidence * 10)]
  
  ##
  ## Step 4: Fit the sex-stratified recalibration model. 
  ##
  
  # helper functions
  link_function = function(...) { log(-log(...)) }
  inverse_link_function = function(...) { exp(-exp(...)) }
  model_alpha = function(...) { coef(...)[1] }
  model_beta = function(...) { coef(...)[2] }
  
  # Aggregate risk in a single table
  aggregate_risk = CPRD[mean_risk, on = .(age_group_start, age_group_end)]
  
  # Fit the recalibration model for males as described in equation 2
  male_recalibration_model = aggregate_risk[,
    lm(link_function(1 - male_expected_risk) ~ link_function(1 - male_predicted_risk))
  ]
  
  # Fit the recalibration model for females as described in equation 2
  female_recalibration_model = aggregate_risk[,
    lm(link_function(1 - female_expected_risk) ~ link_function(1 - female_predicted_risk))
  ]
  
  ##
  ## Step 5: Use the recalibration models to recalibrate each individuals 
  ## predicted 10-year CVD risk.
  ##
  
  # Collate model coefficients into a single table
  recalibration_coefficients = data.table(
    sex=c("male", "female"), 
    alpha=c(model_alpha(male_recalibration_model), model_alpha(female_recalibration_model)),
    beta=c(model_beta(male_recalibration_model), model_beta(female_recalibration_model))
  )
  
  # Recalibrate each individuals predicted risk using the intercept and 
  # beta from the model fit above for each sex as described in equation 3.
  dt[recalibration_coefficients, on = .(sex), 
      recalibrated_risk := 1 - inverse_link_function(alpha + beta * link_function(1 - predicted_risk))]
  
  # Remove recalibrated predicted risk for individuals outside the age 
  # range for which the recalibration is valid
  dt[age < min(CPRD$age_group_start) | age > max(CPRD$age_group_end), recalibrated_risk := NA]

  # Add risk predicted by age and sex stratified CPRD model
  CPRD = melt(CPRD, id.vars=c("age_group_start", "age_group_end"),
              measure.vars=c("male_expected_risk", "female_expected_risk"),
              variable.name="sex", value.name="expected_risk")
  CPRD[, sex := gsub("_expected_risk", "", sex)]
  dt[CPRD, on = .(age >= age_group_start, age <= age_group_end, sex),
      CPRD_risk := expected_risk]
  
  # Return predicted and recalibrated risk for each person.
  return(dt[, .(IID, predicted_risk, recalibrated_risk, CPRD_risk)])
}
