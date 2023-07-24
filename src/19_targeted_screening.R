library(data.table)
library(foreach)
library(survival)
library(boot)

# create output directory
system("mkdir -p analyses/public_health_modelling/targeted_screening")

# Load hypothetical population 
ons_pop <- fread("analyses/public_health_modelling/ONS_hypothetical_100k_pop_by_age_sex.txt")

# Load in predicted risks
pred_risk <- fread("analyses/CVD_score_weighting/CVD_linear_predictors_and_risk.txt")
pred_risk[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# For targeted screening, we will follow up on those at intermediate risk according to SCORE2
score2 <- pred_risk[model == "SCORE2", .(eid, SCORE2_risk=uk_risk)]
pred_risk <- pred_risk[model != "SCORE2", .(eid, sex, age_group, incident_cvd_followup, incident_cvd, model, pred_risk=uk_risk)]
pred_risk <- score2[pred_risk, on = .(eid)]

# Proportions allocated to high risk group and CVD cases therein estimate in bootstrap as with blanket screening
risk_strata <- foreach(this_model = unique(pred_risk$model), .combine=rbind) %:% 
  foreach(this_sex = c("Males", "Females"), .combine=rbind) %:% 
    foreach(this_age_group = sort(unique(pred_risk$age_group)), .combine=rbind) %dopar% {
      # Extract subset of predicted risks to work with in this loop iteration
      this_pred_risk <- pred_risk[model == this_model & sex == gsub("s", "", this_sex) & age_group == this_age_group]

      # Function to pass to bootstrap to compute. 
      phs <- function(dt) {
        # ESC 2021 guidelines use different cut-offs depending on age, we'll use one step down for medium risk
        if (dt$age[1] >= 50) {
          high_risk_threshold <- 0.1    
          medium_risk_threshold <- 0.075
        } else {
          high_risk_threshold <- 0.075
          medium_risk_threshold <- 0.05
        }
      
        # Classification of people into high risk proceeds in two steps: 
        # (1) people meeting the high risk threshold with SCORE2 alone
        # (2) people at medium risk with the SCORE2 prediction are then re-screened using the new model
        high_risk <- dt[
          SCORE2_risk > high_risk_threshold | 
          (SCORE2_risk > medium_risk_threshold & pred_risk > high_risk_threshold)
        ]

        # Get proportion of samples allocated to high risk group and proportion of high risk group that are cases
        if (nrow(high_risk) == 0) {
          pct_high_risk <- 0
          pct_high_risk_cases <- 0
        } else {
          pct_high_risk <- high_risk[,.N]/dt[,.N]
          pct_high_risk_cases <- high_risk[,sum(incident_cvd)]/high_risk[,.N]
        }
        c(pct_high_risk, pct_high_risk_cases)
      }
      
      # Run bootstrap analysis
      surv_cols_idx <- match(c("incident_cvd_followup", "incident_cvd"), names(this_pred_risk))
      phs.boots <- censboot(this_pred_risk, phs, 1000, index=surv_cols_idx) # takes ~ 6 seconds to run
 
      # Extract 95% CIs
      pct_high_risk_95CI <- boot.ci(phs.boots, type="perc", index=1)
      pct_high_risk_cases_95CI <- boot.ci(phs.boots, type="perc", index=2)
    
      # Build table of results to return
      data.table(
        model = this_model, sex = this_sex, age_group = this_age_group, 
        high_risk = phs.boots$t0[1], high_risk_L95 = pct_high_risk_95CI$percent[1, 4], high_risk_U95 = pct_high_risk_95CI$percent[1, 5],
        high_risk_cases = phs.boots$t0[2], high_risk_cases_L95 = pct_high_risk_cases_95CI$percent[1, 4], high_risk_cases_U95 = pct_high_risk_cases_95CI$percent[1, 5]
      )
}
fwrite(risk_strata, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening/ukb_risk_strata_proportions.txt")

# Apply proportions to simulated population
pop_strata <- risk_strata[ons_pop, on = .(sex, age_group)]

pop_strata[, high_risk := floor(N * high_risk)]
pop_strata[, high_risk_L95 := floor(N * high_risk_L95)]
pop_strata[, high_risk_U95 := floor(N * high_risk_U95)]

pop_strata[, high_risk_cases := floor(high_risk * high_risk_cases)]
pop_strata[, high_risk_cases_L95 := floor(high_risk * high_risk_cases_L95)]
pop_strata[, high_risk_cases_U95 := floor(high_risk * high_risk_cases_U95)]

pop_strata[, high_risk_non_cases := high_risk - high_risk_cases]
pop_strata[, high_risk_non_cases_L95 := high_risk_L95 - high_risk_cases_L95]
pop_strata[, high_risk_non_cases_U95 := high_risk_U95 - high_risk_cases_U95]

pop_strata[, low_risk := N - high_risk]
pop_strata[, low_risk_L95 := N - high_risk_U95]
pop_strata[, low_risk_U95 := N - high_risk_L95]

pop_strata[, low_risk_cases := cases - high_risk_cases]
pop_strata[, low_risk_cases_L95 := cases - high_risk_cases_U95]
pop_strata[, low_risk_cases_U95 := cases - high_risk_cases_L95]

pop_strata[, low_risk_non_cases := low_risk - low_risk_cases]
pop_strata[, low_risk_non_cases_L95 := low_risk_L95 - low_risk_cases_L95]
pop_strata[, low_risk_non_cases_U95 := low_risk_U95 - low_risk_cases_U95]

pop_strata <- pop_strata[, .(
  sex, age_group, people=N, cases, non_cases=controls, model,
  high_risk, high_risk_L95, high_risk_U95, low_risk, low_risk_L95, low_risk_U95,
  high_risk_cases, high_risk_cases_L95, high_risk_cases_U95, low_risk_cases, low_risk_cases_L95, low_risk_cases_U95,
  high_risk_non_cases, high_risk_non_cases_L95, high_risk_non_cases_U95, low_risk_non_cases, low_risk_non_cases_L95, low_risk_non_cases_U95
)]
fwrite(pop_strata, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening/ons_pop_stratified.txt")

# Apply targeted screening strategies to general population.
pop_screened <- melt(pop_strata, id.vars=c("sex", "age_group", "model"))
pop_screened <- pop_screened[, .(value=sum(value)), by=.(sex, model, variable)]
pop_screened <- dcast(pop_screened, sex + model ~ variable, value.var="value")

# Events prevented assuming 20% reduction due to statins
pop_screened[, events_prevented := floor(high_risk_cases * 0.2)]
pop_screened[, events_prevented_L95 := floor(high_risk_cases_L95 * 0.2)]
pop_screened[, events_prevented_U95 := floor(high_risk_cases_U95 * 0.2)]

# Number needed to screen per event prevented
pop_screened[, NNS := ceiling(people / events_prevented)]
pop_screened[, NNS_L95 := ceiling(people / events_prevented_L95)]
pop_screened[, NNS_U95 := ceiling(people / events_prevented_U95)]

# Number of statins prescribed per event prevented
pop_screened[, NNT := ceiling(high_risk / events_prevented)]
pop_screened[, NNT_L95 := ceiling(high_risk_L95 / events_prevented_L95)]
pop_screened[, NNT_U95 := ceiling(high_risk_U95 / events_prevented_U95)]

fwrite(pop_screened, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening/population_screening.txt")


