library(data.table)
library(foreach)
library(survival)
library(boot)

# create output directory
system("mkdir -p analyses/public_health_modelling/targeted_screening")

# Load hypothetical population 
ons_pop <- fread("analyses/public_health_modelling/ONS_hypothetical_100k_pop_by_age_sex.txt")

# Load in predicted risks
pred_risk <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")
pred_risk <- pred_risk[score_type == "non-derived"]

# Reformat to wide table so all models have same bootstraps applied
models <- unique(pred_risk[,.(model)])
models <- models[model != "SCORE2"]
models[, colname := paste0("model", .I)]
pred_risk[, colname := model]
pred_risk[models, on = .(model), colname := i.colname]
pred_risk <- dcast(pred_risk, eid + sex + age_group + incident_cvd_followup + incident_cvd ~ colname, value.var="uk_calibrated_risk")

# Determine in each age-group and sex the proportion of participants determined to be high risk" 
# along with the proportion of "high risk" individuals that go on to develop CVD. These proportions 
# will then be applied to the simulated population to estimate impact of screening in a hypothetical 
# population of ~100,000 adults with proportions of age, sex, and 10-year CVD incidence similar to 
# the general UK Population. Proportions will be estimated in a bootstrap procedure so we can get 
# 95% confidence intervals of numbers.
risk_strata <- foreach(this_sex = c("Males", "Females"), .combine=rbind) %:% 
    foreach(this_age_group = sort(unique(pred_risk$age_group)), .combine=rbind) %dopar% {
      # Extract subset of predicted risks to work with in this loop iteration
      this_pred_risk <- pred_risk[sex == gsub("s", "", this_sex) & age_group == this_age_group]

      # Function to pass to bootstrap to compute. 
      boot.stats <- function(dt) {
        # ESC 2021 guidelines use different cut-offs depending on age
        if (dt$age[1] >= 50) {
          high_risk_threshold <- 0.1    
          medium_risk_threshold <- 0.075
        } else {
          high_risk_threshold <- 0.075
          medium_risk_threshold <- 0.05
        }
      
        # Put boostrapped input data into long format to simplify model comparison
        dt <- melt(dt, measure.vars=patterns("^model"), variable.name="model", value.name="uk_calibrated_risk")

        # Build empty results table to fill in
        res <- as.data.table(expand.grid(model=models$colname, metric=c(
          "pct_high_risk_score2", "pct_high_risk_score2_cases",
          "pct_medium_risk_score2", "pct_medium_risk_score2_cases",
          "pct_reclassified", "pct_reclassified_cases"
        )))
        res[, estimate := 0]

        # Compute % of group assigned to high risk and % of cases assigned to high risk for each model
        dt[, high_risk_score2 := ifelse(SCORE2 >= high_risk_threshold, TRUE, FALSE)]
        dt[, medium_risk_score2 := ifelse(!high_risk_score2 & SCORE2 >= medium_risk_threshold, TRUE, FALSE)]
        dt[, reclassified := ifelse(medium_risk_score2 & uk_calibrated_risk >= high_risk_threshold, TRUE, FALSE)]
    
        metric1 <- dt[, .(metric="pct_high_risk_score2", estimate=sum(high_risk_score2)/.N), by=model]
        res[metric1, on = .(model, metric), estimate := i.estimate]

        metric2 <- dt[(high_risk_score2), .(metric="pct_high_risk_score2_cases", estimate=sum(incident_cvd)/.N), by=model]
        res[metric2, on = .(model, metric), estimate := i.estimate]

        metric3 <- dt[, .(metric="pct_medium_risk_score2", estimate=sum(medium_risk_score2)/.N), by=model]
        res[metric3, on = .(model, metric), estimate := i.estimate]

        metric4 <- dt[(medium_risk_score2), .(metric="pct_medium_risk_score2_cases", estimate=sum(incident_cvd)/.N), by=model]
        res[metric4, on = .(model, metric), estimate := i.estimate]
       
        metric5 <- dt[(medium_risk_score2), .(metric="pct_reclassified", estimate=sum(reclassified)/.N), by=model]
        res[metric5, on = .(model, metric), estimate := i.estimate]

        metric6 <- dt[(reclassified), .(metric="pct_reclassified_cases", estimate=sum(incident_cvd)/.N), by=model]
        res[metric6, on = .(model, metric), estimate := i.estimate]

        # return as flat vector, this needs to be mapped after bootstrapping
        res[,estimate]
      }
      
      # Run bootstrap analysis
      surv_cols_idx <- match(c("incident_cvd_followup", "incident_cvd"), names(this_pred_risk))
      this_group_res <- censboot(this_pred_risk, boot.stats, 1000, index=surv_cols_idx) # takes ~ 6 seconds to run

      # Extract table of metrics
      res <- as.data.table(expand.grid(model=models$colname, metric=c(
        "pct_high_risk_score2", "pct_high_risk_score2_cases",
        "pct_medium_risk_score2", "pct_medium_risk_score2_cases",
        "pct_reclassified", "pct_reclassified_cases"
      )))
      res[models, on = .(model=colname), model := i.model]
      res[, estimate := this_group_res$t0]
      res <- cbind(sex=this_sex, age_group=this_age_group, res)

      # Extract table of bootstrap results
      boot_res <- as.data.table(this_group_res$t)
      boot_res[, bootstrap := .I]
      boot_res <- melt(boot_res, id.vars="bootstrap", value.name="estimate")
      boot_res <- boot_res[order(variable)][order(bootstrap)]
      boot_res <- boot_res[, cbind(res[,.(sex, age_group, model, metric)], .SD), by=bootstrap]
      boot_res[, variable := NULL]

      # return combined results
      rbind(res, boot_res, fill=TRUE)
}
risk_strata <- dcast(risk_strata, sex + age_group + model + bootstrap ~ metric, value.var="estimate")
risk_strata[is.na(bootstrap), bootstrap := 0] # simplify sorting and joining
risk_strata <- risk_strata[order(bootstrap)]

fwrite(risk_strata, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening/ukb_risk_strata_proportions_with_bootstraps.txt")

###########

# Apply proportions to simulated population
pop_boot <- risk_strata[ons_pop, on = .(sex, age_group)]

pop_boot[, high_risk_score2 := floor(N * pct_high_risk_score2)]
pop_boot[, medium_risk_score2 := floor(N * pct_medium_risk_score2)]
pop_boot[, low_risk_score2 := N - high_risk_score2 - medium_risk_score2]

pop_boot[, reclassified := floor(medium_risk_score2 * pct_reclassified)]
pop_boot[, not_reclassified := medium_risk_score2 - reclassified]

pop_boot[, high_risk := high_risk_score2 + reclassified]
pop_boot[, low_risk := low_risk_score2 + not_reclassified]

pop_boot[, high_risk_score2_cases := floor(high_risk_score2 * pct_high_risk_score2_cases)]
pop_boot[, high_risk_score2_non_cases := high_risk_score2 - high_risk_score2_cases]

pop_boot[, medium_risk_score2_cases := floor(medium_risk_score2 * pct_medium_risk_score2_cases)]
pop_boot[, medium_risk_score2_non_cases := medium_risk_score2 - medium_risk_score2_cases]

pop_boot[, low_risk_score2_cases := cases - high_risk_score2_cases - medium_risk_score2_cases]
pop_boot[, low_risk_score2_non_cases := low_risk_score2 - low_risk_score2_cases]

pop_boot[, reclassified_cases := floor(reclassified * pct_reclassified_cases)]
pop_boot[, reclassified_non_cases := reclassified - reclassified_cases]

pop_boot[, not_reclassified_cases := medium_risk_score2_cases - reclassified_cases]
pop_boot[, not_reclassified_non_cases := not_reclassified - not_reclassified_cases]

pop_boot[, high_risk_cases := high_risk_score2_cases + reclassified_cases]
pop_boot[, high_risk_non_cases := high_risk_score2_non_cases + reclassified_non_cases]

pop_boot[, low_risk_cases := low_risk_score2_cases + not_reclassified_cases]
pop_boot[, low_risk_non_cases := low_risk_score2_non_cases + not_reclassified_non_cases]


# Reorganise columns
pop_boot <- pop_boot[,.(
  sex, age_group, people=N, cases, non_cases=controls, model, bootstrap,
  high_risk, low_risk, high_risk_cases, low_risk_cases,
  high_risk_non_cases, low_risk_non_cases,
  high_risk_score2, medium_risk_score2, low_risk_score2,
  high_risk_score2_cases, medium_risk_score2_cases, low_risk_score2_cases,
  high_risk_score2_non_cases, medium_risk_score2_non_cases, low_risk_score2_non_cases,
  reclassified, not_reclassified, reclassified_cases, not_reclassified_cases,
  reclassified_non_cases, not_reclassified_non_cases
)]

# Compute screening statistics
pop_boot[, events_prevented := floor(high_risk_cases * 0.2)] # assuming 20% reduction due to statins
pop_boot[, NNS := ceiling(people / events_prevented)] # Number needed to screen per event prevented
pop_boot[, NNT := ceiling(high_risk / events_prevented)] # Number of statins prescribed per event prevented

pop_boot[, events_prevented_score2 := floor(high_risk_score2_cases * 0.2)] # assuming 20% reduction due to statins
pop_boot[, NNS_score2 := ceiling(people / events_prevented_score2)] # Number needed to screen per event prevented
pop_boot[, NNT_score2 := ceiling(high_risk_score2 / events_prevented_score2)] # Number of statins prescribed per event prevented

pop_boot[, delta_high_risk_cases := high_risk_cases - high_risk_score2_cases]
pop_boot[, delta_events_prevented := events_prevented - events_prevented_score2]

# Write out
fwrite(pop_boot, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening/ons_pop_stratified_with_bootstraps.txt")

###########

# Convert to longer format to make aggregation simpler
pop_boot <- melt(pop_boot, id.vars=c("sex", "age_group", "model", "bootstrap"), variable.name="number", value.name="estimate")
pop_boot <- pop_boot[!(number %in% c("events_prevented", "events_prevented_score2", "NNS", "NNS_score2", "NNT", "NNT_score2"))] # not valid when summed

###########

# Get aggregate statistics by age-group
group_agg <- pop_boot[,.(estimate=sum(estimate)), by=.(bootstrap, sex, age_group, model, number)]
group_agg <- dcast(group_agg, bootstrap + sex + age_group + model ~ number, value.var="estimate")

group_agg[, events_prevented := floor(high_risk_cases * 0.2)] # assuming 20% reduction due to statins
group_agg[, NNS := ceiling(people / events_prevented)] # Number needed to screen per event prevented
group_agg[, NNT := ceiling(high_risk / events_prevented)] # Number of statins prescribed per event prevented

group_agg[, events_prevented_score2 := floor(high_risk_score2_cases * 0.2)] # assuming 20% reduction due to statins
group_agg[, NNS_score2 := ceiling(people / events_prevented_score2)] # Number needed to screen per event prevented
group_agg[, NNT_score2 := ceiling(high_risk_score2 / events_prevented_score2)] # Number of statins prescribed per event prevented

group_agg[, delta_high_risk_cases := high_risk_cases - high_risk_score2_cases]
group_agg[, delta_events_prevented := events_prevented - events_prevented_score2]

group_agg <- melt(group_agg, id.vars=c("sex", "age_group", "people", "cases", "non_cases", "model", "bootstrap"), variable.name="number", value.name="estimate")

# Extract estimates and 95% CIs
group_screen_estimates <- group_agg[bootstrap == 0]
group_screen_estimates[, bootstrap := NULL]

group_screen_95ci <- group_agg[bootstrap != 0]
group_screen_95ci <- group_screen_95ci[order(estimate)]
group_screen_95ci <- group_screen_95ci[, cbind(colname=c("L95", "U95"), .SD[c(25, 975)]), by=.(sex, age_group, people, cases, non_cases, model, number)]
group_screen_95ci <- dcast(group_screen_95ci, sex + age_group + people + cases + non_cases + model + number ~ colname, value.var="estimate")

group_screen <- group_screen_estimates[group_screen_95ci, on = .(sex, age_group, people, cases, non_cases, model, number)]

# Write out
fwrite(group_screen, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening/population_screening_by_sex_and_age_group.txt")

###########

# Get aggregate statistics in males and females
sex_agg <- pop_boot[,.(estimate=sum(estimate)), by=.(bootstrap, sex, model, number)]
sex_agg <- dcast(sex_agg, bootstrap + sex + model ~ number, value.var="estimate")

sex_agg[, events_prevented := floor(high_risk_cases * 0.2)] # assuming 20% reduction due to statins
sex_agg[, NNS := ceiling(people / events_prevented)] # Number needed to screen per event prevented
sex_agg[, NNT := ceiling(high_risk / events_prevented)] # Number of statins prescribed per event prevented

sex_agg[, events_prevented_score2 := floor(high_risk_score2_cases * 0.2)] # assuming 20% reduction due to statins
sex_agg[, NNS_score2 := ceiling(people / events_prevented_score2)] # Number needed to screen per event prevented
sex_agg[, NNT_score2 := ceiling(high_risk_score2 / events_prevented_score2)] # Number of statins prescribed per event prevented

sex_agg[, delta_high_risk_cases := high_risk_cases - high_risk_score2_cases]
sex_agg[, delta_events_prevented := events_prevented - events_prevented_score2]

sex_agg <- melt(sex_agg, id.vars=c("sex", "people", "cases", "non_cases", "model", "bootstrap"), variable.name="number", value.name="estimate")

# Extract estimates and 95% CIs
sex_screen_estimates <- sex_agg[bootstrap == 0]
sex_screen_estimates[, bootstrap := NULL]

sex_screen_95ci <- sex_agg[bootstrap != 0]
sex_screen_95ci <- sex_screen_95ci[order(estimate)]
sex_screen_95ci <- sex_screen_95ci[, cbind(colname=c("L95", "U95"), .SD[c(25, 975)]), by=.(sex, people, cases, non_cases, model, number)]
sex_screen_95ci <- dcast(sex_screen_95ci, sex + people + cases + non_cases + model + number ~ colname, value.var="estimate")

sex_screen <- sex_screen_estimates[sex_screen_95ci, on = .(sex, people, cases, non_cases, model, number)]

# Write out
fwrite(sex_screen, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening/population_screening_by_sex.txt")

###########

# Get aggregate statistics in the total population
pop_agg <- pop_boot[,.(estimate=sum(estimate)), by=.(bootstrap, model, number)]
pop_agg <- dcast(pop_agg, bootstrap + model ~ number, value.var="estimate")

pop_agg[, events_prevented := floor(high_risk_cases * 0.2)] # assuming 20% reduction due to statins
pop_agg[, NNS := ceiling(people / events_prevented)] # Number needed to screen per event prevented
pop_agg[, NNT := ceiling(high_risk / events_prevented)] # Number of statins prescribed per event prevented

pop_agg[, events_prevented_score2 := floor(high_risk_score2_cases * 0.2)] # assuming 20% reduction due to statins
pop_agg[, NNS_score2 := ceiling(people / events_prevented_score2)] # Number needed to screen per event prevented
pop_agg[, NNT_score2 := ceiling(high_risk_score2 / events_prevented_score2)] # Number of statins prescribed per event prevented

pop_agg[, delta_high_risk_cases := high_risk_cases - high_risk_score2_cases]
pop_agg[, delta_events_prevented := events_prevented - events_prevented_score2]

pop_agg <- melt(pop_agg, id.vars=c("people", "cases", "non_cases", "model", "bootstrap"), variable.name="number", value.name="estimate")

# Extract estimates and 95% CIs
pop_screen_estimates <- pop_agg[bootstrap == 0]
pop_screen_estimates[, bootstrap := NULL]

pop_screen_95ci <- pop_agg[bootstrap != 0]
pop_screen_95ci <- pop_screen_95ci[order(estimate)]
pop_screen_95ci <- pop_screen_95ci[, cbind(colname=c("L95", "U95"), .SD[c(25, 975)]), by=.(people, cases, non_cases, model, number)]
pop_screen_95ci <- dcast(pop_screen_95ci, people + cases + non_cases + model + number ~ colname, value.var="estimate")

pop_screen <- pop_screen_estimates[pop_screen_95ci, on = .(people, cases, non_cases, model, number)]

# Write out
fwrite(pop_screen, sep="\t", quote=FALSE, file="analyses/public_health_modelling/targeted_screening/population_screening.txt")

