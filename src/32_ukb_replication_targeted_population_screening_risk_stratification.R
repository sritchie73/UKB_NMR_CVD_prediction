library(data.table)
library(foreach)
library(boot)
options(boot.parallel="multicore")
options(boot.ncpus=20) # Takes under 6 hours with 20 cores on sapphire

# create output directory
system("mkdir -p analyses/public_health_modelling")

# Load in predicted risks
pred_risk <- fread("analyses/CVD_weight_training/phase3_CVD_linear_predictors_and_risk.txt")

# Create a mapping of model names to model numbers (which will be used as column names
# during bootstrap)
model_map <- rbind(use.names=FALSE,
  data.table(model_type="", model_name="Risk score", model_colname="model1"),
  data.table("NMR scores", "Risk score + NMR scores", "model2"),
  data.table("Biochemistry", "Risk score + clinical biomarkers", "model3"),
  data.table("PRS", "Risk score + PRSs", "model4"),
  data.table("NMR scores + PRS", "Risk score + NMR scores + PRSs", "model5"),
  data.table("Biochemistry + PRS", "Risk score + clinical biomarkers + PRSs", "model6"),
  data.table("NMR scores + Biochemistry", "Risk score + NMR scores + clinical biomarkers", "model7"),
  data.table("NMR scores + Biochemistry + PRS", "Risk score + NMR scores + clinical biomarkers + PRSs", "model8")
)
pred_risk[model_map, on = .(model_type), model_colname := model_colname]

# Run through all analyses
res <- foreach(this_strategy = c("targeted"), .combine=rbind) %:%
	foreach(this_endpoint = c("cvd", "cvd_narrow"), .combine=rbind) %:%
		foreach(this_score = c("SCORE2", "SCORE2_excl_UKB", "QRISK3"), .combine=rbind) %do% {
			# Extract all scores for this inner loop
			this_dat <- pred_risk[endpoint == this_endpoint & score == this_score]
     
      # Assign risk strata
      if (this_score == "QRISK3") {
        this_dat[, risk_group := fcase(
          uk_calibrated_risk < 0.05, "low",
          uk_calibrated_risk < 0.10, "medium",
          default="high")]
      } else {
        this_dat[, risk_group := fcase(
					age < 50 & uk_calibrated_risk < 0.025, "low",
					age < 50 & uk_calibrated_risk < 0.075, "medium",
					age < 50, "high",
					uk_calibrated_risk < 0.05, "low",
					uk_calibrated_risk < 0.10, "medium",
					default="high")]
      }

      # For targeted screening only use the alternative models for risk stratification
      # in the medium risk group, and only for increasing risk to high risk with targeted
      # rescreening
      this_dat[this_dat[model_type == ""], on = .(eid), pop_screen := i.risk_group]
      this_dat[model_type != "" & pop_screen == "medium" & risk_group == "high", targeted_screen := "high"]
      this_dat[is.na(targeted_screen), targeted_screen := pop_screen]

      # For QRISK3, set the medium risk group back to low: a medium risk group is not defined by
      # the NICE 2023 guidelines, we just created our own one for targeted screening
      if (this_score == "QRISK3") {
        this_dat[risk_group == "medium", risk_group := "low"]
      }

      # We also need to create a risk group column for our null hypothesis testing - whereas in population-wide screening our null hypothesis was
      # risk stratification equivalent to SCORE2 alone, now our null hypothesis is risk stratification equivalent to adding the same number of 
      # additional people using the SCORE2 rankings to the high risk group. Our null hypothesis will be computed on a sex-specific basis since we
      # do sex-specific analyses - and also by age-group to make the comparison fairer (otherwise everyone added to the high risk group in the 
      # comparison will be overwhelmingly male and in the oldest age group).
      null_dat <- this_dat[model_type == "", .(eid, sex, age_group, null_risk_group=risk_group, null_risk=uk_calibrated_risk)]
      null_dat <- null_dat[order(-null_risk)]

      # Get wide format for simultanous bootstrapping of all models
      this_dat <- dcast(this_dat, eid + sex + age_group + incident_cvd + incident_cvd_followup ~ model_colname, value.var="targeted_screen", fill="missing")

      # Missing data in the alternate model is only relevant in the re-screening group
      this_dat[model2 == "missing" & model1 != "medium", model2 := model1]
      this_dat[model3 == "missing" & model1 != "medium", model3 := model1]
      this_dat[model4 == "missing" & model1 != "medium", model4 := model1]
      this_dat[model5 == "missing" & model1 != "medium", model5 := model1]
      this_dat[model6 == "missing" & model1 != "medium", model6 := model1]
      this_dat[model7 == "missing" & model1 != "medium", model7 := model1]
      this_dat[model8 == "missing" & model1 != "medium", model8 := model1]

      # Build null hypothesis risk stratification column for each model 
      for (this_model in model_map[model_type != "", model_colname]) {
        # How many people were additionally assinged to the high risk category when doing the targeted re-screening?
        alt_high_risk <- this_dat[, .SD, .SDcols=c("model1", this_model)]
        setnames(alt_high_risk, c("ref", "alt"))
        alt_high_risk <- alt_high_risk[ref != "high" & alt == "high", .N]
  
        # Adding that same number of people to high risk based on ranking of SCORE2/QRISK3 - while also dropping people who don't have data for 
        # this particular model
        this_missing <- this_dat[,.SD, .SDcols=c("eid", this_model)]
        setnames(this_missing, this_model, "alt")
        this_missing <- this_missing[alt == "missing"]
          
        this_null_dat <- copy(null_dat)
        this_null_dat <- this_null_dat[!this_missing, on = .(eid)]

				this_null_dat[, rank := 0] 
				this_null_dat[null_risk_group != "high", rank := 1:.N]
        this_null_dat[rank <= alt_high_risk, null_risk_group := "high"]

        # add null hypothesis risk grouping for this model to the main this_dat
        this_dat[this_null_dat, on = .(eid), null := null_risk_group]
        this_dat[is.na(null), null := "missing"]
        setnames(this_dat, "null", sprintf("%s_null", this_model))
      }
     
      # The key numbers we need to calculate in each bootstrap resampling are:
      #  (1) Within each age- and sex- group, what proportion of CVD cases are allocated to each risk group?
      #  (2) Within each age- and sex- group, what proportion of non cases are allocated to each risk group?
      # We need to do this (1) for the risk score alone, (2) for the risk score and alternate model in
      # the shared samples with non-missing data for each alternate model, and (3) in the null model specifically
      # set up for the targeted re-screening analyses.

      # The bootstrap function just returns a single unnamed vector, so we need to set up tracking for the statistics
      # we want to compute at each bootstrap
      stat_info <- foreach(this_model_type=model_map$model_type, .combine=rbind) %:% 
        foreach(this_comp = c("ref", "alt", "null"), .combine=rbind) %:%
					foreach(this_sex = sort(unique(this_dat$sex)), .combine=rbind) %:% 
						foreach(this_age_group = sort(unique(this_dat$age_group)), .combine=rbind) %:% 
							foreach(this_risk_group = c("low", "medium", "high"), .combine=rbind) %:%
								foreach(this_status = c("case", "non-case"), .combine=rbind) %do% {
									if (this_score == "QRISK3" & this_risk_group == "medium") return(NULL)
                  if (this_model_type == "" & this_comp != "ref") return(NULL)
									data.table(model_type=this_model_type, comparitor=this_comp, sex=this_sex, age_group=this_age_group, risk_group=this_risk_group, status=this_status)
      }
      stat_info <- model_map[stat_info, on=.(model_type)]

      # Function for the bootstrap analysis to run
      boot_fun <- function(dt) {
        foreach(statIdx = stat_info[,.I], .combine=c) %do% {
          # Filter to complete data for both models being compared
          this_dt <- dt[model1 != "missing"]
          if (stat_info[statIdx, model_type != ""]) {
            this_dt <- dt[dt[[stat_info[statIdx, model_colname]]] != "missing"]
          }

          # Filter to relevant age and sex group
          this_dt <- this_dt[age_group == stat_info[statIdx, age_group]]
          this_dt <- this_dt[sex == stat_info[statIdx, sex]]
   
          # Filter to relevant case status
          if (stat_info[statIdx, status] == "case") {
            this_dt <- this_dt[(incident_cvd)]
          } else {
            this_dt <- this_dt[!(incident_cvd)]
          }

          # How many people across all risk strata?
          total <- this_dt[,.N]
 
          # How many people in the relevant risk strata?
          if (stat_info[statIdx, comparitor] == "ref") {
					 	n <- this_dt[model1 == stat_info[statIdx, risk_group], .N] 
          } else if (stat_info[statIdx, comparitor] == "alt")  {
					 	n <- sum(this_dt[[stat_info[statIdx, model_colname]]] == stat_info[statIdx, risk_group])
          } else {
            n <- sum(this_dt[[stat_info[statIdx, sprintf("%s_null", model_colname)]]] == stat_info[statIdx, risk_group])
          } 

          # Return percentage allocated, or NA if group size is 0 in the bootstrap for some reason
          if (total == 0) { 
            return(NA)
          } else {
						return(n/total)
          }
        } 
      }

			# Run bootstrap analysis
			surv_cols_idx <- match(c("incident_cvd_followup", "incident_cvd"), names(this_dat))
      n_boot <- 1000
			boot_res <- censboot(this_dat, boot_fun, n_boot, index=surv_cols_idx)

			# Extract bootstrap statistics
			boot_stats <- foreach(this_bootstrap = 0:n_boot, .combine=rbind) %do% {
				this_res <- cbind(bootstrap=this_bootstrap, stat_info)
				if(this_bootstrap == 0) {
					this_res[, proportion := boot_res$t0]
				} else {
					this_res[, proportion := boot_res$t[this_bootstrap,]]
				}
				return(this_res)
			}
     
      # add in analysis info and return
      cbind(strategy=this_strategy, endpoint=this_endpoint, score=this_score, boot_stats)
}

# Apply demographic standardisation factors
ons_pop <- fread("analyses/public_health_modelling/ONS_hypothetical_100k_pop_by_age_sex.txt")
ons_pop <- melt(ons_pop, id.vars=c('sex', 'age_group'), measure.vars=c("cases", "controls"), variable.name="status", value.name="N")
ons_pop[, status := ifelse(status == "cases", "case", "non-case")]
res[ons_pop, on = .(sex, age_group, status), N := proportion * N]

# Write out bootstrap statistics
fwrite(res, sep="\t", quote=FALSE, file="analyses/public_health_modelling/replication_standardised_targeted_risk_stratification_bootstraps.txt")

