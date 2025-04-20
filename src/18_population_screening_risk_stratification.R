library(data.table)
library(foreach)
library(boot)
options(boot.parallel="multicore")
options(boot.ncpus=20) # Takes almost 6 hours with 20 cores on icelake

# create output directory
system("mkdir -p analyses/public_health_modelling")

# Load in predicted risks
pred_risk <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Create a mapping of model names to model numbers (which will be used as column names
# during bootstrap)
model_map <- rbind(use.names=FALSE,
  data.table(model_type="", model_name="Risk score", model_colname="model1"),
  data.table("NMR scores", "Risk score + NMR scores", "model2"),
  data.table("Biochemistry", "Risk score + clinical biomarkers", "model3"),
  data.table("PRS", "Risk score + PRSs", "model4"),
  data.table("NMR scores + PRS", "Risk score + NMR scores + PRSs", "model5"),
  data.table("Biochemistry + PRS", "Risk score + clinical biomarkers + PRSs", "model6")
)
pred_risk[model_map, on = .(model_type), model_colname := model_colname]

# Run through all analyses
res <- foreach(this_strategy = c("blanket", "targeted"), .combine=rbind) %:%
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

      # If we're doing targeted screening, only use the alternative models for risk stratification
      # in the medium risk group
      if (this_strategy == "targeted") {
        this_dat[this_dat[model_type == "" & risk_group != "medium"], on = .(eid), risk_group := i.risk_group]
      } 
 
      # For QRISK3, set the medium risk group back to low: a medium risk group is not defined by
      # the NICE 2023 guidelines, we just created our own one for targeted screening
      if (this_score == "QRISK3") {
        this_dat[risk_group == "medium", risk_group := "low"]
      }

      # Get wide format for simultanous bootstrapping of all models
      this_dat <- dcast(this_dat, eid + sex + age_group + incident_cvd + incident_cvd_followup ~ model_colname, value.var="risk_group", fill="missing")

      # The key numbers we need to calculate in each bootstrap resampling are:
      #  (1) Within each age- and sex- group, what proportion of CVD cases are allocated to each risk group?
      #  (2) Within each age- and sex- group, what proportion of non cases are allocated to each risk group?
      # We need to do this (1) for the risk score alone, and (2) for the risk score and alternate model in
      # the shared samples with non-missing data for each alternate model.

      # The bootstrap function just returns a single unnamed vector, so we need to set up tracking for the statistics
      # we want to compute at each bootstrap
      stat_info <- foreach(this_model_type=model_map$model_type, .combine=rbind) %:% 
        foreach(this_comp = c("ref", "alt"), .combine=rbind) %:%
					foreach(this_sex = sort(unique(this_dat$sex)), .combine=rbind) %:% 
						foreach(this_age_group = sort(unique(this_dat$age_group)), .combine=rbind) %:% 
							foreach(this_risk_group = c("low", "medium", "high"), .combine=rbind) %:%
								foreach(this_status = c("case", "non-case"), .combine=rbind) %do% {
									if (this_score == "QRISK3" & this_risk_group == "medium") return(NULL)
                  if (this_model_type == "" & this_comp == "alt") return(NULL)
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
          } else {
					 	n <- sum(this_dt[[stat_info[statIdx, model_colname]]] == stat_info[statIdx, risk_group])
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
fwrite(res, sep="\t", quote=FALSE, file="analyses/public_health_modelling/discovery_standardised_risk_stratification_bootstraps.txt")

