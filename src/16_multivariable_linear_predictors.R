library(data.table)
library(foreach)
library(ggplot2)
source('src/utils/SCORE2.R')
source('src/utils/QRISK3.R')

# Load data
dat <- fread("data/cleaned/analysis_cohort.txt")
coef <- fread("analyses/CVD_weight_training/multivariable_model_weights.txt")

# Load and add predicted NMR scores in the discovery data
nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]

# Add five year age group (downstream analyses uses these)
dat[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Compute linear predictors for main and sensitivity analyses
models <- c("", "NMR scores", "Biochemistry", "PRS", "NMR scores + PRS", "Biochemistry + PRS")
pred_scores <- foreach(this_model=models, .combine=rbind) %:% 
  foreach(this_endpoint=c("cvd", "cvd_narrow"), .combine=rbind) %:%
		foreach(this_score=c("SCORE2", "SCORE2_excl_UKB", "QRISK3"), .combine=rbind) %:%
			foreach(this_sex=c("Male", "Female"), .combine=rbind) %do% {
				# Extract relevant data common to all models
        this_dat <- dat[sex == this_sex, .(eid, age, age_group)]

        # Add in relevant CVD endpoint
        if (this_endpoint == "cvd") {
          this_dat[dat, on = .(eid), c("incident_cvd", "incident_cvd_followup") := .(i.incident_cvd, i.incident_cvd_followup)]
        } else {
          this_dat[dat, on = .(eid), c("incident_cvd", "incident_cvd_followup") := .(i.incident_cvd2, i.incident_cvd2_followup)]
        }
  
				# Extract relevant baseline risk score to start building the linear predictor
        if (this_score == "SCORE2") {
					this_dat[dat, on = .(eid), linear_predictor := SCORE2] 
        } else if (this_score == "SCORE2_excl_UKB") {
					this_dat[dat, on = .(eid), linear_predictor := SCORE2_excl_UKB] 
        } else {
					this_dat[dat, on = .(eid), linear_predictor := QRISK3] 
        }

        # Filter to people with non-missing data (only applies to QRISK3)
        this_dat <- this_dat[!is.na(linear_predictor)]

        # If the model is not just the risk score on its own, we need to add based on the multivariable model fits
        if (this_model != "") {
          # Extract relavent variables
					this_coef <- coef[model == paste("SCORE2 +", this_model) & sex == this_sex]
					var_dat <- dat[sex == this_sex, .SD, .SDcols=c("eid", this_coef[, variable_col])]

					# Filter to samples with all required data
					var_dat <- var_dat[complete.cases(var_dat)]
					this_dat <- this_dat[var_dat[,.(eid)], on = .(eid), nomatch=0] # match row order as well
          var_dat <- var_dat[this_dat[,.(eid)], on = .(eid), nomatch=0] # remove people with missing QRISK3

					# Add relevant coefficients to the linear predictor
					for (this_var in this_coef$variable_col) {
						this_dat[, linear_predictor := linear_predictor + scale(var_dat[[this_var]]) * this_coef[variable_col == this_var, log(HR)]]
					}
        }

        # Add in model information
        model_info <- data.table(endpoint=this_endpoint, sex=this_sex, score=this_score, model_type=this_model)
        if (this_model == "") {
          model_info[, model := this_score]
        } else {
          model_info[, model := paste(this_score, "+", this_model)]
        }
        model_info[, model := gsub("SCORE2_excl_UKB", "SCORE2", model)]
       
        return(cbind(model_info, this_dat))
}

# Compute absolute risk using SCORE2 baseline hazards
pred_scores[score != "QRISK3", uncalibrated_risk := score2_absrisk(sex, linear_predictor)]

# Compute absolute risk recalibrated to low-risk europ4ean region (including UK)
pred_scores[score != "QRISK3", uk_calibrated_risk := score2_recalibration(sex, uncalibrated_risk, "low")]

# And do the same for QRISK3
pred_scores[score == "QRISK3", uk_calibrated_risk := QRISK3_absrisk(sex, linear_predictor)]

# Write out
fwrite(pred_scores, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Check calibration of predicted risks in each five-year age-group, i.e. to see whether calibration scaling
# factors used by SCORE2 (i.e. to the low-risk european region) are still applicable
avg_lp <- pred_scores[endpoint == "cvd", .(mean=mean(linear_predictor)), by=.(score, model_type, sex, age_group)]
avg_lp[, age_group_start := as.integer(gsub("-.*", "", age_group))]
avg_lp[, model_type := paste("Risk score +", model_type)]
avg_lp[model_type == "Risk score + ", model_type := "Risk score"]

g <- ggplot(avg_lp) +
  aes(x=age_group_start, y=mean, color=model_type) +
  facet_grid(sex ~ score) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.6) +
  geom_point(shape=23, size=1.2, fill="white") +
  scale_x_continuous("Age group", breaks=seq(40, 65, by=5), labels=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69")) +
  scale_y_continuous("Average linear predictor") +
  scale_color_manual(values=c(
    "Risk score"="black", "Risk score + NMR scores"="#e41a1c", "Risk score + PRS"="#377eb8", "Risk score + Biochemistry"="#ff7f00",
    "Risk score + NMR scores + PRS"="#4daf4a", "Risk score + Biochemistry + PRS"="#984ea3"
  )) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )

ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/average_predicted_LP_by_age_group.pdf")

avg_lp2 <- avg_lp[model_type != "Risk score"]
avg_lp2[avg_lp[model_type == "Risk score"], on = .(score, sex, age_group), Risk_score := i.mean]

g <- ggplot(avg_lp2) +
  aes(x=mean, y=Risk_score, color=model_type) +
  facet_grid(sex ~ score) +
  geom_abline(intercept=0, slope=1, linetype=2, linewidth=0.6) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.6) +
  geom_point(shape=23, size=1.2, fill="white") +
  scale_y_continuous("Risk score LP (average per 5-year age-group)") +
  scale_x_continuous("New LP (average per 5-year age-group)") +
  scale_color_manual(values=c(
    "Risk score"="black", "Risk score + NMR scores"="#e41a1c", "Risk score + PRS"="#377eb8", "Risk score + Biochemistry"="#ff7f00",
    "Risk score + NMR scores + PRS"="#4daf4a", "Risk score + Biochemistry + PRS"="#984ea3"
  )) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )

ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/new_LP_vs_SCORE2_LP_average_per_age_group.pdf")


