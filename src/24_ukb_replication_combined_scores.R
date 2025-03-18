library(data.table)
library(foreach)
source("src/utils/SCORE2.R")

# Load analaysis cohort
pheno <- fread("data/cleaned/phase3_analysis_cohort.txt")

# Add in NMR scores
nmr_scores <- fread("analyses/nmr_score_training/phase3_nmr_scores.txt")
pheno <- pheno[nmr_scores, on = .(eid, sex)]

# Load combining weights
comb_weights <- fread("analyses/CVD_weight_training/avg_cvd_score_weights_supp_table.txt")
comb_weights <- melt(comb_weights, id.var=c("sex", "variable"), variable.name="score", value.name="weight")
scaling_means <- comb_weights[variable == "offset"]
comb_weights <- comb_weights[variable != "offset"]
comb_weights[scaling_means, on = .(sex, score), centering_offset := i.weight]
setnames(comb_weights, "variable", "model")
comb_weights <- comb_weights[weight != 0]
comb_weights <- comb_weights[order(model)][order(sex)]

# Compute combined scores:
comb_scores <- foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
  foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
		# Extract relevant data
		this_dat <- pheno[sex == this_sex]
		this_dat <- this_dat[,.SD,.SDcols=c("eid", "age", "sex", "incident_cvd", "incident_cvd_followup", "SCORE2_excl_UKB",
			"CAD_metaGRS", "Stroke_metaGRS", "CAD_NMR_score", "Stroke_NMR_score")]
		if (!(this_model %like% "PRSs")) {
			this_dat[, CAD_metaGRS := NULL]
			this_dat[, Stroke_metaGRS := NULL]
		}
		if (!(this_model %like% "NMR")) {
			this_dat[, CAD_NMR_score := NULL]
			this_dat[, Stroke_NMR_score := NULL]
		}
		this_dat <- melt(this_dat, id.vars=c("eid", "age", "sex", "incident_cvd", "incident_cvd_followup"), variable.name="score")

		this_comb_weights <- comb_weights[model == this_model & sex == this_sex]

		# Apply relevant centering/scaling
	  this_dat[this_comb_weights, on = .(score), value := value + i.centering_offset]

		# Apply relavent per score weighting
	  this_dat[this_comb_weights, on = .(score), value := value * i.weight]

		# Sum to get final combined score
		this_comb_score <- this_dat[,.(linear_predictor=sum(value)),by=.(eid, age, incident_cvd, incident_cvd_followup)]

		# Add in information and return
		info <- data.table(sex=this_sex, model=this_model)
		cbind(info, this_comb_score)
	}
}

# Add five year age group (downstream analyses uses these)
comb_scores[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Compute absolute risk using SCORE2 baseline hazards
comb_scores[, uncalibrated_risk := score2_absrisk(sex, linear_predictor)]

# Compute absolute risk recalibrated to low-risk european region (including UK)
comb_scores[, uk_calibrated_risk := score2_recalibration(sex, uncalibrated_risk, "low")]

# Reorganize columns
comb_scores <- comb_scores[, .(eid, sex, age, age_group, incident_cvd, incident_cvd_followup, cvd_prediction_foldid=0,
  model, score_type="non-derived", linear_predictor, uncalibrated_risk, uk_calibrated_risk)]

# Write out
fwrite(comb_scores, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/phase3_CVD_linear_predictors_and_risk.txt")

# Check calibration of predicted risks in each five-year age-group, i.e. to see whether calibration scaling
# factors used by SCORE2 (i.e. to the low-risk european region) are still applicable
avg_lp <- comb_scores[,.(mean=mean(linear_predictor)), by=.(score_type, model, sex, age_group)]
avg_lp[, age_group_start := as.integer(gsub("-.*", "", age_group))]

g <- ggplot(avg_lp[score_type == "non-derived"]) +
  aes(x=age_group_start, y=mean, color=model) +
  facet_wrap(~ sex) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.6) +
  geom_point(shape=23, size=1.2, fill="white") +
  scale_x_continuous("Age group", breaks=seq(40, 65, by=5), labels=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69")) +
  scale_y_continuous("Average linear predictor") +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )

ggsave(g, width=7.2, height=4, file="analyses/CVD_weight_training/phase3_average_predicted_LP_by_age_group.pdf")

avg_lp2 <- avg_lp[model != "SCORE2"]
avg_lp2[avg_lp[model == "SCORE2"], on = .(score_type, sex, age_group), SCORE2 := i.mean]

g <- ggplot(avg_lp2[score_type == "non-derived"]) +
  aes(y=mean, x=SCORE2, color=model) +
  facet_grid(sex ~ model) +
  geom_abline(intercept=0, slope=1, linetype=2, linewidth=0.6) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.6) +
  geom_point(shape=23, size=1.2, fill="white") +
  scale_x_continuous("SCORE2") +
  scale_y_continuous("SCORE2 + NMR and/or PRSs") +
  ggtitle("Points show averages in each 5-year age group") +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="none", plot.title=element_text(size=8)
  )

ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/phase3_new_LP_vs_SCORE2_LP_average_per_age_group.pdf")

