library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(cowplot)
library(scales)
library(RColorBrewer)
source("src/utils/SCORE2.R")

# Make output directory
system("mkdir -p analyses/CVD_weight_training", wait=TRUE)

# Load required data
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "smoking", "sbp", "tchol", "hdl", "incident_cvd_followup", "incident_cvd", "cvd_prediction_foldid", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2_excl_UKB"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

# Load NMR scores from model training
train_scores <- rbind(idcol="type", 
  "non-derived"=fread("analyses/nmr_score_training/non_derived_NMR_scores.txt"),
  "clinical"=fread("analyses/nmr_score_training/clinical_NMR_scores.txt")
)

# Estimate per-score weights to use when combining with SCORE2 using Cox proportional hazards models
cvd_weights <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        # extract training data
        this_dat <- dat[cvd_prediction_foldid != this_test_fold & sex == this_sex]

        # add in respective nmr scores
        this_train_scores <- train_scores[prediction_cv_testfold == this_test_fold & prediction_cv_testfold != cvd_prediction_foldid & sex == this_sex & type == this_score_type]
        this_dat[this_train_scores[endpoint == "CAD"], on = .(eid), CAD_NMR_score := i.linear_predictor]
        this_dat[this_train_scores[endpoint == "Stroke"], on = .(eid), Stroke_NMR_score := i.linear_predictor]

        # Center scores
        scaling_factors <- melt(this_dat, id.vars=c("eid"), measure.vars=c("CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"), variable.name="score")
        scaling_factors <- scaling_factors[, .(offset=-mean(value)), by=score]

        this_dat[, CAD_NMR_score := CAD_NMR_score - mean(CAD_NMR_score)]
        this_dat[, Stroke_NMR_score := Stroke_NMR_score - mean(Stroke_NMR_score)]
        this_dat[, CAD_metaGRS := CAD_metaGRS - mean(CAD_metaGRS)]
        this_dat[, Stroke_metaGRS := Stroke_metaGRS - mean(Stroke_metaGRS)]
  
        # Build model formula
        mf <- "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2)"
        if (this_model %in% c("SCORE2 + NMR scores", "SCORE2 + NMR scores + PRSs")) {
          mf <- paste(mf, "+ CAD_NMR_score + Stroke_NMR_score")
        }
        if (this_model %in% c("SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")) {
          mf <- paste(mf, "+ CAD_metaGRS + Stroke_metaGRS")
        }
   
        # Fit cox model
        cx <- coxph(as.formula(mf), data=this_dat)

        # Extract relevant coefficients
        ci <- confint(cx)
        cf <- as.data.table(coef(summary(cx)), keep.rownames="score")
        cf <- cf[, .(score, weight=coef, L95=ci[,1], U95=ci[,2], pval=`Pr(>|z|)`)]

        # Add scaling factors
        cf <- scaling_factors[cf, on = .(score)]

        # add information
        cbind(test_fold = this_test_fold, sex = this_sex, model = this_model, score_type = this_score_type, cf)
}

# Write out
fwrite(cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/cvd_score_weights_in_5fold_cross_validation.txt")

# Factor levels for plot ordering
cvd_weights[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
cvd_weights[, score := factor(score, levels=c("CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]

# Plot weight consistencies
plot_weights <- function(this_sex, this_score_type) {
  ggplot(cvd_weights[sex == this_sex & score_type == this_score_type]) +
    aes(x = score, y = weight, ymin = L95, ymax = U95, color = paste("Test fold", test_fold)) +
    facet_grid(. ~ model, space="free_x", scales="free_x") +
    geom_errorbar(position=position_dodge(width=0.75), alpha=0.6, width=0) +
    geom_point(shape=19, alpha=0.6, position=position_dodge(width=0.75)) +
    geom_hline(yintercept=0, linetype=2) +
    xlab("") + ylab("Weight") +
    guides(color=guide_legend(title="")) +
    theme_bw() +
    theme(
      axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
      axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
      strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
      panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
      legend.position="bottom", legend.text=element_text(size=6),
      legend.box.margin=margin(-7,-3,0,-3), legend.margin=margin(-3,-3,-3,-3)
    )
}

g1 <- plot_weights("Male", "non-derived")
g2 <- plot_weights("Female", "non-derived")
g <- plot_grid(g1, g2, nrow=2, labels=c("Males  ", "Females"), label_size=10)
ggsave(g, width=7.2, height=4.5, file="analyses/CVD_weight_training/non_derived_NMR_scores_weights.pdf")

g1 <- plot_weights("Male", "clinical")
g2 <- plot_weights("Female", "clinical")
g <- plot_grid(g1, g2, nrow=2, labels=c("Males  ", "Females"), label_size=10)
ggsave(g, width=7.2, height=4.5, file="analyses/CVD_weight_training/clinical_NMR_scores_weights.pdf")

# Compute combined linear predictors in the test datasets and aggregate
pred_scores <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        # extract test data
        this_dat <- dat[cvd_prediction_foldid == this_test_fold & sex == this_sex]

        # Add in relevant NMR scores predicted in withheld test fold
        this_nmr_scores <- train_scores[prediction_cv_testfold == this_test_fold & prediction_cv_testfold == cvd_prediction_foldid & sex == this_sex & type == this_score_type]
        this_nmr_scores[, score := paste0(endpoint, "_NMR_score")]
        this_nmr_scores <- dcast(this_nmr_scores, eid ~ score, value.var="linear_predictor")
        this_dat <- this_nmr_scores[this_dat, on = .(eid)]

        # Melt columns of interest to long format
        lp_cols <- "SCORE2"
        if (this_model %like% "NMR") {
          lp_cols <- c(lp_cols, "CAD_NMR_score", "Stroke_NMR_score")
        } 
        if (this_model %like% "PRS") {
          lp_cols <- c(lp_cols, "CAD_metaGRS", "Stroke_metaGRS")
        }
        this_dat <- melt(this_dat, id.vars=c("eid", "sex", "age", "incident_cvd", "incident_cvd_followup", "cvd_prediction_foldid"), measure.vars=lp_cols, variable.name="score")

        # Add offsets and multiply scores by weights
        this_weights <- cvd_weights[test_fold == this_test_fold & sex == this_sex & model == this_model & score_type == this_score_type]
        this_dat[this_weights, on = .(score), value := (value + i.offset) * i.weight]

        # Sum to create combined linear predictor
        this_dat <- this_dat[, .(linear_predictor=sum(value)), by = .(eid, age, incident_cvd_followup, incident_cvd)]
      
        # add model info
        cbind(score_type = this_score_type, model = this_model, sex = this_sex, test_fold=this_test_fold, this_dat)
}


# Add five year age group (downstream analyses uses these)
pred_scores[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Code factors
pred_scores[, sex := factor(sex, levels=c("Male", "Female"))]
pred_scores[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs", "SCORE2"))]
pred_scores[, age_group := factor(age_group, levels=c("40-44", "45-49", "50-54", "55-59", "60-64", "65-69"))]

# Compute absolute risk using SCORE2 baseline hazards
pred_scores[, uncalibrated_risk := score2_absrisk(sex, linear_predictor)]

# Compute absolute risk recalibrated to low-risk european region (including UK)
pred_scores[, uk_calibrated_risk := score2_recalibration(sex, uncalibrated_risk, "low")]

# Reorganize columns
pred_scores <- pred_scores[, .(eid, sex, age, age_group, incident_cvd, incident_cvd_followup, cvd_prediction_foldid=test_fold,
  model, score_type, linear_predictor, uncalibrated_risk, uk_calibrated_risk)]

# Write out
fwrite(pred_scores, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Build supp table of avg weights and transformations
dt <- cvd_weights[score_type == "non-derived", .(offset=mean(offset), weight=mean(weight)), by=.(model, sex, score)]
dt <- dcast(dt, score + offset + sex ~ model, value.var="weight", fill=0)
dt <- melt(dt, id.vars=c("score", "sex"))
dt <- dcast(dt, sex + variable ~ score, value.var="value")
dt <- dt[order(-sex)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/avg_cvd_score_weights_supp_table.txt")

# Check calibration of predicted risks in each five-year age-group, i.e. to see whether calibration scaling
# factors used by SCORE2 (i.e. to the low-risk european region) are still applicable
avg_lp <- pred_scores[,.(mean=mean(linear_predictor)), by=.(score_type, model, sex, age_group)]
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

ggsave(g, width=7.2, height=4, file="analyses/CVD_weight_training/average_predicted_LP_by_age_group.pdf")

avg_lp2 <- avg_lp[model != "SCORE2"]
avg_lp2[avg_lp[model == "SCORE2"], on = .(score_type, sex, age_group), SCORE2 := i.mean]

g <- ggplot(avg_lp2[score_type == "non-derived"]) +
  aes(x=mean, y=SCORE2, color=model) +
  facet_grid(sex ~ model) +
  geom_abline(intercept=0, slope=1, linetype=2, linewidth=0.6) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.6) +
  geom_point(shape=23, size=1.2, fill="white") +
  scale_y_continuous("SCORE2-LP (average per 5-year age-group)") +
  scale_x_continuous("new LP (average per 5-year age-group)") +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(), 
    strip.text=element_text(size=7, face="bold"), legend.position="none"
  )

ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/new_LP_vs_SCORE2_LP_average_per_age_group.pdf")

# Plot eCDF for the absolute risk predicted for each score
absrisk <- pred_scores[,.(eid, age, age_group, sex, incident_cvd, score_type, model, absrisk=uk_calibrated_risk)]
absrisk[, eCDF := ecdf(absrisk)(absrisk), by=.(score_type, incident_cvd, sex, model)]
absrisk[, case_status := ifelse(incident_cvd, "CVD cases", "Non cases")]

g <- ggplot(absrisk[score_type == "non-derived"]) +
  aes(x=absrisk, y=eCDF, color=model) +
  facet_grid(sex ~ case_status) +
  geom_line(linewidth=0.6) +
  scale_x_continuous("Predicted 10-year CVD risk (calibrated to low-risk region)", labels=scales::percent) +
  scale_y_continuous("Pr(X <= x)") + 
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )

ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/absrisk_CDFs.pdf")

# Create hexbins of absrisk
absrisk2 <- absrisk[model != "SCORE2"]
score2_absrisk <- absrisk[model == "SCORE2" & score_type == "non-derived"]
absrisk2[score2_absrisk, on = .(eid), SCORE2 := i.absrisk]

g1 <- ggplot(absrisk2[score_type == "non-derived" & case_status == "CVD cases"]) +
  aes(x=absrisk, y=SCORE2) +
  facet_grid(sex ~ model) +
  geom_hex() + 
  geom_abline(intercept=0, slope=1, linetype=2) +
  scale_fill_gradientn("CVD cases", colors=brewer.pal(9, "YlOrRd")[2:8], trans="log10") +
  scale_x_continuous("10-year CVD risk predicted by new model", labels=scales::percent) +
  scale_y_continuous("10-year CVD risk predicted by SCORE2", labels=scales::percent) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=6, face="bold"), legend.position="bottom",
    legend.text=element_text(size=6), legend.title=element_text(size=8)
  )

g2 <- ggplot(absrisk2[score_type == "non-derived" & case_status == "Non cases"]) +
  aes(x=absrisk, y=SCORE2) +
  facet_grid(sex ~ model) +
  geom_hex() + 
  geom_abline(intercept=0, slope=1, linetype=2) +
  scale_fill_gradientn("Non-cases", colors=brewer.pal(9, "YlGnBu")[4:8], trans="log10") +
  scale_x_continuous("10-year CVD risk predicted by new model", labels=scales::percent) +
  scale_y_continuous("10-year CVD risk predicted by SCORE2", labels=scales::percent) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=6, face="bold"), legend.position="bottom",
    legend.text=element_text(size=6), legend.title=element_text(size=8)
  )

g <- plot_grid(g1, g2, labels=c("CVD cases", "Non-cases"), label_size=8)
ggsave(g, width=10, height=4.5, file="analyses/CVD_weight_training/absrisk_hexbin_comparison.pdf")

# Create winzorized version with combined absrisk
absrisk3 <- copy(absrisk)
absrisk3[, case_status := "Combined"]
absrisk3 <- rbind(absrisk, absrisk3)
absrisk3[absrisk > 0.3, absrisk := 0.3]
absrisk3[, eCDF := ecdf(absrisk)(absrisk), by=.(score_type, sex, model, case_status)]

g <- ggplot(absrisk3[score_type == "non-derived"]) +
  aes(x=absrisk, y=1-eCDF, color=model) +
  facet_grid(sex ~ case_status) +
  geom_line(linewidth=0.6) +
  scale_x_continuous("Predicted 10-year CVD risk (calibrated to low-risk region)", labels=scales::percent) +
  scale_y_continuous("Pr(X > x)") + 
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )

ggsave(g, width=7.2, height=4, file="analyses/CVD_weight_training/absrisk_CDFs_winzorized.pdf")

# Create version with people younger than 50 shifted +2.5% risk so we can add cut-offs
absrisk4 <- copy(absrisk3)
absrisk4[age < 50, absrisk := absrisk + 0.025]
absrisk4[absrisk > 0.3, absrisk := 0.3]

g <- ggplot(absrisk3[score_type == "non-derived"]) +
  aes(x=absrisk, y=1-eCDF, color=model) +
  facet_grid(sex ~ case_status) +
  geom_line(linewidth=0.6) +
  geom_vline(xintercept=0.1, linetype=2, color="#fc4e2a", linewidth=0.6) +
  geom_vline(xintercept=0.05, linetype=2, color="#feb24c", linewidth=0.6) +
  scale_x_continuous("Predicted 10-year CVD risk (calibrated to low-risk region)", labels=scales::percent) +
  scale_y_continuous("Pr(X > x)") +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  theme_bw() +
  theme(
    axis.text=element_text(size=6), axis.title=element_text(size=8), strip.background=element_blank(),
    strip.text=element_text(size=7, face="bold"), legend.position="bottom",
    legend.text=element_text(size=7), legend.title=element_blank()
  )

ggsave(g, width=7.2, height=4, file="analyses/CVD_weight_training/absrisk_CDFs_with_thresholds.pdf")

