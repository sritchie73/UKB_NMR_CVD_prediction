library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(cowplot)
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
        cf[, score := gsub("_", " ", score)]
        cf <- cf[, .(score, weight=coef, L95=ci[,1], U95=ci[,2], pval=`Pr(>|z|)`)]

        # add information
        cbind(test_fold = this_test_fold, sex = this_sex, model = this_model, score_type = this_score_type, cf)
}

# Write out
fwrite(cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/cvd_score_weights_in_5fold_cross_validation.txt")

# Factor levels for plot ordering
cvd_weights[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
cvd_weights[, score := factor(score, levels=c("CAD NMR score", "Stroke NMR score", "CAD metaGRS", "Stroke metaGRS"))]

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

# Compute combined linear predictors in the test datasets and aggregate, then compare distributions.
# We need to make sure the distributions are similar to the original SCORE2 linear predictor for the
# downstream recalibration scaling factors published in the SCORE2 paper are still valid.
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

        # Multiply scores by weights
        this_weights <- cvd_weights[test_fold == this_test_fold & sex == this_sex & model == this_model & score_type == this_score_type]
        this_weights[, score := gsub(" ", "_", score)]
        this_dat[this_weights, on = .(score), value := value * i.weight]

        # Sum to create combined linear predictor
        this_dat <- this_dat[, .(linear_predictor=sum(value)), by = .(eid, age, incident_cvd_followup, incident_cvd)]
      
        # add model info
        cbind(score_type = this_score_type, model = this_model, sex = this_sex, test_fold=this_test_fold, this_dat)
}

pred_scores[, sex := factor(sex, levels=c("Male", "Female"))]
pred_scores[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs", "SCORE2"))]

pred_score_means <- pred_scores[,.(mean=mean(linear_predictor), sd=sd(linear_predictor)), by=.(score_type, model, sex)]
pred_score_means[, LSD := mean-sd]
pred_score_means[, USD := mean+sd]

g <- ggplot(pred_scores[score_type == "non-derived"]) +
  aes(x=linear_predictor, color=model) +
  facet_wrap(~ sex) +
  geom_density(trim=TRUE) +
  geom_vline(data=pred_score_means[score_type == "non-derived"], aes(xintercept=mean, color=model), linetype=2) +
  geom_rect(data=pred_score_means[score_type == "non-derived"], inherit.aes=FALSE, aes(ymin=-Inf, ymax=Inf, xmin=LSD, xmax=USD, fill=model), alpha=0.2) +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  scale_fill_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  xlab("Linear predictor") +
  ylab("Density") +
  theme_bw() +
  theme(
    axis.title=element_text(size=8), axis.text=element_text(size=6), 
    legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=7), 
    strip.background=element_blank(), strip.text=element_text(size=7, face="bold")
  )
ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/linear_predictor_densities_without_scaling.pdf")

# Compare absolute risks
pred_scores[, uncalibrated_risk := score2_absrisk(sex, linear_predictor)]
pred_scores[, uk_calibrated_risk := score2_recalibration(sex, uncalibrated_risk, "low")]

pred_score_means <- pred_scores[,.(mean=mean(uk_calibrated_risk), sd=sd(uk_calibrated_risk)), by=.(score_type, model, sex)]
pred_score_means[, LSD := mean-sd]
pred_score_means[, USD := mean+sd]

g <- ggplot(pred_scores[score_type == "non-derived"]) +
  aes(x=uk_calibrated_risk, color=model) +
  facet_wrap(~ sex) +
  geom_density(trim=TRUE) +
  geom_vline(data=pred_score_means[score_type == "non-derived"], aes(xintercept=mean, color=model), linetype=2) +
  geom_rect(data=pred_score_means[score_type == "non-derived"], inherit.aes=FALSE, aes(ymin=-Inf, ymax=Inf, xmin=LSD, xmax=USD, fill=model), alpha=0.2) +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  scale_fill_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  xlab("Calibrated 10-year CVD risk") +
  ylab("Density") +
  theme_bw() +
  theme(
    axis.title=element_text(size=8), axis.text=element_text(size=6), 
    legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=7), 
    strip.background=element_blank(), strip.text=element_text(size=7, face="bold")
  )
ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/uk_calibrated_risk_densities_without_scaling.pdf")

# Next calculate scaling factors in 5-fold cross-validation. To do this, we need to calculate the combined linear
# predictors in the training subsets:
train_lps <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        # extract training data
        this_dat <- dat[cvd_prediction_foldid != this_test_fold & sex == this_sex]

        # Add in relevant NMR scores 
        this_nmr_scores <- train_scores[prediction_cv_testfold == this_test_fold & prediction_cv_testfold != cvd_prediction_foldid & sex == this_sex & type == this_score_type]
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

        # Multiply scores by weights
        this_weights <- cvd_weights[test_fold == this_test_fold & sex == this_sex & model == this_model & score_type == this_score_type]
        this_weights[, score := gsub(" ", "_", score)]
        this_dat[this_weights, on = .(score), value := value * i.weight]

        # Sum to create combined linear predictor
        this_dat <- this_dat[, .(linear_predictor=sum(value)), by = .(eid)]
      
        # add model info
        cbind(score_type = this_score_type, model = this_model, sex = this_sex, test_fold = this_test_fold, this_dat)
}

# Now we need to calculate scaling factors to give each linear predictor the same mean and sd as SCORE2
scaling_factors <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        this_lp <- train_lps[score_type == this_score_type & model == this_model & sex == this_sex & test_fold == this_test_fold]
        score2_lp <- train_lps[score_type == this_score_type & model == "SCORE2" & sex == this_sex & test_fold == this_test_fold]

        # Fit relationship between SCORE2-LP and new-LP in the same way as recalibration to make sure resulting absolute risks
        # are calibrated to the SCORE2 low-risk population scaling factors
        comb_lp <- this_lp[score2_lp, on = .(eid), .(eid, SCORE2_LP=i.linear_predictor, new_LP=linear_predictor)]
        comb_lp[dat, on = .(eid), age := i.age]
        comb_lp[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]
        comb_lp <- comb_lp[, .(SCORE2_LP=mean(SCORE2_LP), new_LP=mean(new_LP)), by=age_group]
        l1 <- lm(comb_lp$SCORE2_LP ~ comb_lp$new_LP)
        scaling_factor <- l1$coefficients[2]
        offset <- l1$coefficients[1]

        # Return
        data.table(score_type=this_score_type, model=this_model, sex=this_sex, test_fold=this_test_fold, scale=scaling_factor, offset=offset)
}

# Apply scaling factors to predicted scores
pred_scores[scaling_factors, on = .(score_type, model, sex, test_fold), linear_predictor := linear_predictor * scale + offset]

# Visualise distributions after transformation
pred_score_means <- pred_scores[,.(mean=mean(linear_predictor), sd=sd(linear_predictor)), by=.(score_type, model, sex)]
pred_score_means[, LSD := mean-sd]
pred_score_means[, USD := mean+sd]

g <- ggplot(pred_scores[score_type == "non-derived"]) +
  aes(x=linear_predictor, color=model) +
  facet_wrap(~ sex) +
  geom_density(trim=TRUE) +
  geom_vline(data=pred_score_means[score_type == "non-derived"], aes(xintercept=mean, color=model), linetype=2) +
  geom_rect(data=pred_score_means[score_type == "non-derived"], inherit.aes=FALSE, aes(ymin=-Inf, ymax=Inf, xmin=LSD, xmax=USD, fill=model), alpha=0.2) +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  scale_fill_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  xlab("Linear predictor") +
  ylab("Density") +
  theme_bw() +
  theme(
    axis.title=element_text(size=8), axis.text=element_text(size=6), 
    legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=7), 
    strip.background=element_blank(), strip.text=element_text(size=7, face="bold")
  )
ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/linear_predictor_densities_after_scaling.pdf")

# Add five year age group (downstream analyses uses these)
pred_scores[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Compute absolute risk using SCORE2 baseline hazards
pred_scores[, uncalibrated_risk := score2_absrisk(sex, linear_predictor)]

# Compute absolute risk recalibrated to low-risk european region (including UK)
pred_scores[, uk_calibrated_risk := score2_recalibration(sex, uncalibrated_risk, "low")]

# Reorganize columns
pred_scores <- pred_scores[, .(eid, sex, age, age_group, incident_cvd, incident_cvd_followup, cvd_prediction_foldid=test_fold,
  model, score_type, linear_predictor, uncalibrated_risk, uk_calibrated_risk)]

# Write out
fwrite(pred_scores, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Compare absolute risks
pred_score_means <- pred_scores[,.(mean=mean(uk_calibrated_risk), sd=sd(uk_calibrated_risk)), by=.(score_type, model, sex)]
pred_score_means[, LSD := mean-sd]
pred_score_means[, USD := mean+sd]

g <- ggplot(pred_scores[score_type == "non-derived"]) +
  aes(x=uk_calibrated_risk, color=model) +
  facet_wrap(~ sex) +
  geom_density(trim=TRUE) +
  geom_vline(data=pred_score_means[score_type == "non-derived"], aes(xintercept=mean, color=model), linetype=2) +
  geom_rect(data=pred_score_means[score_type == "non-derived"], inherit.aes=FALSE, aes(ymin=-Inf, ymax=Inf, xmin=LSD, xmax=USD, fill=model), alpha=0.2) +
  scale_color_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  scale_fill_manual(values=c("SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRSs"="#377eb8", "SCORE2 + NMR scores + PRSs"="#4daf4a")) +
  xlab("Calibrated 10-year CVD risk") +
  ylab("Density") +
  theme_bw() +
  theme(
    axis.title=element_text(size=8), axis.text=element_text(size=6), 
    legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=7), 
    strip.background=element_blank(), strip.text=element_text(size=7, face="bold")
  )
ggsave(g, width=7.2, height=5, file="analyses/CVD_weight_training/uk_calibrated_risk_densities_after_scaling.pdf")

# Build supp table of avg weights and transformations
dt <- cvd_weights[score_type == "non-derived", .(weight=mean(weight)), by=.(model, sex, score)]
dt <- dcast(dt, model + sex ~ score, value.var="weight")
dt2 <- scaling_factors[score_type == "non-derived", .(scale=mean(scale), offset=mean(offset)), by=.(model, sex)]
dt <- dt[dt2, on = .(model, sex)]
dt <- dt[order(-sex)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/avg_cvd_score_weights_supp_table.txt")

