library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(cowplot)
source("src/utils/SCORE2.R")
source('src/utils/mean_pvalue.R')
source('src/utils/format_pval.R')

# Make output directory
system("mkdir -p analyses/CVD_weight_training", wait=TRUE)

# Load required data
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "age", "incident_cvd_followup", "incident_cvd", "cvd_prediction_foldid", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2_excl_UKB"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

# Add in test scores
test_scores <- fread("analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")
setnames(test_scores, gsub("NMR", "non_derived_NMR", names(test_scores)))
dat <- dat[test_scores, on = .(eid)]

test_scores <- fread("analyses/nmr_score_training/aggregate_test_clinical_NMR_scores.txt")
setnames(test_scores, gsub("NMR", "clinical_NMR", names(test_scores)))
dat <- dat[test_scores, on = .(eid)]

# Standardise PRSs
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]

# Estimate weights using Cox proportional hazards models
cvd_weights <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        # extract training data
        this_dat <- dat[cvd_prediction_foldid != this_test_fold & sex == this_sex]
  
        # Build model formula
        mf <- "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2)"
        if (this_model %in% c("SCORE2 + NMR scores", "SCORE2 + NMR scores + PRSs")) {
          if (this_score_type == "non-derived") {
            mf <- paste(mf, "+ CAD_non_derived_NMR_score + Stroke_non_derived_NMR_score")
          } else {
            mf <- paste(mf, "+ CAD_clinical_NMR_score + Stroke_clinical_NMR_score")
          }
        }
        if (this_model %in% c("SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")) {
          mf <- paste(mf, "+ CAD_metaGRS + Stroke_metaGRS")
        }
   
        # Fit cox model
        cx <- coxph(as.formula(mf), data=this_dat)

        # Extract relevant coefficients
        ci <- confint(cx)
        cf <- as.data.table(coef(summary(cx)), keep.rownames="score")
        cf[, score := gsub("_non_derived", "", score)]
        cf[, score := gsub("_clinical", "", score)]
        cf[, score := gsub("_", " ", score)]
        cf <- cf[, .(score, weight=coef, L95=ci[,1], U95=ci[,2], pval=`Pr(>|z|)`)]

        # add information
        cbind(test_fold = this_test_fold, sex = this_sex, model = this_model, score_type = this_score_type, cf)
}

# Compute average weights
avg_cvd_weights <- cvd_weights[, .(
  weight=round(mean(weight), digits=3),
  L95=round(mean(L95), digits=3),
  U95=round(mean(U95), digits=3),
  pval=format_pval(mean_pvalue(pval, weight))
), by=.(sex, model, score_type, score)]

# Write out
fwrite(cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/cvd_score_weights_in_5fold_cross_validation.txt")
fwrite(avg_cvd_weights, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/average_cvd_score_weights.txt")

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

# Extract linear predictors in each test fold for downstream analyses
pred_scores <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        # extract test data
        this_dat <- dat[cvd_prediction_foldid == this_test_fold & sex == this_sex]

        # Melt columns of interest to long format for summing
        if (this_score_type == "non-derived") {
          setnames(this_dat, gsub("_non_derived", "", names(this_dat)))
        } else {
          setnames(this_dat, gsub("_clinical", "", names(this_dat)))
        }
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
      
        # Sum to obtain linear predictor
        this_dat <- this_dat[, .(linear_predictor=sum(value)), by=.(eid, sex, age, incident_cvd, incident_cvd_followup, cvd_prediction_foldid)]
   
        # add model info
        cbind(score_type = this_score_type, model = this_model, this_dat)
}

# Add five year age group (downstream analyses uses these)
pred_scores[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]

# Compute absolute risk using SCORE2 baseline hazards
pred_scores[, uncalibrated_risk := score2_absrisk(sex, linear_predictor)]

# Compute absolute risk recalibrated to low-risk european region (including UK)
pred_scores[, uk_calibrated_risk := score2_recalibration(sex, uncalibrated_risk, "low")]

# Reorganize columns
pred_scores <- pred_scores[, .(eid, sex, age, age_group, incident_cvd, incident_cvd_followup, cvd_prediction_foldid,
  model, score_type, linear_predictor, uncalibrated_risk, uk_calibrated_risk)]

# Write out
fwrite(pred_scores, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Format and output representative weights for supp table
dt <- cvd_weights[score_type == "non-derived", .(weight=mean(weight)), by=.(model, sex, score)]
dt[, sex := factor(sex, levels=c("Male", "Female"))]
dt[, score := factor(score, levels=c("CAD NMR score", "Stroke NMR score", "CAD metaGRS", "Stroke metaGRS"))]
dt[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dt <- dcast(dt, model + score ~ sex, value.var="weight")
prs_scaling <- fread("data/standardised/prs_scaling_factors.txt")
prs_scaling[, PRS := gsub("_", " ", PRS)]
dt[prs_scaling, on = .(score=PRS), c("mean", "sd") := .(i.mean, i.sd)]
dt <- dt[,.(model, score, mean, sd, Male, Female)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/CVD_weight_training/avg_cvd_score_weights_supp_table.txt")

# Compute weighted versions of NMR scores, so that we can get estimates of HRs per SD
weighted_NMR <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %:%
      foreach(this_test_fold = 1:5, .combine=rbind) %do% {
        # extract test data
        this_dat <- dat[cvd_prediction_foldid == this_test_fold & sex == this_sex]

        # Melt columns of interest to long format for summing
        if (this_score_type == "non-derived") {
          setnames(this_dat, gsub("_non_derived", "", names(this_dat)))
        } else {
          setnames(this_dat, gsub("_clinical", "", names(this_dat)))
        }
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

        # Extract NMR scores 
        this_dat <- dcast(this_dat[score %like% "NMR"], eid + sex  ~ score, value.var="value")
      
        # add model info
        cbind(score_type = this_score_type, model = this_model, this_dat)
}

# Estimate independent contributions of NMR scores (and PRS) to CVD prediction
dat <- fread("data/cleaned/analysis_cohort.txt")
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]
weighted_NMR[, CAD_NMR_score := scale(CAD_NMR_score), by=.(sex, model, score_type)]
weighted_NMR[, Stroke_NMR_score := scale(Stroke_NMR_score), by=.(sex, model, score_type)]

cvd_hrs <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
        # Extract the subset of the data to run the analysis in
        this_dat <- dat[sex == this_sex]
        if (this_model %like% "NMR") {
          this_scores <- weighted_NMR[model == this_model & sex == this_sex & score_type == this_score_type, .(eid, CAD_NMR_score, Stroke_NMR_score)]
          this_dat <- this_dat[this_scores, on = .(eid)]
        }

        # Build model formula
        mf <- "Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2_excl_UKB)"
        if (this_model %like% "NMR") {
          mf <- paste(mf, "+ CAD_NMR_score + Stroke_NMR_score")
        }
        if (this_model %like% "PRS") {
          mf <- paste(mf, "+ CAD_metaGRS + Stroke_metaGRS")
        }
   
        # Fit cox model
        cx <- coxph(as.formula(mf), data=this_dat)

        # Extract relevant coefficients
        ci <- confint(cx)
        cf <- as.data.table(coef(summary(cx)), keep.rownames="score")
        cf[, score := gsub("_", " ", score)]
        cf <- cf[, .(score, HR=`exp(coef)`, L95=exp(ci[,1]), U95=exp(ci[,2]), pval=`Pr(>|z|)`)]

        # add information
        cbind(sex = this_sex, model = this_model, score_type = this_score_type, cf)
}

# Write out
fwrite(cvd_hrs, sep="\t", quote=FALSE, file="analyses/test/score_HRs_independent_of_score2.txt")

# Factor levels for plot ordering
cvd_hrs[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
cvd_hrs[, score := factor(score, levels=c("CAD NMR score", "Stroke NMR score", "CAD metaGRS", "Stroke metaGRS"))]

# Plot HRs
plot_hrs <- function(this_sex, this_score_type) {
  ggplot(cvd_hrs[sex == this_sex & score_type == this_score_type]) +
    aes(x = score, y = HR, ymin = L95, ymax = U95) +
    facet_grid(. ~ model, space="free_x", scales="free_x") +
    geom_errorbar(width=0) +
    geom_point(shape=23, fill="white") +
    geom_hline(yintercept=1, linetype=2) +
    xlab("") + ylab("Hazard Ratio (95% CI)") +
    guides(color=guide_legend(title="")) +
    theme_bw() +
    theme(
      axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
      axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
      strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
      panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()
    )
}

g1 <- plot_hrs("Male", "non-derived")
g2 <- plot_hrs("Female", "non-derived")
g <- plot_grid(g1, g2, nrow=2, labels=c("Males  ", "Females"), label_size=10)
ggsave(g, width=7.2, height=4.5, file="analyses/test/non_derived_score_HRs_independent_of_score2.pdf")

g1 <- plot_hrs("Male", "clinical")
g2 <- plot_hrs("Female", "clinical")
g <- plot_grid(g1, g2, nrow=2, labels=c("Males  ", "Females"), label_size=10)
ggsave(g, width=7.2, height=4.5, file="analyses/test/clinical_score_HRs_independent_of_score2.pdf")
  
# Now estimate HRs independently of individual SCORE2 risk factors
dat[, age := (age - 60)/5]
dat[is.na(smoking), smoking := FALSE]
dat[, smoking := factor(smoking, levels=c(FALSE, TRUE))]
dat[, sbp := (sbp - 120)/20]
dat[, tchol := (tchol - 6)/1]
dat[, hdl := (hdl - 1.3)/0.5]

cvd_hrs2 <- foreach(this_score_type = c("non-derived", "clinical"), .combine=rbind) %:%
  foreach(this_model = c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:%
    foreach(this_sex = c("Male", "Female"), .combine=rbind) %do% {
        # Extract the subset of the data to run the analysis in
        this_dat <- dat[sex == this_sex]
        if (this_model %like% "NMR") {
          this_scores <- weighted_NMR[model == this_model & sex == this_sex & score_type == this_score_type, .(eid, CAD_NMR_score, Stroke_NMR_score)]
          this_dat <- this_dat[this_scores, on = .(eid)]
        }

        # Build model formula
        mf <- "Surv(incident_cvd_followup, incident_cvd) ~ age*smoking + age*sbp + age*tchol + age*hdl"
        if (this_model %like% "NMR") {
          mf <- paste(mf, "+ CAD_NMR_score + Stroke_NMR_score")
        }
        if (this_model %like% "PRS") {
          mf <- paste(mf, "+ CAD_metaGRS + Stroke_metaGRS")
        }
   
        # Fit cox model
        cx <- coxph(as.formula(mf), data=this_dat)

        # Extract relevant coefficients
        ci <- confint(cx)
        cf <- as.data.table(coef(summary(cx)), keep.rownames="score")
        cf[, score := gsub("_", " ", score)]
        cf <- cf[, .(score, HR=`exp(coef)`, L95=exp(ci[,1]), U95=exp(ci[,2]), pval=`Pr(>|z|)`)]

        # add information
        cbind(sex = this_sex, model = this_model, score_type = this_score_type, cf)
}
setnames(cvd_hrs2, "score", "coefficient")
cvd_hrs2[, coefficient := fcase(
  coefficient == "age", "Age (per 5-year increase)",
  coefficient == "smokingTRUE", "Smoking (current vs. other)",
  coefficient == "sbp", "SBP (per 20 mm Hg increase)",
  coefficient == "tchol", "Total cholesterol (per 1.0 mmol/L increase)",
  coefficient == "hdl", "HDL cholesterol (per 0.5 mmol/L increase)",
  coefficient == "age:smokingTRUE", "Age x smoking",
  coefficient == "age:sbp", "Age x SBP",
  coefficient == "age:tchol", "Age x total cholesterol",
  coefficient == "age:hdl", "Age x HDL cholesterol",
  coefficient == "CAD NMR score", "CAD NMR score (per SD increase)",
  coefficient == "Stroke NMR score", "Stroke NMR score (per SD increase)",
  coefficient == "CAD metaGRS", "CAD metaGRS (per SD increase)",
  coefficient == "Stroke metaGRS", "Stroke metaGRS (per SD increase)"
)]

# Write out
fwrite(cvd_hrs2, sep="\t", quote=FALSE, file="analyses/test/score_HRs_independent_of_score2_components.pdf")

# Factor levels for plot ordering
cvd_hrs2[, model := factor(model, levels=c("SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
cvd_hrs2[, coefficient := factor(coefficient, levels=c(
  "Age (per 5-year increase)", "Smoking (current vs. other)", "SBP (per 20 mm Hg increase)", 
  "Total cholesterol (per 1.0 mmol/L increase)", "HDL cholesterol (per 0.5 mmol/L increase)", 
  "Age x smoking", "Age x SBP", "Age x total cholesterol",  "Age x HDL cholesterol", 
  "CAD NMR score (per SD increase)", "Stroke NMR score (per SD increase)", "CAD metaGRS (per SD increase)",
  "Stroke metaGRS (per SD increase)"))]
cvd_hrs2[, color_anno := fcase(
  coefficient %like% "NMR", "NMR scores",
  coefficient %like% "GRS", "PRSs",
  default = "Component of SCORE2"
)]
cvd_hrs2[, color_anno := factor(color_anno, levels=c("Component of SCORE2", "NMR scores", "PRSs"))]
cvd_hrs2[, sex := factor(paste0(sex, "s"), levels=c("Males", "Females"))]

# Plot HRs
plot_hrs <- function(this_score_type) {
  ggplot(cvd_hrs2[score_type == this_score_type & model == "SCORE2 + NMR scores + PRSs"]) +
    aes(x = coefficient, y = HR, ymin = L95, ymax = U95, color=color_anno) +
    facet_grid(. ~ sex, space="free_x", scales="free_x") +
    geom_hline(yintercept=1, linetype=2) +
    geom_errorbar(width=0) +
    geom_point(shape=23, fill="white") +
    scale_color_manual(values=c("Component of SCORE2"="black", "NMR scores"="#762a83", "PRSs"="#1b7837")) +
    xlab("") + ylab("Hazard Ratio (95% CI)") +
    guides(color=guide_legend(title="")) +
    theme_bw() +
    theme(
      axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
      axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
      strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
      panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
      legend.text=element_text(size=6), legend.position="bottom"
    )
}

g <- plot_hrs("non-derived")
ggsave(g, width=7.2, height=4.5, file="analyses/test/non_derived_score_HRs_independent_of_score2_components.pdf")

g <- plot_hrs("clinical")
ggsave(g, width=7.2, height=4.5, file="analyses/test/clinical_score_HRs_independent_of_score2_components.pdf")

# Make supp figure
ggdt <- rbind(idcol="model",
  "Risk factors"=cvd_hrs2[score_type == "non-derived" & model == "SCORE2 + NMR scores + PRSs", .(sex, coefficient, HR, L95, U95, pval)],
  "SCORE2"=cvd_hrs[score_type == "non-derived" & model == "SCORE2 + NMR scores + PRSs", .(sex=paste0(sex, "s"), coefficient=score, HR, L95, U95, pval)]
)
ggdt[model == "SCORE2", coefficient := paste(coefficient, "(per SD increase)")]
ggdt[, coefficient := factor(coefficient, levels=c(
  "Age (per 5-year increase)", "Smoking (current vs. other)", "SBP (per 20 mm Hg increase)",
  "Total cholesterol (per 1.0 mmol/L increase)", "HDL cholesterol (per 0.5 mmol/L increase)",
  "Age x smoking", "Age x SBP", "Age x total cholesterol",  "Age x HDL cholesterol",
  "CAD NMR score (per SD increase)", "Stroke NMR score (per SD increase)", "CAD metaGRS (per SD increase)",
  "Stroke metaGRS (per SD increase)"))]
ggdt[, sex := factor(sex, levels=c("Males", "Females"))]
ggdt[, model := factor(model, levels=c("SCORE2", "Risk factors"))]
ggdt[, color_anno := fcase(
  coefficient %like% "NMR", "NMR scores",
  coefficient %like% "GRS", "PRSs",
  default = "Component of SCORE2"
)]
ggdt[, color_anno := factor(color_anno, levels=c("Component of SCORE2", "NMR scores", "PRSs"))]

g1 <- ggplot(ggdt[model == "SCORE2"]) +
  aes(x = coefficient, y = HR, ymin = L95, ymax = U95, color=color_anno) +
  facet_grid(sex ~ .) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("Component of SCORE2"="black", "NMR scores"="#762a83", "PRSs"="#1b7837")) +
  scale_y_continuous("Hazard Ratio (95% CI)") +
  guides(color=guide_legend(title="")) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
    axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5), axis.title.x=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.position="bottom"
  )

g2 <- ggplot(ggdt[model == "Risk factors"]) +
  aes(x = coefficient, y = HR, ymin = L95, ymax = U95, color=color_anno) +
  facet_grid(sex ~ .) +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0) +
  geom_point(shape=23, fill="white") +
  scale_color_manual(values=c("Component of SCORE2"="black", "NMR scores"="#762a83", "PRSs"="#1b7837")) +
  scale_y_continuous("Hazard Ratio (95% CI)") +
  guides(color=guide_legend(title="")) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8),
    axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5), axis.title.x=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.position="bottom"
  )

g <- plot_grid(g1, g2, nrow=1, rel_widths=c(1, 3), align="hv")
ggsave(g, width=7.2, height=5, file="analyses/test/multivariable_HRs_NMR_scores_and_PRSs.pdf")


# Format table for supp
dt <- dcast(ggdt, model + coefficient ~ sex, value.var=c("HR", "L95", "U95", "pval"))
dt <- dt[, .(model, coefficient, HR_Males, L95_Males, U95_Males, pval_Males, HR_Females, L95_Females, U95_Females, pval_Females)]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/multivariable_HRs_NMR_scores_and_PRSs.txt")




