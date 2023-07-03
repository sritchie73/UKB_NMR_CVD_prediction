library(data.table)
library(foreach)
library(ukbnmr)
library(ggplot2)
library(forcats)
library(ggrastr)
library(ggpp)
library(ggh4x)
library(hexbin)
library(cowplot)
library(survival)
source("src/utils/SCORE2.R")
source("src/utils/cox_test.R")
source("src/utils/score_cindex.R")

options("ggrastr.default.dpi" = 300) 

# Load dataset with NMR biomarkers imputed
dat <- fread("data/imputed/analysis_cohort.txt")

# Load scaling factors for NMR data
nmr_scaling <- fread("data/standardised/nmr_scaling_factors.txt")

# Load NMR score coefficients
nmr_coef <- fread("analyses/nmr_score_training/best_fits_coef.txt")

# Build long table of NMR coefficients and age interactions standardised using the given scaling factors
nmr_long <- melt(dat, id.vars=c("eid", "sex", "age", "prediction_cv_foldid"), measure.vars=ukbnmr::nmr_info[Type == "Non-derived", Biomarker])
nmr_long[nmr_scaling, on = .(variable=biomarker), value := (value - i.mean)/i.sd]
nmr_long[, age := (age - 60)/5]
nmr_age_interactions <- copy(nmr_long)
nmr_age_interactions[, value := value * age]
nmr_age_interactions[, variable := paste0("age:", variable)]
nmr_long <- rbind(nmr_long, nmr_age_interactions)

# Compute NMR biomarker scores
models <- unique(nmr_coef[,.(prediction_cv_testfold, endpoint, lambda.fit)])
nmr <- foreach(midx = models[,.I], .combine=rbind) %do% {
  this_info <- cbind(score="NMR", models[midx])
  this_coef <- nmr_coef[this_info, on = .(prediction_cv_testfold, endpoint, lambda.fit)]
  this_score <- nmr_long[this_coef, on = .(sex, variable=coef)][, .(linear_predictor=sum(beta * value)), by=.(eid, sex, prediction_cv_foldid)]
  cbind(this_info, this_score)
}

# Compute average score across all 5 training runs
avg_nmr_coef <- fread("analyses/nmr_score_training/avg_best_coefs.txt")
models <- unique(models[,.(endpoint, lambda.fit)])
avg_nmr <- foreach(midx = models[,.I], .combine=rbind) %do% {
  this_info <- cbind(score="NMR", models[midx])
  this_coef <- avg_nmr_coef[this_info, on = .(endpoint, lambda.fit)]
  this_score <- nmr_long[this_coef, on = .(sex, variable=coef)][, .(linear_predictor=sum(beta * value)), by=.(eid, sex, prediction_cv_foldid)]
  cbind(this_info, this_score)
}

# Assess consistency of NMR biomarker scores across training runs
score_comp <- foreach(this_test_fold = c(1:5, "avg"), .combine=rbind) %do% {
  if (this_test_fold == "avg") {
    this_x <- avg_nmr[, .(eid, prediction_cv_foldid, sex, endpoint, lambda.fit,
        x_model="avg", x_score=linear_predictor)]
  } else {
    this_x <- nmr[prediction_cv_testfold == this_test_fold,  
      .(eid, prediction_cv_foldid, sex, endpoint, lambda.fit, 
        x_model=prediction_cv_testfold, x_score=linear_predictor)]
  } 
  this_y <- nmr[, .(eid, prediction_cv_foldid, sex, endpoint, lambda.fit,
      y_model=prediction_cv_testfold, y_score=linear_predictor)]
  this_y <- rbind(this_y, avg_nmr[,
    .(eid, prediction_cv_foldid, sex, endpoint, lambda.fit,
      y_model="avg", y_score=linear_predictor)])
  merge(this_x, this_y, by=c("eid", "prediction_cv_foldid", "sex", "endpoint", "lambda.fit"), all=FALSE)
}

# Color annotation for plots
score_comp[, color_anno := fcase(
  x_model != prediction_cv_foldid & y_model != prediction_cv_foldid, "Training samples",
  x_model == prediction_cv_foldid, "X-axis test samples", 
  y_model == prediction_cv_foldid, "Y-axis test samples" 
)]

# Create scatterplots of scores
score_scatter <- function(score_comp) {
  # Better facet names
  facet_anno <- function(vars) { 
    sapply(vars, function(ii) {
      if (ii == "avg") {
        "Average"
      } else {
        paste("Test fold", ii) 
      }
    })
  }

  # correlation annoations
  cor_anno <- score_comp[,.(label=sprintf("r=%.2f\np=%.2f", 
    cor(x_score, y_score), 
    cor(x_score, y_score, method="spearman")
  )), by=.(x_model, y_model)]

  ggplot(score_comp) +
    aes(x=x_score, y=y_score, color=color_anno) +
    rasterise(geom_point(data=score_comp[x_model != prediction_cv_foldid & y_model != prediction_cv_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model == prediction_cv_foldid & y_model != prediction_cv_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model != prediction_cv_foldid & y_model == prediction_cv_foldid], alpha=0.5, shape=19, size=0.1)) +
    rasterise(geom_point(data=score_comp[x_model == prediction_cv_foldid & y_model == prediction_cv_foldid], alpha=0.5, shape=19, size=0.1)) +
    scale_color_manual(name="", values=c("Training samples"="#252525", "X-axis test samples"="#d94801", "Y-axis test samples"="#006d2c")) +
    guides(color = guide_legend(override.aes = list(size = 1))) +
    geom_abline(intercept=0, slope=1, linetype=2, color="red") +
    geom_text_npc(data=cor_anno, npcx="left", npcy="top", aes(label=label), color="red", size=6*0.352777778) +
    facet_grid2(fct_rev(y_model) ~ x_model, scales="free", independent="all", labeller=as_labeller(facet_anno)) +
    xlab("NMR score") + ylab("NMR score") +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(), 
      panel.grid=element_blank(), legend.text=element_text(size=6)
    )
}
g1 <- score_scatter(score_comp[endpoint == "CAD" & sex == "Male" & lambda.fit == "lambda.min"])
g2 <- score_scatter(score_comp[endpoint == "CAD" & sex == "Female" & lambda.fit == "lambda.min"])
g3 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Male" & lambda.fit == "lambda.min"])
g4 <- score_scatter(score_comp[endpoint == "Stroke" & sex == "Female" & lambda.fit == "lambda.min"])

ggsave(g1, width=7, height=6, file="analyses/nmr_score_training/CAD_male_score_consistency.pdf")
ggsave(g2, width=7, height=6, file="analyses/nmr_score_training/CAD_female_score_consistency.pdf")
ggsave(g3, width=7, height=6, file="analyses/nmr_score_training/Stroke_male_score_consistency.pdf")
ggsave(g4, width=7, height=6, file="analyses/nmr_score_training/Stroke_female_score_consistency.pdf")

# Compare average score to test score
avg_vs_test_scatter <- function(lambda) {
  cor_anno <- score_comp[x_model == "avg" & y_model == prediction_cv_foldid & lambda.fit == lambda,
    .(label=sprintf("r=%.2f\np=%.2f",
      cor(x_score, y_score),
      cor(x_score, y_score, method="spearman")
    )), by=.(endpoint, sex)]

  ggplot(score_comp[x_model == "avg" & y_model == prediction_cv_foldid & lambda.fit == lambda]) +
    aes(x=x_score, y=y_score, color=paste("Test fold", y_model)) +
    rasterise(geom_point(alpha=0.5, shape=19, size=0.1)) +
    guides(color = guide_legend(title="", override.aes = list(size = 1))) +
    geom_abline(intercept=0, slope=1, linetype=2, color="red") +
    geom_text_npc(data=cor_anno, npcx="left", npcy="top", aes(label=label), color="red", size=6*0.352777778) +
    facet_grid2(fct_rev(sex) ~ endpoint, scales="free", independent="all") +
    xlab("NMR score from coefficient averages") + ylab("NMR scores aggregated across test folds") +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
      panel.grid=element_blank(), legend.text=element_text(size=6), legend.title=element_text(size=8)
    )
}

g <- avg_vs_test_scatter("lambda.min")
ggsave(g, width=7, height=5.5, file="analyses/nmr_score_training/aggregate_test_vs_coef_average.pdf")

# Write out aggregate scores
test_scores <- nmr[lambda.fit == "lambda.min" & prediction_cv_testfold == prediction_cv_foldid]
test_scores <- dcast(test_scores, eid  ~ endpoint, value.var="linear_predictor")
setnames(test_scores, c("CAD", "Stroke"), c("CAD_NMR_score", "Stroke_NMR_score"))
fwrite(test_scores, sep="\t", quote=FALSE, file="analyses/nmr_score_training/aggregate_test_NMR_scores.txt")

# Sanity check hazard ratios/C-index in conventional Cox models
if ("CAD_NMR_score" %in% names(dat)) dat[, c("CAD_NMR_score", "Stroke_NMR_score") := NULL]
dat <- dat[test_scores, on = .(eid)]

cx_list <- list(
  "NMR_CAD_males"=cox.test("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_NMR_score", "incident_cad", dat[sex == "Male"]),
  "NMR_CAD_females"=cox.test("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_NMR_score", "incident_cad", dat[sex == "Female"]),
  "NMR_CAD_strata"=cox.test("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2 + CAD_NMR_score", "incident_cad", dat),
  "NMR_Stroke_males"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_NMR_score", "incident_stroke", dat[sex == "Male"]),
  "NMR_Stroke_females"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_NMR_score", "incident_stroke", dat[sex == "Female"]),
  "NMR_Stroke_strata"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2 + Stroke_NMR_score", "incident_stroke", dat)
)

hrs <- rbindlist(idcol="model", lapply(cx_list, `[[`, "coefficients"))
hrs[, c("score", "endpoint", "sex") := tstrsplit(model, "_")]

cinds <- rbind(idcol="model",
  "SCORE2_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2", dat[sex == "Male"]),
  "SCORE2_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2", dat[sex == "Female"]),
  "SCORE2_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2", dat),
  "SCORE2_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2", dat[sex == "Male"]),
  "SCORE2_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2", dat[sex == "Female"]),
  "SCORE2_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2", dat),
  "NMR_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_NMR_score", dat[sex == "Male"]),
  "NMR_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_NMR_score", dat[sex == "Female"]),
  "NMR_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2 + CAD_NMR_score", dat),
  "NMR_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_NMR_score", dat[sex == "Male"]),
  "NMR_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_NMR_score", dat[sex == "Female"]),
  "NMR_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2 + Stroke_NMR_score", dat)
)
cinds[, c("score", "endpoint", "sex") := tstrsplit(model, "_")]

avg_scores <- dcast(avg_nmr[lambda.fit == "lambda.min"], eid ~ endpoint, value.var="linear_predictor")
setnames(avg_scores, c("CAD", "Stroke"), c("CAD_NMR_score", "Stroke_NMR_score"))

dat[, c("CAD_NMR_score", "Stroke_NMR_score") := NULL]
dat <- dat[avg_scores, on = .(eid)]

cx_list2 <- list(
  "NMR_CAD_males"=cox.test("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_NMR_score", "incident_cad", dat[sex == "Male"]),
  "NMR_CAD_females"=cox.test("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_NMR_score", "incident_cad", dat[sex == "Female"]),
  "NMR_CAD_strata"=cox.test("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2 + CAD_NMR_score", "incident_cad", dat),
  "NMR_Stroke_males"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_NMR_score", "incident_stroke", dat[sex == "Male"]),
  "NMR_Stroke_females"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_NMR_score", "incident_stroke", dat[sex == "Female"]),
  "NMR_Stroke_strata"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2 + Stroke_NMR_score", "incident_stroke", dat)
)

hrs2 <- rbindlist(idcol="model", lapply(cx_list2, `[[`, "coefficients"))
hrs2[, c("score", "endpoint", "sex") := tstrsplit(model, "_")]

cinds2 <- rbind(idcol="model",
  "SCORE2_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2", dat[sex == "Male"]),
  "SCORE2_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2", dat[sex == "Female"]),
  "SCORE2_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2", dat),
  "SCORE2_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2", dat[sex == "Male"]),
  "SCORE2_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2", dat[sex == "Female"]),
  "SCORE2_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2", dat),
  "NMR_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_NMR_score", dat[sex == "Male"]),
  "NMR_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_NMR_score", dat[sex == "Female"]),
  "NMR_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2 + CAD_NMR_score", dat),
  "NMR_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_NMR_score", dat[sex == "Male"]),
  "NMR_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_NMR_score", dat[sex == "Female"]),
  "NMR_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2 + Stroke_NMR_score", dat)
)
cinds2[, c("score", "endpoint", "sex") := tstrsplit(model, "_")]

hrs <- rbind(idcol="type", "aggregate"=hrs, "average"=hrs2)
cinds <- rbind(idcol="type", "aggregate"=cinds, "average"=cinds2)

hrs[, coefficient := factor(coefficient, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score"))]
hrs[, sex := factor(sex, levels=c("males", "females", "strata"))]

cinds[score == "NMR", score := sprintf("SCORE2 +\n%s_NMR_score", endpoint)]
cinds[, score := factor(score, levels=c("SCORE2", "SCORE2 +\nCAD_NMR_score", "SCORE2 +\nStroke_NMR_score"))]
cinds[, sex := factor(sex, levels=c("males", "females", "strata"))]

g1 <- ggplot(hrs) +
  aes(x=coefficient, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_wrap(~ type + endpoint, nrow=1, scales="free") +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) + 
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("") + ylab("Hazard Ratio (95% CI) per SD") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), axis.text.x=element_text(size=6, angle=90, hjust=1),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g1, width=7, height=3, file="analyses/nmr_score_training/hazard_ratio_sanity_check.pdf")

g2 <- ggplot(cinds) +
  aes(x=score, y=C.index, ymin=L95, ymax=U95, color=sex) + 
  facet_wrap(~ type + endpoint, nrow=1, scales="free") +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) +
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("") + ylab("C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), axis.text.x=element_text(size=6, angle=90, hjust=1),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g2, width=7, height=3, file="analyses/nmr_score_training/cindex_sanity_check.pdf")

score2_cinds <- cinds[score == "SCORE2"]
score2_cinds[, score := sprintf("SCORE2 +\n%s_NMR_score", endpoint)]
cinds[score2_cinds, on = .(type, score, endpoint, sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]

g3 <- ggplot(cinds[!is.na(deltaC)]) +
  aes(x=score, y=deltaC, ymin=deltaC.L95, ymax=deltaC.U95, color=sex) + 
  facet_wrap(~ type + endpoint, nrow=1, scales="free") +
  geom_hline(yintercept=0, linetype=2) +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) +
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("") + ylab("Delta C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), axis.text.x=element_text(size=6, angle=90, hjust=1),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g3, width=7, height=3, file="analyses/nmr_score_training/delta_cindex_sanity_check.pdf")

fwrite(hrs, sep="\t", quote=FALSE, file="analyses/nmr_score_training/hazard_ratio_sanity_check.txt")
cinds[, score := gsub("\\\n", " ", score)]
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/nmr_score_training/cindex_sanity_check.txt")

