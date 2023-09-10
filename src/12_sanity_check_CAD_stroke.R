library(data.table)
library(survival)
library(ggplot2)
source("src/utils/cox_test.R")
source("src/utils/score_cindex.R")

# Load in CAD and Stroke follow-up data and SCORE2
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "incident_cad_followup", "incident_cad", "incident_stroke_followup", "incident_stroke", "SCORE2_excl_UKB"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

# Add in test scores
test_scores <- fread("analyses/nmr_score_training/aggregate_test_non_derived_NMR_scores.txt")
setnames(test_scores, gsub("NMR", "non_derived_NMR", names(test_scores)))
dat <- dat[test_scores, on = .(eid)]

test_scores <- fread("analyses/nmr_score_training/aggregate_test_clinical_NMR_scores.txt")
setnames(test_scores, gsub("NMR", "clinical_NMR", names(test_scores)))
dat <- dat[test_scores, on = .(eid)]

# Compute C-indices from scores directly as linear predictors, (e.g., SCORE2, or SCORE2 + CAD_NMR_score, or SCORE2 + Stroke_NMR_score)
cinds <- rbind(idcol="model",
  "_SCORE2_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2", dat[sex == "Male"]),
  "_SCORE2_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2", dat[sex == "Female"]),
  "_SCORE2_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2", dat),
  "_SCORE2_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2", dat[sex == "Male"]),
  "_SCORE2_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2", dat[sex == "Female"]),
  "_SCORE2_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2", dat),
  "non-derived_NMR_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_non_derived_NMR_score", dat[sex == "Male"]),
  "non-derived_NMR_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_non_derived_NMR_score", dat[sex == "Female"]),
  "non-derived_NMR_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2 + CAD_non_derived_NMR_score", dat),
  "non-derived_NMR_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_non_derived_NMR_score", dat[sex == "Male"]),
  "non-derived_NMR_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_non_derived_NMR_score", dat[sex == "Female"]),
  "non-derived_NMR_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2 + Stroke_non_derived_NMR_score", dat),
  "clinical_NMR_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_clinical_NMR_score", dat[sex == "Male"]),
  "clinical_NMR_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_clinical_NMR_score", dat[sex == "Female"]),
  "clinical_NMR_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2 + CAD_clinical_NMR_score", dat),
  "clinical_NMR_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_clinical_NMR_score", dat[sex == "Male"]),
  "clinical_NMR_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_clinical_NMR_score", dat[sex == "Female"]),
  "clinical_NMR_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2 + Stroke_clinical_NMR_score", dat)
)
cinds[, c("type", "score", "endpoint", "sex") := tstrsplit(model, "_")]
cinds <- cinds[order(sex)][order(endpoint)]

# Compute hazard ratios by fitting Cox proportional hazards models to scores as dependent variables
# These are first standardised across the whole dataset, so that HRs and per SD, and that SD is comparable
# when comparing males to females to stratified models
dat[, SCORE2 := scale(SCORE2)]
dat[, CAD_non_derived_NMR_score := scale(CAD_non_derived_NMR_score)]
dat[, Stroke_non_derived_NMR_score := scale(Stroke_non_derived_NMR_score)]
dat[, CAD_clinical_NMR_score := scale(CAD_clinical_NMR_score)]
dat[, Stroke_clinical_NMR_score := scale(Stroke_clinical_NMR_score)]

cx_list <- list(
  "non-derived_NMR_CAD_males"=cox.test("Surv(incident_cad_followup, incident_cad) ~ offset(SCORE2) + CAD_non_derived_NMR_score", "incident_cad", dat[sex == "Male"]),
  "non-derived_NMR_CAD_females"=cox.test("Surv(incident_cad_followup, incident_cad) ~ offset(SCORE2) + CAD_non_derived_NMR_score", "incident_cad", dat[sex == "Female"]),
  "non-derived_NMR_CAD_strata"=cox.test("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + offset(SCORE2) + CAD_non_derived_NMR_score", "incident_cad", dat),
  "non-derived_NMR_Stroke_males"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ offset(SCORE2) + Stroke_non_derived_NMR_score", "incident_stroke", dat[sex == "Male"]),
  "non-derived_NMR_Stroke_females"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ offset(SCORE2) + Stroke_non_derived_NMR_score", "incident_stroke", dat[sex == "Female"]),
  "non-derived_NMR_Stroke_strata"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + offset(SCORE2) + Stroke_non_derived_NMR_score", "incident_stroke", dat),
  "clinical_NMR_CAD_males"=cox.test("Surv(incident_cad_followup, incident_cad) ~ offset(SCORE2) + CAD_clinical_NMR_score", "incident_cad", dat[sex == "Male"]),
  "clinical_NMR_CAD_females"=cox.test("Surv(incident_cad_followup, incident_cad) ~ offset(SCORE2) + CAD_clinical_NMR_score", "incident_cad", dat[sex == "Female"]),
  "clinical_NMR_CAD_strata"=cox.test("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + offset(SCORE2) + CAD_clinical_NMR_score", "incident_cad", dat),
  "clinical_NMR_Stroke_males"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ offset(SCORE2) + Stroke_clinical_NMR_score", "incident_stroke", dat[sex == "Male"]),
  "clinical_NMR_Stroke_females"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ offset(SCORE2) + Stroke_clinical_NMR_score", "incident_stroke", dat[sex == "Female"]),
  "clinical_NMR_Stroke_strata"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + offset(SCORE2) + Stroke_clinical_NMR_score", "incident_stroke", dat)
)

hrs <- rbindlist(idcol="model", lapply(cx_list, `[[`, "coefficients"))
hrs[, c("type", "score", "endpoint", "sex") := tstrsplit(model, "_")]
hrs <- hrs[order(sex)][order(endpoint)]

# Repeat, using the scores computed from coefficient averages
# (what would be distributed on publication as the representative score, but unsuitable 
# for testing in UK Biobank due to overfitting to training samples)
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "incident_cad_followup", "incident_cad", "incident_stroke_followup", "incident_stroke", "SCORE2_excl_UKB"))
setnames(dat, "SCORE2_excl_UKB", "SCORE2")

test_scores <- fread("analyses/nmr_score_training/coef_avg_non_derived_NMR_scores.txt")
setnames(test_scores, gsub("NMR", "non_derived_NMR", names(test_scores)))
dat <- dat[test_scores, on = .(eid)]

test_scores <- fread("analyses/nmr_score_training/coef_avg_clinical_NMR_scores.txt")
setnames(test_scores, gsub("NMR", "clinical_NMR", names(test_scores)))
dat <- dat[test_scores, on = .(eid)]

# Compute C-indices from scores directly as linear predictors, (e.g., SCORE2, or SCORE2 + CAD_NMR_score, or SCORE2 + Stroke_NMR_score)
cinds2 <- rbind(idcol="model",
  "_SCORE2_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2", dat[sex == "Male"]),
  "_SCORE2_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2", dat[sex == "Female"]),
  "_SCORE2_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2", dat),
  "_SCORE2_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2", dat[sex == "Male"]),
  "_SCORE2_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2", dat[sex == "Female"]),
  "_SCORE2_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2", dat),
  "non-derived_NMR_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_non_derived_NMR_score", dat[sex == "Male"]),
  "non-derived_NMR_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_non_derived_NMR_score", dat[sex == "Female"]),
  "non-derived_NMR_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2 + CAD_non_derived_NMR_score", dat),
  "non-derived_NMR_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_non_derived_NMR_score", dat[sex == "Male"]),
  "non-derived_NMR_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_non_derived_NMR_score", dat[sex == "Female"]),
  "non-derived_NMR_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2 + Stroke_non_derived_NMR_score", dat),
  "clinical_NMR_CAD_males"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_clinical_NMR_score", dat[sex == "Male"]),
  "clinical_NMR_CAD_females"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ SCORE2 + CAD_clinical_NMR_score", dat[sex == "Female"]),
  "clinical_NMR_CAD_strata"=score_cindex("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + SCORE2 + CAD_clinical_NMR_score", dat),
  "clinical_NMR_Stroke_males"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_clinical_NMR_score", dat[sex == "Male"]),
  "clinical_NMR_Stroke_females"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ SCORE2 + Stroke_clinical_NMR_score", dat[sex == "Female"]),
  "clinical_NMR_Stroke_strata"=score_cindex("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + SCORE2 + Stroke_clinical_NMR_score", dat)
)
cinds2[, c("type", "score", "endpoint", "sex") := tstrsplit(model, "_")]
cinds2 <- cinds2[order(sex)][order(endpoint)]

# Compute hazard ratios by fitting Cox proportional hazards models to scores as dependent variables
# These are first standardised across the whole dataset, so that HRs and per SD, and that SD is comparable
# when comparing males to females to stratified models
dat[, CAD_non_derived_NMR_score := scale(CAD_non_derived_NMR_score)]
dat[, Stroke_non_derived_NMR_score := scale(Stroke_non_derived_NMR_score)]
dat[, CAD_clinical_NMR_score := scale(CAD_clinical_NMR_score)]
dat[, Stroke_clinical_NMR_score := scale(Stroke_clinical_NMR_score)]

cx_list2 <- list(
  "non-derived_NMR_CAD_males"=cox.test("Surv(incident_cad_followup, incident_cad) ~ offset(SCORE2) + CAD_non_derived_NMR_score", "incident_cad", dat[sex == "Male"]),
  "non-derived_NMR_CAD_females"=cox.test("Surv(incident_cad_followup, incident_cad) ~ offset(SCORE2) + CAD_non_derived_NMR_score", "incident_cad", dat[sex == "Female"]),
  "non-derived_NMR_CAD_strata"=cox.test("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + offset(SCORE2) + CAD_non_derived_NMR_score", "incident_cad", dat),
  "non-derived_NMR_Stroke_males"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ offset(SCORE2) + Stroke_non_derived_NMR_score", "incident_stroke", dat[sex == "Male"]),
  "non-derived_NMR_Stroke_females"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ offset(SCORE2) + Stroke_non_derived_NMR_score", "incident_stroke", dat[sex == "Female"]),
  "non-derived_NMR_Stroke_strata"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + offset(SCORE2) + Stroke_non_derived_NMR_score", "incident_stroke", dat),
  "clinical_NMR_CAD_males"=cox.test("Surv(incident_cad_followup, incident_cad) ~ offset(SCORE2) + CAD_clinical_NMR_score", "incident_cad", dat[sex == "Male"]),
  "clinical_NMR_CAD_females"=cox.test("Surv(incident_cad_followup, incident_cad) ~ offset(SCORE2) + CAD_clinical_NMR_score", "incident_cad", dat[sex == "Female"]),
  "clinical_NMR_CAD_strata"=cox.test("Surv(incident_cad_followup, incident_cad) ~ strata(sex) + offset(SCORE2) + CAD_clinical_NMR_score", "incident_cad", dat),
  "clinical_NMR_Stroke_males"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ offset(SCORE2) + Stroke_clinical_NMR_score", "incident_stroke", dat[sex == "Male"]),
  "clinical_NMR_Stroke_females"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ offset(SCORE2) + Stroke_clinical_NMR_score", "incident_stroke", dat[sex == "Female"]),
  "clinical_NMR_Stroke_strata"=cox.test("Surv(incident_stroke_followup, incident_stroke) ~ strata(sex) + offset(SCORE2) + Stroke_clinical_NMR_score", "incident_stroke", dat)
)

hrs2 <- rbindlist(idcol="model", lapply(cx_list2, `[[`, "coefficients"))
hrs2[, c("type", "score", "endpoint", "sex") := tstrsplit(model, "_")]
hrs2 <- hrs2[order(sex)][order(endpoint)]

# Combine results
cinds <- rbind(idcol="type2", "aggregate"=cinds, "average"=cinds2)
hrs <- rbind(idcol="type2", "aggregate"=hrs, "average"=hrs2)

# Show and compare hazard ratios
g <- ggplot(hrs) +
  aes(x=coefficient, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_wrap(~ endpoint + type2 + type, nrow=1, scales="free_x") +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) + 
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("NMR score") + ylab("Hazard Ratio (95% CI) per SD") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), 
    axis.text.x=element_blank(), axis.title.x=element_text(size=8), axis.ticks.x=element_blank(),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g, width=7, height=3, file="analyses/nmr_score_training/hazard_ratio_sanity_check.pdf")

# Show change in C-index relative to SCORE2
s1 <- cinds[score == "SCORE2"]
s2 <- cinds[score == "SCORE2"]
s1[, model := gsub("_SCORE2_", "non-derived_NMR_", model)]
s2[, model := gsub("_SCORE2_", "clinical_NMR_", model)]
score2_cinds <- rbind(s1, s2)
cinds[score2_cinds, on = .(type2, model), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]

g <- ggplot(cinds[!is.na(deltaC)]) +
  aes(x=score, y=deltaC, ymin=deltaC.L95, ymax=deltaC.U95, color=sex) + 
  facet_wrap(~ endpoint + type2 + type, nrow=1, scales="free_x") +
  geom_hline(yintercept=0, linetype=2) +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) +
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("NMR score") + ylab("Delta C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), 
    axis.text.x=element_blank(), axis.title.x=element_text(size=8), axis.ticks.x=element_blank(),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g, width=7, height=3, file="analyses/nmr_score_training/delta_cindex_sanity_check.pdf")

# Write out tables
hrs[, model := NULL]
cinds[, model := NULL]
fwrite(hrs, sep="\t", quote=FALSE, file="analyses/nmr_score_training/hazard_ratio_sanity_check.txt")
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/nmr_score_training/cindex_sanity_check.txt")

