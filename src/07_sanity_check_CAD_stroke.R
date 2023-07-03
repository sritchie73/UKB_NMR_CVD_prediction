library(data.table)
library(survival)
library(ggplot2)
source("src/utils/cox_test.R")
source("src/utils/score_cindex.R")

# Load in CAD and Stroke follow-up data and SCORE2
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "incident_cad_followup", "incident_cad", "incident_stroke_followup", "incident_stroke", "SCORE2"))

# Add in test scores
test_scores <- fread("analyses/nmr_score_training/aggregate_test_NMR_scores.txt")
dat <- dat[test_scores, on = .(eid)]

# Compute C-indices from scores directly as linear predictors, (e.g., SCORE2, or SCORE2 + CAD_NMR_score, or SCORE2 + Stroke_NMR_score)
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

# Compute hazard ratios by fitting Cox proportional hazards models to scores as dependent variables
# These are first standardised across the whole dataset, so that HRs and per SD, and that SD is comparable
# when comparing males to females to stratified models
dat[, SCORE2 := scale(SCORE2)]
dat[, CAD_NMR_score := scale(CAD_NMR_score)]
dat[, Stroke_NMR_score := scale(Stroke_NMR_score)]

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

# Repeat, using the scores compute from coefficient averages (what would be distributed on
# publication as the representative score, but unsuitable for testing in UK Biobank due to 
# overfitting to training samples)
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "incident_cad_followup", "incident_cad", "incident_stroke_followup", "incident_stroke", "SCORE2"))
test_scores <- fread("analyses/nmr_score_training/coef_avg_NMR_scores.txt")
dat <- dat[test_scores, on = .(eid)]

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

dat[, SCORE2 := scale(SCORE2)]
dat[, CAD_NMR_score := scale(CAD_NMR_score)]
dat[, Stroke_NMR_score := scale(Stroke_NMR_score)]

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

# Combine results from aggregate and average scores
hrs <- rbind(idcol="type", "aggregate"=hrs, "average"=hrs2)
cinds <- rbind(idcol="type", "aggregate"=cinds, "average"=cinds2)

hrs[, coefficient := factor(coefficient, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score"))]
hrs[, sex := factor(sex, levels=c("males", "females", "strata"))]

cinds[score == "NMR", score := sprintf("SCORE2 +\n%s_NMR_score", endpoint)]
cinds[, score := factor(score, levels=c("SCORE2", "SCORE2 +\nCAD_NMR_score", "SCORE2 +\nStroke_NMR_score"))]
cinds[, sex := factor(sex, levels=c("males", "females", "strata"))]

# Show and compare hazard ratios
g <- ggplot(hrs) +
  aes(x=coefficient, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_wrap(~ type + endpoint, nrow=1, scales="free_x") +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) + 
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("") + ylab("Hazard Ratio (95% CI) per SD") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), 
    axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g, width=7, height=3, file="analyses/nmr_score_training/hazard_ratio_sanity_check.pdf")

# Show without SCORE2
g <- ggplot(hrs[coefficient != "SCORE2"]) +
  aes(x=coefficient, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_wrap(~ type + endpoint, nrow=1, scales="free_x") +
  geom_hline(yintercept=1, linetype=2) +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) + 
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("") + ylab("Hazard Ratio (95% CI) per SD") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), 
    axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g, width=7, height=3, file="analyses/nmr_score_training/hazard_ratio_sanity_check_drop_score2.pdf")

# Show C-index
g <- ggplot(cinds) +
  aes(x=score, y=C.index, ymin=L95, ymax=U95, color=sex) + 
  facet_wrap(~ type + endpoint, nrow=1, scales="free_x") +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) +
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("") + ylab("C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), 
    axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g, width=7, height=3, file="analyses/nmr_score_training/cindex_sanity_check.pdf")

# Show change in C-index relative to SCORE2
score2_cinds <- cinds[score == "SCORE2"]
score2_cinds[, score := sprintf("SCORE2 +\n%s_NMR_score", endpoint)]
cinds[score2_cinds, on = .(type, score, endpoint, sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]

g <- ggplot(cinds[!is.na(deltaC)]) +
  aes(x=score, y=deltaC, ymin=deltaC.L95, ymax=deltaC.U95, color=sex) + 
  facet_wrap(~ type + endpoint, nrow=1, scales="free_x") +
  geom_hline(yintercept=0, linetype=2) +
  geom_errorbar(width=0, position=position_dodge(width=0.8)) +
  geom_point(shape=23, position=position_dodge(width=0.8)) +
  scale_color_manual("Sex", values=c("males"="#e41a1c", "females"="#377eb8", "strata"="#4daf4a")) +
  xlab("") + ylab("Delta C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=6), axis.title.y=element_text(size=8), 
    axis.text.x=element_text(size=6, angle=90, hjust=1, vjust=0.5),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    legend.text=element_text(size=6), legend.title=element_text(size=8),
    legend.position="bottom"
  )
ggsave(g, width=7, height=3, file="analyses/nmr_score_training/delta_cindex_sanity_check.pdf")

# Write out tables
fwrite(hrs, sep="\t", quote=FALSE, file="analyses/nmr_score_training/hazard_ratio_sanity_check.txt")
cinds[, score := gsub("\\\n", " ", score)]
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/nmr_score_training/cindex_sanity_check.txt")

