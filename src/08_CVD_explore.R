library(data.table)
library(foreach)
library(survival)
library(ggplot2)
library(forcats)
library(ggrastr)
library(ggpp)
library(ggh4x)
source("src/utils/cox_test.R")

options("ggrastr.default.dpi" = 300)

# Make output directory
system("mkdir -p analyses/CVD_score_weighting/explore", wait=TRUE)

# Load in scores and CVD data
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "incident_cvd_followup", "incident_cvd", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2"))
nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]

# Compute hazard ratios by fitting Cox proportional hazards models to scores as dependent variables
# These are first standardised across the whole dataset, so that HRs and per SD, and that SD is comparable
# when comparing males to females to stratified models
dat[, SCORE2 := scale(SCORE2)]
dat[, CAD_NMR_score := scale(CAD_NMR_score)]
dat[, Stroke_NMR_score := scale(Stroke_NMR_score)]
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]

# Sanity check HRs
cx_list <- list(
  "NMR_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Male"]),
  "NMR_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Female"]),
  "NMR_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + SCORE2 + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat),
  "PRS_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"]),
  "PRS_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"]),
  "PRS_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + SCORE2 + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat),
  "both_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"]),
  "both_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ SCORE2 + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"]),
  "both_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + SCORE2 + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat)
)

hrs <- rbindlist(idcol="model", lapply(cx_list, `[[`, "coefficients"))
hrs[, c("score", "endpoint", "sex") := tstrsplit(model, "_")]

hrs[, coefficient := factor(coefficient, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]
hrs[, score := factor(score, levels=c("NMR", "PRS", "both"))]
hrs[, sex := factor(sex, levels=c("males", "females", "strata"))]

g <- ggplot(hrs) +
  aes(x=coefficient, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_wrap(~ score, nrow=1, scales="free_x") +
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
ggsave(g, width=7, height=3, file="analyses/CVD_score_weighting/explore/hazard_ratio_sanity_check.pdf")

g <- ggplot(hrs[coefficient != "SCORE2"]) +
  aes(x=coefficient, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_wrap(~ score, nrow=1, scales="free_x") +
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
ggsave(g, width=7, height=3, file="analyses/CVD_score_weighting/explore/hazard_ratio_sanity_check_drop_score2.pdf")

# Try SCORE2 as offset
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "incident_cvd_followup", "incident_cvd", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2"))
nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]

dat[, CAD_NMR_score := scale(CAD_NMR_score)]
dat[, Stroke_NMR_score := scale(Stroke_NMR_score)]
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]

# Sanity check HRs
cx_list2 <- list(
  "NMR_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Male"]),
  "NMR_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Female"]),
  "NMR_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat),
  "PRS_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"]),
  "PRS_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"]),
  "PRS_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat),
  "both_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"]),
  "both_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"]),
  "both_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat)
)

hrs2 <- rbindlist(idcol="model", lapply(cx_list2, `[[`, "coefficients"))
hrs2[, c("score", "endpoint", "sex") := tstrsplit(model, "_")]

hrs2[, coefficient := factor(coefficient, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]
hrs2[, score := factor(score, levels=c("NMR", "PRS", "both"))]
hrs2[, sex := factor(sex, levels=c("males", "females", "strata"))]

hrs <- rbind(idcol="type", "SCORE2 as dependent variable"=hrs, "SCORE2 as offset"=hrs2)

g <- ggplot(hrs[coefficient != "SCORE2"]) +
  aes(x=coefficient, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_wrap(~ type + score, nrow=1, scales="free_x") +
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
ggsave(g, width=7, height=3, file="analyses/CVD_score_weighting/explore/hazard_ratio_offset_vs_dependent.pdf")

# Compare scores derived from coefficient averages, which are most likely overfit
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "incident_cvd_followup", "incident_cvd", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2"))
nmr_scores <- fread("analyses/nmr_score_training/coef_avg_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]

dat[, CAD_NMR_score := scale(CAD_NMR_score)]
dat[, Stroke_NMR_score := scale(Stroke_NMR_score)]
dat[, CAD_metaGRS := scale(CAD_metaGRS)]
dat[, Stroke_metaGRS := scale(Stroke_metaGRS)]

# Sanity check HRs
cx_list3 <- list(
  "NMR_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Male"]),
  "NMR_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat[sex == "Female"]),
  "NMR_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score", "incident_cvd", dat),
  "PRS_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"]),
  "PRS_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"]),
  "PRS_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2) + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat),
  "both_CVD_males"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Male"]),
  "both_CVD_females"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat[sex == "Female"]),
  "both_CVD_strata"=cox.test("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + offset(SCORE2) + CAD_NMR_score + Stroke_NMR_score + CAD_metaGRS + Stroke_metaGRS", "incident_cvd", dat)
)

hrs3 <- rbindlist(idcol="model", lapply(cx_list3, `[[`, "coefficients"))
hrs3[, c("score", "endpoint", "sex") := tstrsplit(model, "_")]

hrs3[, coefficient := factor(coefficient, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]
hrs3[, score := factor(score, levels=c("NMR", "PRS", "both"))]
hrs3[, sex := factor(sex, levels=c("males", "females", "strata"))]
hrs3[, type := "SCORE2 as dependent variable"]

hrs <- rbind(idcol="score_type", "aggregate"=hrs, "average"=hrs3)

g <- ggplot(hrs[coefficient != "SCORE2" & type == "SCORE2 as dependent variable" & score != "PRS"]) +
  aes(x=coefficient, y=HR, ymin=L95, ymax=U95, color=sex) +
  facet_wrap(~ score_type + score, nrow=1, scales="free_x") +
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
ggsave(g, width=7, height=3, file="analyses/CVD_score_weighting/explore/hazard_ratio_aggregate_vs_average.pdf")

# Write out table of all hazard ratios
fwrite(hrs, sep="\t", quote=FALSE, file="analyses/CVD_score_weighting/explore/hazard_ratios.txt")

# Compare scores to each other
dat <- fread("data/cleaned/analysis_cohort.txt", select=c("eid", "sex", "CAD_metaGRS", "Stroke_metaGRS", "SCORE2"))
nmr_scores <- fread("analyses/nmr_score_training/aggregate_test_NMR_scores.txt")
dat <- dat[nmr_scores, on = .(eid)]
dat <- melt(dat, id.vars=c("eid", "sex"), variable.name="score")

score_comp <- foreach(this_rn = unique(dat$score), .combine=rbind) %:%
  foreach(this_cn = unique(dat$score), .combine=rbind) %do% {
    this_x <- dat[score == this_rn, .(eid, sex, x_score=score, x_value=value)]
    this_y <- dat[score == this_cn, .(eid, sex, y_score=score, y_value=value)]
    merge(this_x, this_y, by=c("eid", "sex"))
}
score_comp[, x_score := factor(x_score, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]
score_comp[, y_score := factor(y_score, levels=c("SCORE2", "CAD_NMR_score", "Stroke_NMR_score", "CAD_metaGRS", "Stroke_metaGRS"))]

score_scatter <- function(score_comp) {
  cor_anno <- score_comp[,.(label=sprintf("r=%.2f\np=%.2f",
    cor(x_value, y_value),
    cor(x_value, y_value, method="spearman")
  )), by=.(x_score, y_score)]

  ggplot(score_comp) +
    aes(x=x_value, y=y_value) +
    rasterise(geom_point(alpha=0.5, shape=19, size=0.1)) +
    geom_smooth(method="lm", color="red", linetype=2) + 
    geom_text_npc(data=cor_anno, npcx="left", npcy="top", aes(label=label), color="red", size=6*0.352777778) +
    facet_grid2(fct_rev(y_score) ~ x_score, scales="free", independent="all") +
    xlab("") + ylab("") +
    theme_bw() +
    theme(
      axis.text=element_text(size=6), axis.title=element_text(size=8),
      strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
      panel.grid=element_blank(), legend.text=element_text(size=6)
    )
}

g1 <- score_scatter(score_comp[sex == "Male"])
g2 <- score_scatter(score_comp[sex == "Female"])

ggsave(g1, width=7, height=7, file="analyses/CVD_score_weighting/explore/score_compare_males.pdf")
ggsave(g2, width=7, height=7, file="analyses/CVD_score_weighting/explore/score_compare_females.pdf")







