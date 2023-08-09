library(data.table)
library(foreach)
library(ggplot2)
library(ggstance)
source('src/utils/score_cindex.R')

# Make output directory
system("mkdir -p analyses/test")

# Load required data
dat <- fread("analyses/CVD_weight_training/CVD_linear_predictors_and_risk.txt")

# Compute C-indices
cinds <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %:%
  foreach(this_model=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %:% 
    foreach(this_score_type=c("non-derived", "clinical"), .combine=rbind) %do% {
      this_dat <- dat[model == this_model & score_type == this_score_type]
      this_info <- data.table(sex=this_sex,  model=this_model, score_type=this_score_type)
      if (this_sex == "Sex-stratified") {
        this_cind <- score_cindex("Surv(incident_cvd_followup, incident_cvd) ~ strata(sex) + linear_predictor", this_dat)
      } else {
        this_cind <- score_cindex("Surv(incident_cvd_followup, incident_cvd) ~ linear_predictor", this_dat[sex == gsub("s$", "", this_sex)])
      }
      cbind(this_info, this_cind)
}

# Compute delta C-index from SCORE2
ref <- cinds[model == "SCORE2"]
cinds[ref, on = .(sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]
cinds[model == "SCORE2", c("deltaC", "deltaC.L95", "deltaC.U95") := NA]

# Compute % improvement over SCORE2
cinds[ref, on = .(sex), c("pct_change", "pct.L95", "pct.U95") := .(deltaC/(i.C.index - 0.5)*100, deltaC.L95/(i.C.index - 0.5)*100, deltaC.U95/(i.C.index - 0.5)*100)]
cinds[model == "SCORE2", c("pct_change", "pct.L95", "pct.U95") := NA]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/cindices.txt")

# Encode factors for plotting
cinds[, model := factor(model, levels=rev(c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")))]
cinds[, sex := factor(sex, levels=c("Males", "Females", "Sex-stratified"))]
ref[, model := NULL]
ref[, sex := factor(sex, levels=c("Males", "Females", "Sex-stratified"))]

# Plot full biomarker scores
g <- ggplot(cinds[score_type == "non-derived"]) +
  aes(x=C.index, xmin=L95, xmax=U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1, scales="free_x") +
  geom_vline(data=ref, aes(xintercept=C.index), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") + 
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8", "Sex-stratified"="#006d2c")) +
  ylab("") + xlab("C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=2, file="analyses/test/cindices_non_derived.pdf")

# Plot change in C-index
g <- ggplot(cinds[model != "SCORE2" & score_type == "non-derived"]) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8", "Sex-stratified"="#006d2c")) +
  ylab("") + xlab("ΔC-index over SCORE2 (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=1.8, file="analyses/test/delta_cindices_non_derived.pdf", device=cairo_pdf)

# Plot % improvement
g <- ggplot(cinds[model != "SCORE2" & score_type == "non-derived"]) +
  aes(x=pct_change, xmin=pct.L95, xmax=pct.U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8", "Sex-stratified"="#006d2c")) +
  ylab("") + xlab("% improvement in C-index over SCORE2 (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=1.8, file="analyses/test/pct_improvement_cindices_non_derived.pdf")

# Create plots showing clinical biomarker scores have worse performance
cinds[, score_type := fcase(
  !(model %like% "NMR"), "Model does not include NMR biomarker scores",
  score_type == "clinical", "NMR biomarker scores trained from 21 clinically accredited biomarkers",
  score_type == "non-derived", "NMR biomarker scores trained from 106 fractionated biomarkers"
)]
cinds[, score_type := factor(score_type, levels=c(
  "Model does not include NMR biomarker scores", "NMR biomarker scores trained from 106 fractionated biomarkers",
  "NMR biomarker scores trained from 21 clinically accredited biomarkers"
))]
cinds <- unique(cinds)

g <- ggplot(cinds) +
  aes(x=C.index, xmin=L95, xmax=U95, y=model, color=score_type) +
  facet_wrap(~ sex, nrow=1, scales="free_x") +
  geom_vline(data=ref, aes(xintercept=C.index), linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.6)) +
  geom_point(shape=23, size=2, fill="white", position=position_dodgev(height=0.6)) + 
  scale_color_manual(values=c(
    "Model does not include NMR biomarker scores"="black",
    "NMR biomarker scores trained from 106 fractionated biomarkers"="#4575b4",
    "NMR biomarker scores trained from 21 clinically accredited biomarkers"="#f46d43"
  )) +
  ylab("") + xlab("C-index (95% CI)") +
  guides(color=guide_legend(title="", nrow=3)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.text=element_text(size=6, color="black") 
  )
ggsave(g, width=7.2, height=4, file="analyses/test/cindices_with_clinical_scores.pdf")

g <- ggplot(cinds[model != "SCORE2"]) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=model, color=score_type) +
  facet_wrap(~ sex, nrow=1) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.6)) +
  geom_point(shape=23, size=2, fill="white", position=position_dodgev(height=0.6)) + 
  scale_color_manual(values=c(
    "Model does not include NMR biomarker scores"="black",
    "NMR biomarker scores trained from 106 fractionated biomarkers"="#4575b4",
    "NMR biomarker scores trained from 21 clinically accredited biomarkers"="#f46d43"
  )) +
  ylab("") + xlab("ΔC-index over SCORE2 (95% CI)") +
  guides(color=guide_legend(title="", nrow=3)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.text=element_text(size=6, color="black") 
  )
ggsave(g, width=7.2, height=3, file="analyses/test/delta_cindices_with_clinical_scores.pdf", device=cairo_pdf)

g <- ggplot(cinds[model != "SCORE2"]) +
  aes(x=pct_change, xmin=pct.L95, xmax=pct.U95, y=model, color=score_type) +
  facet_wrap(~ sex, nrow=1) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.6)) +
  geom_point(shape=23, size=2, fill="white", position=position_dodgev(height=0.6)) + 
  scale_color_manual(values=c(
    "Model does not include NMR biomarker scores"="black",
    "NMR biomarker scores trained from 106 fractionated biomarkers"="#4575b4",
    "NMR biomarker scores trained from 21 clinically accredited biomarkers"="#f46d43"
  )) +
  ylab("") + xlab("% improvement in C-index over SCORE2 (95% CI)") +
  guides(color=guide_legend(title="", nrow=3)) +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="bottom", legend.text=element_text(size=6, color="black") 
  )
ggsave(g, width=7.2, height=3, file="analyses/test/pct_improvement_cindices_with_clinical_scores.pdf", device=cairo_pdf)

# Create formatted table for manuscript
dt <- cinds[score_type != "NMR biomarker scores trained from 21 clinically accredited biomarkers"]
dt[, score_type := NULL]
dt[, SE := NULL]
dt[, pct_change := pct_change / 100]
dt[, pct.L95 := pct.L95 / 100]
dt[, pct.U95 := pct.U95 / 100]
fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/cindices_for_supp.txt")

