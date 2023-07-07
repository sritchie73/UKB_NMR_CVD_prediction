library(data.table)
library(foreach)
library(ggplot2)
source('src/utils/score_cindex.R')

# Make output directory
system("mkdir -p analyses/test")

# Load required data
dat <- fread("analyses/CVD_score_weighting/CVD_linear_predictors_and_risk.txt")

# Compute C-indices
cinds <- foreach(this_sex=c("Males", "Females"), .combine=rbind) %:%
  foreach(this_model=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"), .combine=rbind) %do% {
    cbind(sex=this_sex, model=this_model, score_cindex("Surv(incident_cvd_followup, incident_cvd) ~ linear_predictor", dat[sex == gsub("s$", "", this_sex) & model == this_model]))
}

# Compute delta C-index from SCORE2
ref <- cinds[model == "SCORE2"]
cinds[ref, on = .(sex), c("deltaC", "deltaC.L95", "deltaC.U95") := .(C.index - i.C.index, L95 - i.C.index, U95 - i.C.index)]
cinds[model == "SCORE2", c("deltaC", "deltaC.L95", "deltaC.U95") := NA]

# Compute % improvement over SCORE2
cinds[ref, on = .(sex), c("pct_change", "pct.L95", "pct.U95") := .(C.index/i.C.index*100-100, L95/i.C.index*100-100, U95/i.C.index*100-100)]
cinds[model == "SCORE2", c("pct_change", "pct.L95", "pct.U95") := NA]

# Write out
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/cindices.txt")

# Plot
cinds[, model := factor(model, levels=rev(c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")))]
cinds[, sex := factor(sex, levels=c("Males", "Females"))]
ref[, model := NULL]
ref[, sex := factor(sex, levels=c("Males", "Females"))]
g <- ggplot(cinds) +
  aes(x=C.index, xmin=L95, xmax=U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1, scales="free_x") +
  geom_vline(data=ref, aes(xintercept=C.index), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") + 
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8")) +
  ylab("") + xlab("C-index (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=2, file="analyses/test/cindices.pdf")

# Plot change in C-index
g <- ggplot(cinds[model != "SCORE2"]) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8")) +
  ylab("") + xlab("ΔC-index over SCORE2 (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=1.8, file="analyses/test/delta_cindices.pdf", device=cairo_pdf)

# Plot % improvement
g <- ggplot(cinds[model != "SCORE2"]) +
  aes(x=pct_change, xmin=pct.L95, xmax=pct.U95, y=model, color=sex) +
  facet_wrap(~ sex, nrow=1) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_color_manual("Sex", values=c("Males"="#e41a1c", "Females"="#377eb8")) +
  ylab("") + xlab("% improvement in C-index over SCORE2 (95% CI)") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.text=element_text(size=8, face="bold"), strip.background=element_blank(),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=1.8, file="analyses/test/pct_improvement_cindices.pdf")

