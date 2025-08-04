library(data.table)
library(foreach)
library(forcats)
library(ggplot2)
library(ggstance)
library(patchwork)

# Load and collate results
res <- rbind(idcol="cohort",
  "discovery"=fread("analyses/test/discovery_cindices_with_assay_scores.txt"),
  "replication"=fread("analyses/test/replication_cindices_with_assay_scores.txt"),
  "pooled"=fread("analyses/test/pooled_cindices_with_assay_scores.txt")
)

# Build main figure 3A
ggdt <- res[cohort == "pooled" & model_sex == "Sex-stratified" & endpoint == "cvd" & score == "SCORE2_excl_UKB"]
ggdt[, model_name := fcase(
  model == "SCORE2 + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "SCORE2 + Assay scores", "SCORE2 + Clinical assay scores for CHD and IS",
  model == "SCORE2 + PRS", "SCORE2 + PRSs for CHD and IS",
  model == "SCORE2 + NMR scores + Assay scores", "SCORE2 + NMR scores for CHD and IS + Clinical assay scores for CHD and IS",
  model == "SCORE2 + NMR scores + PRS", "SCORE2 + NMR scores for CHD and IS + PRSs for CHD and IS",
  model == "SCORE2 + Assay scores + PRS", "SCORE2 + Clinical assay scores for CHD and IS + PRSs for CHD and IS",
  model == "SCORE2 + NMR scores + Assay scores + PRS", "SCORE2 + NMR scores for CHD and IS + Clinical assay scores for CHD and IS + PRSs for CHD and IS"
)]
ggdt[,model_name := factor(model_name, levels=model_name)]

g <- ggplot(ggdt) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=fct_rev(model_name), color=model) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_x_continuous("Change in C-index (95% CI)", limits=c(0, 0.03), breaks=c(0, 0.01, 0.02, 0.03)) +
	scale_color_manual(values=c(
    "SCORE2"="black", "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRS"="#377eb8", "SCORE2 + Assay scores"="#ff7f00",
    "SCORE2 + NMR scores + PRS"="#4daf4a", "SCORE2 + Assay scores + PRS"="#984ea3", "SCORE2 + NMR scores + Assay scores"="#f781bf",
    "SCORE2 + NMR scores + Assay scores + PRS"="#a65628"
	)) +
  ylab("") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none"
  )
ggsave(g, width=7.2, height=2.2, file="analyses/test/main_delta_cindices_with_assay_scores.pdf")

