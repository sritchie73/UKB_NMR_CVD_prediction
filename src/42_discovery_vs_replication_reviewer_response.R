library(data.table)
library(foreach)
library(forcats)
library(ggplot2)
library(scales)
library(ggstance)
library(patchwork)

# Build C-index figure
cinds <- rbind(idcol="cohort",
  "discovery"=fread("analyses/test/discovery_cindices.txt"),
  "replication"=fread("analyses/test/replication_cindices.txt")
)

ggdt <- cinds[model_sex == "Sex-stratified" & endpoint == "cvd" & score == "SCORE2_excl_UKB"]
ggdt[, model_name := fcase(
  model == "SCORE2 + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "SCORE2 + Biochemistry", "SCORE2 + 11 clinical chemistry biomarkers",
  model == "SCORE2 + PRS", "SCORE2 + PRSs for CHD and IS",
  model == "SCORE2 + NMR scores + PRS", "SCORE2 + NMR scores for CHD and IS + PRSs for CHD and IS",
  model == "SCORE2 + Biochemistry + PRS", "SCORE2 + 11 clinical chemistry biomarkers + PRSs for CHD and IS"
)]
ggdt[,model_name := factor(model_name, levels=model_name)]
ggdt[,cohort := factor(cohort, levels=c("discovery", "replication"))]

g <- ggplot(ggdt) +
  aes(x=deltaC, xmin=deltaC.L95, xmax=deltaC.U95, y=fct_rev(model_name), color=model) +
  facet_grid(cohort ~ .) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, size=2, fill="white") +
  scale_x_continuous("Change in C-index (95% CI)", limits=c(0, 0.03), breaks=c(0, 0.01, 0.02, 0.03)) +
  scale_color_manual(values=c(
    "SCORE2 + NMR scores"="#e41a1c", "SCORE2 + PRS"="#377eb8", "SCORE2 + Biochemistry"="#ff7f00",
    "SCORE2 + NMR scores + PRS"="#4daf4a", "SCORE2 + Biochemistry + PRS"="#984ea3"
  )) +
  ylab("") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none", strip.background=element_blank(), strip.text=element_blank()
  )
ggsave(g, width=5.3, height=4, file="analyses/test/discovery_vs_replication_delta_cindices.pdf")


#  Build NRI figure
nri <- rbind(idcol="cohort",
  "discovery"=fread("analyses/test/discovery_categorical_nri_estimates.txt"),
  "replication"=fread("analyses/test/replication_categorical_nri_estimates.txt")
)

ggdt <- nri[model_sex == "Sex-stratified" & endpoint == "cvd" & score == "SCORE2_excl_UKB"]
ggdt <- ggdt[metric %in% c("NRI+", "NRI-")]
ggdt[, model_name := fcase(
  model == "SCORE2 + NMR scores", "SCORE2 + NMR scores for CHD and IS",
  model == "SCORE2 + Biochemistry", "SCORE2 + 11 clinical chemistry biomarkers",
  model == "SCORE2 + PRS", "SCORE2 + PRSs for CHD and IS",
  model == "SCORE2 + NMR scores + PRS", "SCORE2 + NMR scores for CHD and IS + PRSs for CHD and IS",
  model == "SCORE2 + Biochemistry + PRS", "SCORE2 + 11 clinical chemistry biomarkers + PRSs for CHD and IS"
)]
ggdt[,model_name := factor(model_name, levels=unique(model_name))]
ggdt[,cohort := factor(cohort, levels=c("discovery", "replication"))]

g <- ggplot(ggdt) +
  aes(x=Estimate, xmin=L95, xmax=U95, y=fct_rev(model_name), color=metric) +
  facet_grid(cohort ~ .) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0, position=position_dodgev(height=0.3)) +
  geom_point(shape=23, size=2, fill="white", position=position_dodgev(height=0.3)) +
  scale_color_manual(values=c("NRI+"="#c51b7d", "NRI-"="#4d9221")) +
  scale_x_continuous("Categorical NRI, % reclassified (95% CI)", labels=percent, limits=c(-0.05, 0.17)) +
  ylab("") +
  theme_bw() +
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    legend.position="none", strip.background=element_blank(), strip.text=element_blank()
  )
ggsave(g, width=5.3, height=4, file="analyses/test/discovery_vs_replication_categorical_nri.pdf")

