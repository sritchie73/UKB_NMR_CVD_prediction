library(data.table)
library(ggplot2)
library(ggh4x)
library(cowplot)

# Load results stratified by sex
blanket <- fread("analyses/public_health_modelling/blanket_screening/population_screening_by_sex.txt")
targeted <- fread("analyses/public_health_modelling/targeted_screening/population_screening_by_sex.txt")

# Combine
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)

# Format columns of interest
shownum <- function(num) { format(num, big.mark=",", trim=TRUE) }
dat[, text := sprintf("%s (%s-%s)", shownum(estimate), shownum(L95), shownum(U95))]
dat[, people := shownum(people)]
dat[, cases := shownum(cases)]
dat <- dcast(dat, sex + people + cases + strategy + model ~ number, value.var="text")
dat <- dat[, .(sex, people, cases, strategy, model, high_risk, high_risk_cases, events_prevented, NNS, NNT, delta_high_risk_cases, delta_events_prevented)]

# Impose ordering
dat[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dat[, sex := factor(sex, levels=c("Males", "Females"))]
dat[, strategy := factor(strategy, levels=c("blanket", "targeted"))]
dat <- dat[order(model)][order(strategy)][order(sex)]

# Write out
fwrite(dat, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_summary_by_sex.txt")

# Also plot
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)
dat[, metric := fcase(
  number == "high_risk", "People classified as high risk",
  number == "high_risk_cases", "CVD cases classified as high risk",
  number == "events_prevented", "CVD events prevented by statin prescription",
  number == "NNS", "Number needed to screen to prevent 1 CVD event",
  number == "NNT", "Number of statins prescribed per CVD event prevented",
  default = NA
)]
dat <- dat[!is.na(metric)]
dat[, metric := factor(metric, levels=unique(metric))]
dat[, strategy := fcase(
  strategy == "blanket", "Population screening",
  strategy == "targeted", "Targeted screening"
)]
dat[, model := factor(model, levels=rev(c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs")))]

ref <- dat[model == "SCORE2"]
ref[, strategy := NULL]

g1 <- ggplot(dat[sex == "Males"]) + 
  aes(x=estimate, xmin=L95, xmax=U95, y=model, color=metric) +
  facet_grid(strategy ~ metric, scales="free", labeller=labeller(strategy=label_wrap_gen(width=10), metric=label_wrap_gen(width=20))) +
  geom_vline(data=ref[sex == "Males"], aes(xintercept=estimate), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white", size=1.2) +
  scale_x_continuous("Number (95% CI)", n.breaks=4) +
  theme_bw() + 
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=6, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), 
    legend.position="none"
  )

g2 <- ggplot(dat[sex == "Females"]) + 
  aes(x=estimate, xmin=L95, xmax=U95, y=model, color=metric) +
  facet_grid(strategy ~ metric, scales="free", labeller=labeller(strategy=label_wrap_gen(width=10), metric=label_wrap_gen(width=20))) +
  geom_vline(data=ref[sex == "Females"], aes(xintercept=estimate), linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(shape=23, fill="white", size=1.2) +
  scale_x_continuous("Number (95% CI)", n.breaks=4) +
  theme_bw() + 
  theme(
    axis.text.y=element_text(size=8, color="black"), axis.title.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    strip.background=element_blank(), strip.text=element_text(size=6, face="bold"),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), 
    legend.position="none"
  )

g <- plot_grid(g1, g2, labels=c("Males  ", "Females"), label_size=10, ncol=1)
ggsave(g, width=7.2, height=4, file="analyses/public_health_modelling/screening_comparison_by_sex.pdf")

# Load results not stratified by sex
blanket <- fread("analyses/public_health_modelling/blanket_screening/population_screening.txt")
targeted <- fread("analyses/public_health_modelling/targeted_screening/population_screening.txt")

# Combine
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)

# Format columns of interest
dat[, text := sprintf("%s (%s-%s)", shownum(estimate), shownum(L95), shownum(U95))]
dat[, people := shownum(people)]
dat[, cases := shownum(cases)]
dat <- dcast(dat, people + cases + strategy + model ~ number, value.var="text")
dat <- dat[, .(people, cases, strategy, model, high_risk, high_risk_cases, events_prevented, NNS, NNT, delta_high_risk_cases, delta_events_prevented)]

# Impose ordering
dat[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dat[, strategy := factor(strategy, levels=c("blanket", "targeted"))]
dat <- dat[order(model)][order(strategy)]

# Write out
fwrite(dat, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_summary.txt")

