library(data.table)
library(ggplot2)
library(ggh4x)
library(cowplot)

# Load results
blanket <- fread("analyses/public_health_modelling/blanket_screening/population_screening.txt")
targeted <- fread("analyses/public_health_modelling/targeted_screening/population_screening.txt")

# Combine
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)

# Format columns of interest
shownum <- function(num) { format(num, big.mark=",", trim=TRUE) }
dat <- dat[,.(
  strategy, model, sex, 
  high_risk = sprintf("%s (%s-%s)", shownum(high_risk), shownum(high_risk_L95), shownum(high_risk_U95)),
  cases = sprintf("%s (%s-%s)", shownum(high_risk_cases), shownum(high_risk_cases_L95), shownum(high_risk_cases_U95)),
  prevented = sprintf("%s (%s-%s)", shownum(events_prevented), shownum(events_prevented_L95), shownum(events_prevented_U95)),
  NNS = sprintf("%s (%s-%s)", shownum(NNS), shownum(NNS_L95), shownum(NNS_U95)),
  NNT = sprintf("%s (%s-%s)", shownum(NNT), shownum(NNT_L95), shownum(NNT_U95))
)]

# Impose ordering
dat[, model := factor(model, levels=c("SCORE2", "SCORE2 + NMR scores", "SCORE2 + PRSs", "SCORE2 + NMR scores + PRSs"))]
dat[, sex := factor(sex, levels=c("Males", "Females"))]
dat[, strategy := factor(strategy, levels=c("blanket", "targeted"))]
dat <- dat[order(sex)][order(model)][order(strategy)]

dat <- dcast(dat, strategy + model ~ sex, value.var=c("high_risk", "cases", "prevented", "NNS", "NNT"))
dat <- dat[, .(strategy, model, high_risk_Males, cases_Males,  prevented_Males, NNS_Males, NNT_Males,
  high_risk_Females, cases_Females, prevented_Females, NNS_Females, NNT_Females)]

# Write out
fwrite(dat, sep="\t", quote=FALSE, file="analyses/public_health_modelling/screening_summary.txt")

# Also plot
dat <- rbind(idcol="strategy", "blanket"=blanket, "targeted"=targeted)
dat <- rbind(idcol="metric",
  "People classified as high risk"=dat[,.(strategy, sex, model, estimate=high_risk, L95=high_risk_L95, U95=high_risk_U95)],
  "CVD cases classified as high risk"=dat[,.(strategy, sex, model, estimate=high_risk_cases, L95=high_risk_cases_L95, U95=high_risk_cases_U95)],
  "CVD events prevented by statin prescription"=dat[,.(strategy, sex, model, estimate=events_prevented, L95=events_prevented_L95, U95=events_prevented_U95)],
  "Number needed to screen to prevent 1 CVD event"=dat[,.(strategy, sex, model, estimate=NNS, L95=NNS_L95, U95=NNS_U95)],
  "Number of statins prescribed per CVD event prevented"=dat[,.(strategy, sex, model, estimate=NNT, L95=NNT_L95, U95=NNT_U95)]
)
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
ggsave(g, width=7.2, height=4, file="analyses/public_health_modelling/screening_comparison.pdf")






