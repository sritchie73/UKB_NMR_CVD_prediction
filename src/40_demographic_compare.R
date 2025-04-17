library(data.table)
library(ggplot2)
library(patchwork)

# Compute numbers and CVD cases per 100,000 in each age-group/sex in the different UKB cohorts used
ukb <- rbind(idcol="cohort", fill=TRUE,
  "discovery"=fread("data/cleaned/analysis_cohort.txt"),
  "replication"=fread("data/cleaned/phase3_analysis_cohort.txt")
)
ukb_pooled <- copy(ukb)
ukb_pooled[, cohort := "pooled"]
ukb <- rbind(ukb, ukb_pooled)

ukb[, age_group := sprintf("%s-%s", age %/% 5 * 5, age %/% 5 * 5 + 4)]
ukb <- ukb[, .(N=.N, cases=sum(incident_cvd)), by=.(sex, age_group, cohort)]

ukb[, N := as.numeric(N)]
ukb[, cases := as.numeric(cases)]

totals <- ukb[,.(total=sum(N)),by=cohort]
ukb[totals, on = .(cohort), N := N/total * 100000]
ukb[totals, on = .(cohort), cases := cases/total * 100000]
ukb <- ukb[order(age_group)][order(sex)]

# Load numbers per 100,000 from ONS/CPRD
ons <- fread("analyses/public_health_modelling/ONS_hypothetical_100k_pop_by_age_sex.txt")
ons[, controls := NULL]
ons <- ons[order(age_group)][order(sex)]

# Combine into one table and write out
comp <- rbind("UK population"=ons, "UK Biobank"=ukb, idcol="population", fill=TRUE)
fwrite(comp, sep="\t", quote=FALSE, file="analyses/public_health_modelling/demographic_compare.txt")

# Create plot comparing UKB (pooled cohort, discovery and replication are similar) to ONS/CPRD
ggdt <- comp[population == "UK population" | cohort == "pooled"]
ggdt[sex == "Male", cases := -cases]
ggdt[sex == "Male", N := -N]

g1 <- ggplot(ggdt[population == "UK Biobank"]) +
  aes(y=age_group, x=N, fill=sex) +
  geom_col() +
  geom_vline(xintercept=0, linetype=2) + 
  scale_x_continuous(name="N per 100,000", limits=c(-13500, 13500), breaks=c(-10000, -5000, 0, 5000, 10000), labels=abs) +
  scale_fill_manual(values=c("Male"="#9970ab", "Female"="#5aae61")) +
  ylab("Age-group") +
  ggtitle("UK Biobank") +
  theme_bw() + 
  theme(
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), plot.title=element_text(size=8, face="bold"),
    axis.text=element_text(size=6), axis.title=element_text(size=8), legend.position="none"
  )

g2 <- ggplot(ggdt[population == "UK population"]) +
  aes(y=age_group, x=N, fill=sex) +
  geom_col() +
  geom_vline(xintercept=0, linetype=2) + 
  scale_x_continuous(name="N per 100,000", limits=c(-13500, 13500), breaks=c(-10000, -5000, 0, 5000, 10000), labels=abs) +
  scale_fill_manual(values=c("Male"="#9970ab", "Female"="#5aae61")) +
  ylab("Age-group") +
  ggtitle("UK primary care population") +
  theme_bw() + 
  theme(
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), plot.title=element_text(size=8, face="bold"),
    axis.text=element_text(size=6), axis.title=element_text(size=8), legend.position="none"
  )

g3 <- ggplot(ggdt[population == "UK Biobank"]) +
  aes(y=age_group, x=cases, fill=sex) +
  geom_col() +
  geom_vline(xintercept=0, linetype=2) + 
  scale_x_continuous(name="CVD per 100,000", limits=c(-1250, 1250), breaks=c(-1000, -500, 0, 500, 1000), labels=abs) +
  scale_fill_manual(values=c("Male"="#9970ab", "Female"="#5aae61")) +
  ylab("Age-group") +
  ggtitle("UK Biobank") +
  theme_bw() + 
  theme(
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), plot.title=element_text(size=8, face="bold"),
    axis.text=element_text(size=6), axis.title=element_text(size=8), legend.position="none"
  )

g4 <- ggplot(ggdt[population == "UK population"]) +
  aes(y=age_group, x=cases, fill=sex) +
  geom_col() +
  geom_vline(xintercept=0, linetype=2) + 
  scale_x_continuous(name="CVD per 100,000", limits=c(-1250, 1250), breaks=c(-1000, -500, 0, 500, 1000), labels=abs) +
  scale_fill_manual(values=c("Male"="#9970ab", "Female"="#5aae61")) +
  ylab("Age-group") +
  ggtitle("UK primary care population") +
  theme_bw() + 
  theme(
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), plot.title=element_text(size=8, face="bold"),
    axis.text=element_text(size=6), axis.title=element_text(size=8), legend.position="none"
  )

g <- (g1 + g3) / (g2 + g4)

ggsave(g, width=7.2, height=5, file="analyses/public_health_modelling/demographic_compare.pdf")

