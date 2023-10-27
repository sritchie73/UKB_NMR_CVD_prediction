library(data.table)
library(foreach)
library(doMC)
library(survival)
library(ggplot2)

system("mkdir -p analyses/test/")

# Load data
dat <- fread("data/cleaned/analysis_cohort.txt")

# Build data.table of models to test
model_info <- foreach(this_sex=c("Males", "Females", "Sex-stratified"), .combine=rbind) %:%
  foreach(this_model=c("SCORE2", "SCORE2_excl_UKB"), .combine=rbind) %do% {
    data.table(sex=this_sex,  model=this_model)
}

# Convert sex to integer for strata
strata_num <- dat[,.GRP, by=sex]
dat[strata_num, on = .(sex), sex_int := i.GRP]

# Filter to relevant columns for testing
dat <- dat[,.(eid, sex, sex_int, incident_cvd, incident_cvd_followup, SCORE2, SCORE2_excl_UKB)]

# Comput C-index and jacknife SE for each model
cinds <- foreach(midx = model_info[,.I], .combine=rbind) %do% {
  this_minfo <- model_info[midx]
  if (this_minfo$sex == "Sex-stratified") {
    this_y <- Surv(dat[["incident_cvd_followup"]], dat[["incident_cvd"]])
    this_x <- dat[[this_minfo$model]]
    cf <- survival::concordancefit(y = this_y, x = this_x, strata = dat[["sex_int"]], reverse = TRUE, timefix = TRUE)
  } else {
    this_dt <- dat[sex == gsub("s$", "", this_minfo$sex)]
    this_y <- Surv(this_dt[["incident_cvd_followup"]], this_dt[["incident_cvd"]])
    this_x <- this_dt[[this_minfo$model]]
    cf <- survival::concordancefit(y = this_y, x = this_x, reverse = TRUE, timefix = TRUE)
  }
  cbind(this_minfo, C.index=cf$concordance, SE=sqrt(cf$var)) # jackknife SE computed same way as in coxph()
}

# Compute 95% CI
cinds[, L95 := C.index - qnorm(1-(0.05/2))*SE]
cinds[, U95 := C.index + qnorm(1-(0.05/2))*SE]

# organise columns and write out
cinds <- cinds[, .(sex, SCORE2_method = ifelse(model == "SCORE2", "Weights derived from all datasets", "Weights derived excluding UK Biobank"), C.index, SE, L95, U95)]
fwrite(cinds, sep="\t", quote=FALSE, file="analyses/test/cindex_by_SCORE2_method.txt")

# Plot
cinds[, sex := factor(sex, levels=c("Sex-stratified", "Males", "Females"))]
cinds[, SCORE2_method := factor(SCORE2_method, levels=c("Weights derived excluding UK Biobank", "Weights derived from all datasets"))]
g <- ggplot(cinds) + 
  aes(x=C.index, xmin=L95, xmax=U95, y=SCORE2_method, color=SCORE2_method) +
  facet_wrap(~ sex, scales="free_x") +
  geom_errorbarh(height=0) +
  geom_point(shape=18) +
  scale_color_manual(values=c("Weights derived excluding UK Biobank"="#2166ac", "Weights derived from all datasets"="#b2182b")) +
  xlab("C-index (95% CI)") +
  guides(color=guide_legend(title="SCORE2 model", reverse=TRUE)) +
  theme_bw() +
  theme(
    axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(),
    axis.text.x=element_text(size=6), axis.title.x=element_text(size=8),
    panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(),
    strip.background=element_blank(), strip.text=element_text(size=8, face="bold"), 
    legend.text=element_text(size=6), legend.title=element_text(size=8)
  )
ggsave(g, width=7.2, height=1.5, file="analyses/test/cindex_by_SCORE2_method.pdf")


