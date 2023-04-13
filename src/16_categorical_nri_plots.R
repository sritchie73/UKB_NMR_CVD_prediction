library(data.table)
library(ggplot2)

# Load NRI estimates computed from recalibrated risk
nri_estimates <- fread("analyses/public_health_modelling/net_reclassification/nri_estimates.txt")

# Load continuous NRI computed directly from Cox models
nri_before_recalibration <- fread("analyses/test/nri_estimates.txt")

# Compare 
nri_before_recalibration <- rbind(idcol="recalibration", fill=TRUE,
  "After"=nri_estimates[nri_type == "Continuous NRI"],
  "Before"=nri_before_recalibration
)

g <- ggplot(nri_before_recalibration) + 
  aes(x=Estimate, xmin=Lower, xmax=Upper, color=recalibration, y=model) +
  geom_errorbarh(height=0, position=position_dodge(width=0.6)) +
  geom_point(size=2, position=position_dodge(width=0.6)) +
  facet_grid(~ metric, scales="free") +
  scale_color_manual(name="Risk recalibration", values=c("Before"="black", "After"="red")) +
  ylab("") + xlab("Continous NRI metric") +
  theme_bw() + theme(legend.position="bottom")
ggsave(g, width=13, height=7, units="in", file="analyses/public_health_modelling/net_reclassification/continuous_nri_before_and_after_recalibration.pdf")


# Load model information
model_info <- fread("analyses/test/model_fit_information.txt")

# NRI estimates compare how case/control are reclassified when changing from one model 
# (e.g. conventional risk factors alone) to another (e.g. conventional risk factors +
#  selected biomarkers).
#
# We've compute various comparisons, using risk thresholds of 0/5/10% (UK guidelines),
# 0/5/7.5% (American guidelines), and continuous NRI.
gg_dt <- nri_estimates[metric %in% c("NRI+", "NRI-")]
gg_dt[, compared_to := "Conventional RF"]
gg_dt[, biomarkers := gsub(".*to ", "", model)]
gg_dt[, biomarkers := gsub(" \\(.*\\)", "", biomarkers)]
gg_dt[model_info[lambda == "lambda.min"], on = .(biomarkers=name), long_name := ifelse(lambda == "Best model", i.long_name, NA_character_)]
gg_dt[biomarkers == "PGS", long_name := "Conventional RF + PGS"]
gg_dt[biomarkers == "PGS + NMR" & lambda == "Best model", 
      long_name := model_info[lambda == "lambda.min" & (PGS) & name == "NMR", long_name]]
gg_dt[biomarkers == "PGS + NMR" & lambda == "Best model with fewest features",
      long_name := model_info[lambda == "lambda.1se" & (PGS) & name == "NMR", long_name]]
gg_dt[biomarkers %like% "CRP", long_name := paste("Conventional RF +", biomarkers)]
gg_dt[, risk_categories := gsub(" Categorical NRI", "", nri_type)]
gg_dt[, risk_categories := factor(risk_categories, levels=c("Continuous NRI", "ACC/AHA 2019", "NICE 2014"))]
gg_dt[, nri_group := ifelse(metric == "NRI+", "cases", "controls")]
gg_dt[, biomarkers := factor(biomarkers, levels=rev(c("CRP", "NMR", "PGS", "PGS + CRP", "PGS + NMR")))]
gg_dt[, lambda := factor(lambda, levels=rev(c("No feature selection", "Best model")))]
gg_dt <- gg_dt[order(lambda)][order(biomarkers)]
gg_dt[, long_name := factor(long_name, levels=unique(long_name))]

# curate dataset for showing Conventional RF + PGS intercept line
gg_dt2 <- gg_dt[biomarkers == "PGS" & nri_type != "Continuous NRI"]

# plot NRI estimates
g <- ggplot(gg_dt[nri_type != "Continuous NRI"]) +
  aes(x = Estimate*100, xmin=Lower*100, xmax=Upper*100, y=long_name, color=lambda, shape=biomarkers) +
  geom_vline(xintercept=0, linetype=2) +
  geom_vline(data=gg_dt2, aes(xintercept=Estimate*100), color="#d95f02", linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_grid(compared_to ~ risk_categories + nri_group) +
  scale_shape_manual(name="Biomarkers", values=c("CRP"=23, "NMR"=24, "PGS"=17, "PGS + CRP"=18, "PGS + NMR"=17)) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="#d95f02", "Best model"="#7570b3")) +
  ylab("") + xlab("Categorical NRI, % reclassified (95% CI)") +
  theme_bw() + theme(legend.position="bottom")
ggsave(g, width=13, height=7, units="in", file="analyses/public_health_modelling/net_reclassification/nri_compare.pdf")

# Build wide table for supp
wide <- dcast(gg_dt[lambda != "Best model with fewest features"], long_name + biomarkers + compared_to ~ nri_type + metric, value.var=c("Estimate", "Lower", "Upper"))
wide <- wide[order(-biomarkers)][order(compared_to)]
wide <- wide[, .(long_name, biomarkers, compared_to, 
  `ACCAHA_NRI+`=`Estimate_ACC/AHA 2019 Categorical NRI_NRI+`, `ACCAHA_NRI+_L95`=`Lower_ACC/AHA 2019 Categorical NRI_NRI+`, `ACCAHA_NRI+_U95`=`Upper_ACC/AHA 2019 Categorical NRI_NRI+`,
  `ACCAHA_NRI-`=`Estimate_ACC/AHA 2019 Categorical NRI_NRI-`, `ACCAHA_NRI-_L95`=`Lower_ACC/AHA 2019 Categorical NRI_NRI-`, `ACCAHA_NRI-_U95`=`Upper_ACC/AHA 2019 Categorical NRI_NRI-`,
  `NICE_NRI+`=`Estimate_NICE 2014 Categorical NRI_NRI+`, `NICE_NRI+_L95`=`Lower_NICE 2014 Categorical NRI_NRI+`, `NICE_NRI+_U95`=`Upper_NICE 2014 Categorical NRI_NRI+`,
  `NICE_NRI-`=`Estimate_NICE 2014 Categorical NRI_NRI-`, `NICE_NRI-_L95`=`Lower_NICE 2014 Categorical NRI_NRI-`, `NICE_NRI-_U95`=`Upper_NICE 2014 Categorical NRI_NRI-`
)]

# Obtain sample, cases, and reclassification numbers
reclassified <- fread("analyses/public_health_modelling/net_reclassification/nri_reclassified.txt")
reclassified[, compared_to := "Conventional RF"]
reclassified[, biomarkers := gsub(".*to ", "", model)]
reclassified[, biomarkers := gsub(" \\(.*\\)", "", biomarkers)]
reclassified[model_info[lambda == "lambda.min"], on = .(biomarkers=name), long_name := ifelse(lambda == "Best model", i.long_name, NA_character_)]
reclassified[biomarkers == "PGS", long_name := "Conventional RF + PGS"]
reclassified[biomarkers == "PGS + NMR" & lambda == "Best model", 
      long_name := model_info[lambda == "lambda.min" & (PGS) & name == "NMR", long_name]]
reclassified[biomarkers == "PGS + NMR" & lambda == "Best model with fewest features",
      long_name := model_info[lambda == "lambda.1se" & (PGS) & name == "NMR", long_name]]
reclassified[biomarkers %like% "CRP", long_name := paste("Conventional RF +", biomarkers)]
reclassified[, risk_categories := gsub(" Categorical NRI", "", nri_type)]
reclassified[, risk_categories := factor(risk_categories, levels=c("Continuous NRI", "ACC/AHA 2019", "NICE 2014"))]
reclassified[, biomarkers := factor(biomarkers, levels=rev(c("CRP", "NMR", "PGS", "PGS + CRP", "PGS + NMR")))]
reclassified[, lambda := factor(lambda, levels=rev(c("No feature selection", "Best model")))]
reclassified <- reclassified[order(lambda)][order(biomarkers)]
reclassified[, long_name := factor(long_name, levels=unique(long_name))]
reclassified[, Controls := All - Cases]

# Add to wide table
wide[reclassified, on = .(long_name), c("Samples", "Cases") := .(i.Total_Samples, i.Total_Cases)]
wide[reclassified[Old == ">= 0.1", .(Cases=sum(Cases)), by=long_name], on = .(long_name), NICE_Old_Case_Correct := i.Cases]
wide[reclassified[New == ">= 0.1", .(Cases=sum(Cases)), by=long_name], on = .(long_name), NICE_New_Case_Correct := i.Cases]
wide[reclassified[Old == "< 0.1" & New == ">= 0.1", .(Cases=sum(Cases)), by=long_name], on = .(long_name), NICE_Case_Reclassified := i.Cases]
wide[reclassified[Old == "< 0.1", .(Controls=sum(Controls)), by=long_name], on = .(long_name), NICE_Old_Control_Correct := i.Controls]
wide[reclassified[New == "< 0.1", .(Controls=sum(Controls)), by=long_name], on = .(long_name), NICE_New_Control_Correct := i.Controls]
wide[reclassified[Old == ">= 0.1" & New == "< 0.1", .(Controls=sum(Controls)), by=long_name], on = .(long_name), NICE_Control_Reclassified := i.Controls]
wide[reclassified[Old == ">= 0.075", .(Cases=sum(Cases)), by=long_name], on = .(long_name), ACCAHA_Old_Case_Correct := i.Cases]
wide[reclassified[New == ">= 0.075", .(Cases=sum(Cases)), by=long_name], on = .(long_name), ACCAHA_New_Case_Correct := i.Cases]
wide[reclassified[Old == "< 0.075" & New == ">= 0.075", .(Cases=sum(Cases)), by=long_name], on = .(long_name), ACCAHA_Case_Reclassified := i.Cases]
wide[reclassified[Old == "< 0.075", .(Controls=sum(Controls)), by=long_name], on = .(long_name), ACCAHA_Old_Control_Correct := i.Controls]
wide[reclassified[New == "< 0.075", .(Controls=sum(Controls)), by=long_name], on = .(long_name), ACCAHA_New_Control_Correct := i.Controls]
wide[reclassified[Old == ">= 0.075" & New == "< 0.075", .(Controls=sum(Controls)), by=long_name], on = .(long_name), ACCAHA_Control_Reclassified := i.Controls]

# Organise columns
wide <- wide[, .(long_name, biomarkers, compared_to, Samples, Cases,
  NICE_Old_Case_Correct, NICE_New_Case_Correct, NICE_Case_Reclassified, 
  `NICE_NRI+`, `NICE_NRI+_L95`, `NICE_NRI+_U95`,
  NICE_Old_Control_Correct, NICE_New_Control_Correct, NICE_Control_Reclassified, 
  `NICE_NRI-`, `NICE_NRI-_L95`, `NICE_NRI-_U95`,
  ACCAHA_Old_Case_Correct, ACCAHA_New_Case_Correct, ACCAHA_Case_Reclassified, 
  `ACCAHA_NRI+`, `ACCAHA_NRI+_L95`, `ACCAHA_NRI+_U95`,
  ACCAHA_Old_Control_Correct, ACCAHA_New_Control_Correct, ACCAHA_Control_Reclassified, 
  `ACCAHA_NRI-`, `ACCAHA_NRI-_L95`, `ACCAHA_NRI-_U95`)]

fwrite(wide, sep="\t", quote=FALSE, file="analyses/public_health_modelling/net_reclassification/nri_compare.txt")
