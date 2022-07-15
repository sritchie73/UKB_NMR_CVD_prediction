library(data.table)
library(ggplot2)

# Load data C-indices
cinds <- fread("analyses/test/C_indices.txt")
nri <- fread("analyses/test/nri_estimates.txt")
model_info <- fread("analyses/test/model_fit_information.txt")

# Make nice for plotting
cinds <- cinds[order(C.index)]
cinds[, long_name := factor(long_name, levels=rev(unique(long_name)))]
cinds[, lambda := fcase(
  is.na(lambda) | lambda == "", "No feature selection",
  lambda == "lambda.min", "Best model",
  lambda == "lambda.1se", "Best model with fewest features")]
cinds[, PGS := ifelse(PGS, "With CAD metaGRS + Stroke metaGRS", "Without PGS")]
cinds[, PGS := factor(PGS, levels=c("Without PGS", "With CAD metaGRS + Stroke metaGRS"))]

nri <- nri[metric %in% c("NRI+", "NRI-")]
nri[, compared_to := "Conventional RF"]
nri[, biomarkers := gsub(".*to ", "", model)]
nri[, biomarkers := gsub(" \\(.*\\)", "", biomarkers)]
nri[model_info[lambda == "lambda.min"], on = .(biomarkers=name), long_name := ifelse(lambda == "Best model", i.long_name, NA_character_)]
nri[model_info[lambda == "lambda.1se"], on = .(biomarkers=name), long_name := ifelse(lambda == "Best model with fewest features", i.long_name, x.long_name)]
nri[biomarkers == "PGS", long_name := "Conventional RF + PGS"]
nri[biomarkers == "PGS + Assays" & lambda == "Best model",
      long_name := model_info[lambda == "lambda.min" & (PGS) & name == "Assays", long_name]]
nri[biomarkers == "PGS + Assays" & lambda == "Best model with fewest features",
      long_name := model_info[lambda == "lambda.1se" & (PGS) & name == "Assays", long_name]]
nri[biomarkers == "PGS + NMR" & lambda == "Best model",
      long_name := model_info[lambda == "lambda.min" & (PGS) & name == "NMR", long_name]]
nri[biomarkers == "PGS + NMR" & lambda == "Best model with fewest features",
      long_name := model_info[lambda == "lambda.1se" & (PGS) & name == "NMR", long_name]]
nri[biomarkers == "PGS + NMR + Assays" & lambda == "Best model",
      long_name := model_info[lambda == "lambda.min" & (PGS) & name == "NMR + Assays", long_name]]
nri[biomarkers == "PGS + NMR + Assays" & lambda == "Best model with fewest features",
      long_name := model_info[lambda == "lambda.1se" & (PGS) & name == "NMR + Assays", long_name]]
nri[biomarkers %like% "CRP", long_name := paste("Conventional RF +", biomarkers)]
nri[, nri_group := ifelse(metric == "NRI+", "cases", "controls")]
nri[, biomarkers := factor(biomarkers, levels=rev(c("CRP", "NMR", "Assays", "NMR + Assays",
  "PGS", "PGS + CRP", "PGS + NMR", "PGS + Assays", "PGS + NMR + Assays")))]
nri[, lambda := factor(lambda, levels=rev(c("No feature selection", "Best model with fewest features", "Best model")))]
nri <- nri[order(lambda)][order(biomarkers)]
nri[, long_name := factor(long_name, levels=unique(long_name))]
nri[, PGS := ifelse(biomarkers %like% "PGS", "With CAD metaGRS + Stroke metaGRS", "Without PGS")]
nri[, PGS := factor(PGS, levels=c("Without PGS", "With CAD metaGRS + Stroke metaGRS"))]

# Question 1: do lambda.min models outperform lambda.1se models in the test data? 
g <- ggplot(cinds) +
  aes(x = C.index, xmin = L95, xmax = U95, y = long_name, color=lambda, shape=name) +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_wrap(~ PGS, nrow=2, scales="free_y") +
  scale_shape_manual(name="Biomarkers", values=c("Conventional RF"=17, "CRP"=23, "Assays"=23, "NMR"=24, "NMR + Assays"=22)) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="red", "Best model"="#00B050",
                     "Best model with fewest features"="#0070C0")) +
  labs(x = "C-index (95% confidence interval)", y = "") +
  theme_bw() +
  theme(legend.position="right", legend.box="vertical", axis.text.y=element_text(size=8))
ggsave(g, width=14, height=5, units="in", file="analyses/test/cindex_compare_lambda_min_1se.pdf")

g <- ggplot(nri) +
  aes(x = Estimate*100, xmin=Lower*100, xmax=Upper*100, y=biomarkers, color=lambda, shape=biomarkers) +
  geom_vline(xintercept=0, linetype=2) +
  geom_vline(data=nri[biomarkers == "PGS"], aes(xintercept=Estimate*100), color="red", linetype=2) +
  geom_errorbarh(height=0, position=position_dodge(width=0.6)) +
  geom_point(size=2, fill="white", position=position_dodge(width=0.6)) +
  facet_grid( ~ nri_group, scales="free") +
  scale_shape_manual(name="Biomarkers", values=c(
    "CRP"=23, "Assays"=23, "NMR"=24, "NMR + Assays"=22, "PGS"=17,
    "PGS + CRP"=18, "PGS + Assays"=18, "PGS + NMR"=17, "PGS + NMR + Assays"=15
  )) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="red", "Best model"="#00B050",
                     "Best model with fewest features"="#0070C0")) +
  ylab("") + xlab("Continuous NRI, % reclassified (95% CI)") +
  theme_bw() + theme(legend.position="bottom", legend.box="vertical")
ggsave(g, width=13, height=7, units="in", file="analyses/test/nri_compare_lambda_min_1se.pdf")

# Lambda.min models all outperfrom lambda.1se, so we will proceed with those
cinds <- cinds[lambda != "Best model with fewest features"]
nri <- nri[lambda != "Best model with fewest features"]

# Plot change in C-index relative to conventional risk factors model
conv_rf_cind <- cinds[name == "Conventional RF" & lambda == "No feature selection" & PGS == "Without PGS", C.index]
conv_rf_pgs_cind <- cinds[name == "Conventional RF" & lambda == "No feature selection" & PGS == "With CAD metaGRS + Stroke metaGRS", C.index]
g <- ggplot(cinds) +
  aes(x = deltaC, xmin = deltaC.L95, xmax = deltaC.U95, y = long_name, color=lambda, shape=name) +
  geom_vline(xintercept=0, linetype=2) +
  geom_vline(xintercept=conv_rf_pgs_cind - conv_rf_cind, linetype=2, color="orange") +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_wrap(~ PGS, nrow=2, scales="free_y") +
  scale_shape_manual(name="Biomarkers", values=c("Conventional RF"=17, "CRP"=23, "Assays"=23, "NMR"=24, "NMR + Assays"=22)) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="#d95f02", "Best model"="#7570b3")) +
  scale_x_continuous(name="Change in C-index (95% confidence interval)", sec.axis=sec_axis(~ . + conv_rf_cind, name="C-index (95% confidence interval)")) +
  ylab("") +
  theme_bw() +
  theme(legend.position="right", legend.box="vertical", axis.text.y=element_text(size=8))
ggsave(g, width=10, height=5, units="in", file="analyses/test/delta_cindex_compare.pdf")

g <- ggplot(nri) +
  aes(x = Estimate*100, xmin=Lower*100, xmax=Upper*100, y=long_name, color=lambda, shape=biomarkers) +
  geom_vline(xintercept=0, linetype=2) +
  geom_vline(data=nri[biomarkers == "PGS", .(Estimate, nri_group)], aes(xintercept=Estimate*100), color="red", linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_grid(PGS ~ nri_group, scales="free") +
  scale_shape_manual(name="Biomarkers", values=c(
    "CRP"=23, "Assays"=23, "NMR"=24, "NMR + Assays"=22, "PGS"=17,
    "PGS + CRP"=18, "PGS + Assays"=18, "PGS + NMR"=17, "PGS + NMR + Assays"=15
  )) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="#d95f02", "Best model"="#7570b3")) +
  ylab("") + 
  xlab("Continuous NRI, % reclassified (95% CI)") +
  theme_bw() + 
  theme(legend.position="bottom", legend.box="vertical", axis.text.y=element_text(size=8))
ggsave(g, width=10, height=5, units="in", file="analyses/test/nri_compare.pdf")

# Create table for supp
dt <- cinds[order(PGS), .(long_name, Samples, Cases, C.index, L95, U95, deltaC, deltaC.L95, deltaC.U95)]
dt[nri[nri_group == "cases"], on  = .(long_name), 
  c("NRI.samples", "NRI.cases", "CaseNRI", "CaseNRI.L95", "CaseNRI.U95") := 
  .(i.samples, i.cases, Estimate, Lower, Upper)]
dt[nri[nri_group == "controls"], on  = .(long_name), 
  c("ControlNRI", "ControlNRI.L95", "ControlNRI.U95") := 
  .(Estimate, Lower, Upper)]

fwrite(dt, sep="\t", quote=FALSE, file="analyses/test/cind_nri_compare.txt")


