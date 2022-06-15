library(data.table)
library(ggplot2)

# Load in C-indices
cinds <- fread("analyses/test/C_indices.txt")

# Compare C-indices for all models - key question here is do the lambda.min models outperform the lambda.1se
# models in the test data?
ggdt <- cinds[order(C.index)]
ggdt[, long_name := factor(long_name, levels=rev(unique(long_name)))]
ggdt[, lambda := fcase(
  is.na(lambda) | lambda == "", "No feature selection",
  lambda == "lambda.min", "Best model",
  lambda == "lambda.1se", "Best model with fewest features")]
ggdt[, PGS := ifelse(PGS, "With CAD metaGRS + Stroke metaGRS", "Without PGS")]
ggdt[, PGS := factor(PGS, levels=c("Without PGS", "With CAD metaGRS + Stroke metaGRS"))]

g <- ggplot(ggdt) +
  aes(x = C.index, xmin = L95, xmax = U95, y = long_name, color=lambda, shape=name) +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_wrap(~ PGS, nrow=2, scales="free_y") +
  scale_shape_manual(name="Biomarkers", values=c("Conventional RF"=17, "CRP"=23, "GlycA"=24, "Blood"=23, "Nightingale"=24, "Blood + Nightingale"=22)) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="red", "Best model"="#00B050",
                     "Best model with fewest features"="#0070C0")) +
  labs(x = "C-index (95% confidence interval)", y = "") +
  theme_bw() +
  theme(legend.position="right", legend.box="vertical", axis.text.y=element_text(size=8))
ggsave(g, width=14, height=5, units="in", file="analyses/test/cindex_compare_lambda_min_1se.pdf")

# Lambda.min models all outperfrom lambda.1se, so we will proceed with those
ggdt <- ggdt[lambda != "Best model with fewest features"]

g <- ggplot(ggdt) +
  aes(x = C.index, xmin = L95, xmax = U95, y = long_name, color=lambda, shape=name) +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_wrap(~ PGS, nrow=2, scales="free_y") +
  scale_shape_manual(name="Biomarkers", values=c("Conventional RF"=17, "CRP"=23, "GlycA"=24, "Blood"=23, "Nightingale"=24, "Blood + Nightingale"=22)) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="#d95f02", "Best model"="#7570b3")) +
  labs(x = "C-index (95% confidence interval)", y = "") +
  theme_bw() +
  theme(legend.position="right", legend.box="vertical", axis.text.y=element_text(size=8))
ggsave(g, width=14, height=5, units="in", file="analyses/test/cindex_compare.pdf")

# Plot change in C-index relative to conventional risk factors model
g <- ggplot(ggdt) +
  aes(x = deltaC, xmin = deltaC.L95, xmax = deltaC.U95, y = long_name, color=lambda, shape=name) +
  geom_vline(xintercept=0, linetype=2) +
  geom_vline(xintercept=ggdt[name == "Conventional RF" & lambda == "No feature selection" & PGS == "With CAD metaGRS + Stroke metaGRS", deltaC], linetype=2, color="orange") +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_wrap(~ PGS, nrow=2, scales="free_y") +
  scale_shape_manual(name="Biomarkers", values=c("Conventional RF"=17, "CRP"=23, "GlycA"=24, "Blood"=23, "Nightingale"=24, "Blood + Nightingale"=22)) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="#d95f02", "Best model"="#7570b3")) +
  labs(x = "Change in C-index (95% confidence interval)", y = "") +
  theme_bw() +
  theme(legend.position="right", legend.box="vertical", axis.text.y=element_text(size=8))
ggsave(g, width=14, height=5, units="in", file="analyses/test/delta_cindex_compare.pdf")

# Also relative to conventional risk factors + PGS
g <- ggplot(ggdt[PGS == "With CAD metaGRS + Stroke metaGRS"]) +
  aes(x = deltaC.PGS, xmin = deltaC.PGS.L95, xmax = deltaC.PGS.U95, y = long_name, color=lambda, shape=name) +
  geom_vline(xintercept=0, linetype=2) +
  geom_errorbarh(height=0) +
  geom_point(size=2, fill="white") +
  facet_wrap(~ PGS, nrow=2, scales="free_y") +
  scale_shape_manual(name="Biomarkers", values=c("Conventional RF"=17, "CRP"=23, "GlycA"=24, "Blood"=23, "Nightingale"=24, "Blood + Nightingale"=22)) +
  scale_color_manual(name="Feature selection", values=c("No feature selection"="#d95f02", "Best model"="#7570b3")) +
  labs(x = "Change in C-index (95% confidence interval)", y = "") +
  theme_bw() +
  theme(legend.position="right", legend.box="vertical", axis.text.y=element_text(size=8))
ggsave(g, width=14, height=5, units="in", file="analyses/test/delta_cindex_from_PGS.pdf")
