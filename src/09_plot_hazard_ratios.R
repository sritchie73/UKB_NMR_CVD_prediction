library(data.table)
library(ggplot2)

# Make output directory
system("mkdir -p analyses/test/hazard_ratios/")

# Load hazard ratios
hrs <- fread("analyses/test/hazard_ratios.txt")

# Load model info
model_info <- fread("analyses/test/model_fit_information.txt")

# Iterate through models and plot hazard ratios for coefficients:
for (midx in model_info[,.I]) {
  this_model <- model_info[midx]

  hr_dt <- hrs[long_name == this_model$long_name]
  hr_dt <- hr_dt[order(-abs(logHR))]
  hr_dt[, coef_name := factor(coef_name, levels=unique(coef_name))]
  hr_dt[coef_type == "Dataset-specific covariate", coef_type := "Covariate"]
  hr_dt[coef_type == "Polygenic Risk Score", coef_type := "PGS"]

  g <- ggplot(hr_dt) +
    aes(x=coef_name, y = HR, ymin = L95, ymax = U95) +
    geom_errorbar(width=0) +
    geom_point(size=2, shape=18) +
    geom_hline(yintercept = 1, linetype = 2) +
    facet_grid(. ~ coef_type, scales="free_x", space="free_x") +
    ylab("Hazard Ratio (95% CI)") +
    xlab("") +
    theme_bw() +
    theme(legend.position="bottom", legend.box="vertical", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

  ggsave(g, width=13, height=5, units="in", file=sprintf("analyses/test/hazard_ratios/%s%s%s.pdf", 
    tolower(gsub(" \\+? ?", "_", this_model[, name])),
    ifelse(this_model[,PGS], "_PGS", ""),
    ifelse(this_model[,lambda] == "", "", paste0("_", this_model[,lambda]))
  ))
}


