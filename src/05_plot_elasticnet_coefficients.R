library(data.table)
library(ggplot2)

# Make output directory
system("mkdir -p analyses/train/hazard_ratios/")

# Load hazard ratios
hrs <- fread("analyses/train/lasso_coefficients.txt")

# Load model info
model_info <- fread("analyses/train/cox_lasso_models.txt")

# Iterate through models and plot hazard ratios for coefficients:
for (midx in model_info[,.I]) {
  this_model <- model_info[midx]
  
  hr_dt <- hrs[this_model, on = .(name, lambda), nomatch=0]
  
  if (nrow(hr_dt) > 0) {
		hr_dt <- hr_dt[order(-abs(beta))]
		hr_dt[, coef_name := factor(coef_name, levels=unique(coef_name))]
		hr_dt[coef_type == "Dataset-specific covariate", coef_type := "Covariate"]

		g <- ggplot(hr_dt) +
			aes(x=coef_name, y = exp(beta)) +
			geom_point(size=2, shape=18) +
			geom_hline(yintercept = 1, linetype = 2) +
			facet_grid(. ~ coef_type, scales="free_x", space="free_x") +
			ylab("Hazard Ratio (95% CI)") +
			xlab("") +
			theme_bw() +
			theme(legend.position="bottom", legend.box="vertical", axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

		ggsave(g, width=13, height=5, units="in", file=sprintf("analyses/train/hazard_ratios/%s%s.pdf", 
			tolower(gsub(" \\+? ?", "_", this_model[, name])),
			ifelse(this_model[,lambda] == "", "", paste0("_", this_model[,lambda]))
		))
  }
}

# Output sheet of coefficients in wide format
hr_wide <- dcast(hrs, coef_name + coef_type ~ name + lambda, value.var="beta")
fwrite(hr_wide, sep="\t", quote=FALSE, file="analyses/train/hazard_ratios/all_coef_wide.txt")


