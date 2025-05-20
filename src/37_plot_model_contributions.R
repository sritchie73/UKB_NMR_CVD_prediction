library(data.table)
library(ggplot2)

res <- fread("analyses/CVD_weight_training/multivariable_model_weights.txt")

# Filter to subset of models to plot
res <- res[model %in% c("SCORE2 + PRS", "SCORE2 + NMR scores", "SCORE2 + Biochemistry")]

# Codify factors for plot ordering
res[, sex := factor(sex, levels=c("Male", "Female"))]
res[, model := factor(model, levels=c("SCORE2 + Biochemistry", "SCORE2 + NMR scores", "SCORE2 + PRS"))]
res[, variable_name := factor(variable_name, levels=res[,.(Z=mean(qnorm(HR.pval))), by=.(variable_name, model)][order(-abs(Z)), variable_name])]

# Categorize significance
res[, HR.pval.cat := ifelse(HR.pval < 0.05, "P < 0.05", "P > 0.05")]
res[, HR.pval.cat := factor(HR.pval.cat, levels=c("P < 0.05", "P > 0.05"))]

# Make plot
g <- ggplot(res) + 
  aes(x=variable_name, y=HR, ymin=HR.L95, ymax=HR.U95, color=sex, fill=HR.pval.cat) +
  facet_grid(~ model, scales="free_x", space="free_x") +
  geom_hline(yintercept=1, linetype=2) + 
  geom_errorbar(width=0, position=position_dodge(width=0.3)) +
  geom_point(shape=23, position=position_dodge(width=0.3)) +
  scale_color_manual(values=c("Male"="#d53e4f", "Female"="#3288bd")) +
  scale_fill_manual(values=c("P < 0.05"="white", "P > 0.05"="#d53e4f")) +
  ylab("Hazard Ratio (95% CI)") +
  theme_bw() +
  theme(
    axis.title.y=element_text(size=8), axis.text.y=element_text(size=6),
    axis.title.x=element_blank(), axis.text.x=element_text(size=8, color="black", angle=90, vjust=0.5, hjust=1),
    panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank(),
    strip.background=element_blank(), strip.text=element_blank(),
    legend.position="bottom", legend.title=element_blank(), legend.text=element_text(size=7)
  )
ggsave(g, width=7.2, height=4, file="analyses/CVD_weight_training/multivariable_model_weights.pdf")



