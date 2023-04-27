library(data.table)
library(ggplot2)

# Load c-indices
cind <- fread("analyses/test/univariate/C_indices.txt")

# Build table of baseline cindices for each panel
conv_rf_c <- fread("analyses/test/C_indices.txt")
ref <- conv_rf_c[!(PGS) & name == "Conventional RF"]
ref <- ref[,.(type=c("PRS", "assays", "NMR"), C.index)]

ref <- rbind(cind[type == "demographics" & name == "age+sex", .(type = "risk_factors", C.index)], ref)

# Partition and order
cind <- cind[order(C.index)]
cind[, display_name := sprintf("%s%s", display_name, c("", ".2")[1:.N]), by=display_name]
cind[, display_name := factor(display_name, levels=unique(display_name))]

# Build plot
g <- ggplot(cind, aes(x=C.index, xmin=L95, xmax=U95, y=display_name)) +
  geom_vline(data=ref, aes(xintercept=C.index), linetype=2, colour="red") +
  geom_point(shape=18, size=1) +
  geom_errorbarh(height=0, alpha=0.8) +
  xlab("C-index (95% CI)") + ylab("") +
  facet_grid(type ~ ., space="free", scales="free") + 
  theme_bw() +
  theme(
    axis.text=element_text(size=6), 
    axis.title=element_text(size=7),
    strip.text=element_text(size=0),
    strip.background=element_blank(),
    panel.grid.major.y=element_blank()
  )

ggsave(g, width=2, height=28.8, file="analyses/test/univariate/cindex_compare.pdf")

# Collate table
cind[, type := factor(type, levels=c("demographics", "risk_factors", "PRS", "assays", "NMR"))]
cind <- cind[order(type)]
cind <- cind[name != "PRS"]
cind <- rbind(cind[type != "assays" & type != "NMR"], cind[type == "assays" | type == "NMR"][order(-C.index)][order(type)])

cind[, strata := ifelse(name %in% c("sex", "age"), "None", "Sex")]
cind[, display_name := gsub("\\.2", "", display_name)]
cind[, adjusted_for := fcase(
  type == "demographics", "",
  type == "risk_factors", "Age",
  default = "Conventional risk factors")]
cind[, reference_model := fcase(
  type == "demographics", "",
  type == "risk_factors", "Age + sex",
  default = "Conventional risk factors")]

hrs <- fread("analyses/test/univariate/hazard_ratios.txt")
cind[hrs[name == var | name == "age+sex"], on = .(name), c("HR", "HR.L95", "HR.U95", "Pvalue") := .(i.HR, i.L95, i.U95, i.Pvalue)]

cind <- cind[,.(group=type, strata, adjusted_for, variable=display_name, Samples, Cases, 
  HR, HR.L95, HR.U95, Pvalue, C.index, L95, U95, reference_model, deltaC, deltaC.L95, deltaC.U95)]

fwrite(cind, sep="\t", file="analyses/test/univariate/cind_compare.txt")


  




