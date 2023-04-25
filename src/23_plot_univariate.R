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


