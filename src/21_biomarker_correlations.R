library(data.table)
library(ukbnmr)
library(ggplot2)
library(ggrastr)
library(WGCNA)

# Load training data
dat <- fread("data/processed/training/processed_training_data.txt")

# Extract matrix of biomarker data
dat <- dat[,.SD,.SDcols=c("eid", ukbnmr::nmr_info$Biomarker)]
dat <- as.matrix(dat, rownames="eid")

# Set non-finite values to missing (i.e. 0/0 or x/0 ratios)
dat[!is.na(dat) & !is.finite(dat)] <- NA

# Compute correlation coefficients.
dat <- cor(dat, method="pearson")

# Get distance using TOMdist.
# This computes distance between biomarkers as a function of both their
# correlation with each other, and similarity of correlation to all other
# biomarkers
dist <- TOMdist(abs(dat)^6, TOMType="unsigned")
clust <- hclust(as.dist(dist), method="average")
biomarker_order <- colnames(dat)[clust$order]

dat <- dat[biomarker_order, biomarker_order]

# Melt to long for plotting
dat <- as.data.table(dat, keep.rownames="rn")
dat <- melt(dat, id.vars="rn", variable.name="cn", value.name="pearson")
dat[,rn := factor(rn, levels=biomarker_order)]
dat[,cn := factor(cn, levels=biomarker_order)]

g <- ggplot(dat) +
  aes(x=rn, y=cn, fill=pearson) + 
  rasterize(geom_raster(), dpi=1200) +
  scale_x_discrete(name="", expand=expansion(add = 1.3)) +
  scale_y_discrete(name="", expand=expansion(add = 1.3)) +
  scale_fill_gradient2(name="Pearson correlation", low="#313695", mid="white", high="#a50026", limits=c(-1,1)) +
  coord_fixed() + theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
        legend.position="bottom", axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=1), axis.text.x=element_text(size=1, angle=90, hjust=1, vjust=0.5),
        legend.title=element_text(size=7), legend.text=element_text(size=6), legend.key.height=unit(0.4, "cm"),
        plot.background=element_rect(size=0.5)
  ) + guides(fill = guide_colorbar(title.position="top", title.hjust=0.5))
ggsave(width=3.6, height=4, units="in", file="analyses/train/NMR_biomarker_correlations.pdf")

# Do the same for the test dataset (in the same order)
dat <- fread("data/processed/test/processed_test_data.txt")
dat <- dat[(complete_data)]

# Extract matrix of biomarker data
dat <- dat[,.SD,.SDcols=c("eid", ukbnmr::nmr_info$Biomarker)]
dat <- as.matrix(dat, rownames="eid")

# Set non-finite values to missing (i.e. 0/0 or x/0 ratios)
dat[!is.na(dat) & !is.finite(dat)] <- NA

# Compute correlation coefficients.
dat <- cor(dat, use="pairwise.complete.obs", method="pearson")

# order the same way as training data
dat <- dat[biomarker_order, biomarker_order]

# Melt to long for plotting
dat <- as.data.table(dat, keep.rownames="rn")
dat <- melt(dat, id.vars="rn", variable.name="cn", value.name="pearson")
dat[,rn := factor(rn, levels=biomarker_order)]
dat[,cn := factor(cn, levels=biomarker_order)]

g <- ggplot(dat) +
  aes(x=rn, y=cn, fill=pearson) +
  rasterize(geom_raster(), dpi=1200) +
  scale_x_discrete(name="", expand=expansion(add = 1.3)) +
  scale_y_discrete(name="", expand=expansion(add = 1.3)) +
  scale_fill_gradient2(name="Pearson correlation", low="#313695", mid="white", high="#a50026", limits=c(-1,1)) +
  coord_fixed() + theme_bw() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="bottom", axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        axis.text.y=element_text(size=1), axis.text.x=element_text(size=1, angle=90, hjust=1, vjust=0.5),
        legend.title=element_text(size=7), legend.text=element_text(size=6), legend.key.height=unit(0.4, "cm"),
        plot.background=element_rect(size=0.5)
  ) + guides(fill = guide_colorbar(title.position="top", title.hjust=0.5))
ggsave(width=3.6, height=4, units="in", file="analyses/test/NMR_biomarker_correlations.pdf")


