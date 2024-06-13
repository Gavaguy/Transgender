require(vegan)
require(OTUtable)
require(magrittr)
require(pals)
require(pairwiseAdonis)
require(mgcv)

# Load the metadata
Data <- read.csv(file = "path/to/meta_data.csv")

# Standardize column names
colnames(x = Data) <- tolower(x = colnames(Data))

# Load the pathways data
metab <- read.csv(file = "path/to/pathways.csv", row.names = 1)

# Attach metadata for easy access
attach(what = Data, warn.conflicts = FALSE)

# Compute Bray-Curtis distance
bray_dist <- dist(x = t(metab))

# Perform beta dispersion analysis
mod <- betadisper(bray_dist, paste0(time, transition))
anova(mod)

# Perform pairwise PERMANOVA
set.seed(123)
pairwise.adonis2(bray_dist ~ time * transition, data = Data, strata = "couple", nperm = 9999)

# Perform PCoA
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 3)  # k = 2 for 2D biplot
eigvals <- pcoa$eig
total_variance <- sum(eigvals)

# Calculate the variance explained by PCoA 1 and PCoA 2
var_explained_pcoa1 <- (eigvals[1] / total_variance) * 100
var_explained_pcoa2 <- (eigvals[2] / total_variance) * 100

# biplot

par(ps = 10)

plot(0, type = 'n', axes = FALSE, ann = FALSE)

plot(0, 0, xlim = range(pcoa$points[, 1]), ylim = range(pcoa$points[, 2]), type = "n",
     bty = "l",
     xlab = paste("PC1 (", round(var_explained_pcoa1, 2), "%)", sep = ""),
     ylab = paste("PC2 (", round(var_explained_pcoa2, 2), "%)", sep = ""))

abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

ordiellipse(ord = pcoa, groups = paste(time, transition), col = brewer.paired(n = 4), draw = "polygon", border = 0)

points(pcoa$points[, 1], pcoa$points[, 2], pch = 16, col = brewer.paired(n = 4), cex = 1.2)

legend("bottomright", legend = unique(paste(time, transition)), pch = 21, pt.bg = brewer.paired(n = 4))

aov(pcoa$points[, 1] ~ time * transition + Error(couple / time), data = Data) %>%
  summary()

# Paired t-tests
Before <- which(transition %in% "MTF" & time %in% "BL")
After <- which(transition %in% "MTF" & time %in% "M3")
t.test(pcoa$points[, 1][Before], pcoa$points[, 1][After], paired = TRUE)

Before <- which(transition %in% "FTM" & time %in% "BL")
After <- which(transition %in% "FTM" & time %in% "M3")
t.test(pcoa$points[, 1][Before], pcoa$points[, 1][After], paired = TRUE)

par(mfrow = c(1, 2), xpd = TRUE, ps = 10)

# Boxplot for PC1
boxplot(pcoa$points[, 1] ~ time * transition, outline = FALSE,
        col = brewer.paired(n = 4), xlab = NA, notch = TRUE,
        ylab = paste("PC1 (", round(var_explained_pcoa1, 2), "%)", sep = ""), frame = FALSE)

segments(x0 = 1, x1 = 2, y0 = 0.0016)
text(x = 1.5, y = 0.00165, labels = "**")
segments(x0 = 2, x1 = 3, y0 = 0.00185)
text(x = 2.5, y = 0.00195, labels = "ns")
segments(x0 = 3, x1 = 4, y0 = 0.0018)
text(x = 3.5, y = 0.0019, labels = "(*)")
segments(x0 = 1, x1 = 4, y0 = -0.00125)
text(x = 2.5, y = -0.00115, labels = "ns")

# Boxplot for PC2
boxplot(pcoa$points[, 2] ~ time * transition, outline = FALSE,
        col = brewer.paired(n = 4), xlab = NA, notch = TRUE,
        ylab = paste("PC2 (", round(var_explained_pcoa2, 2), "%)", sep = ""), frame = FALSE)

# dev.off()
