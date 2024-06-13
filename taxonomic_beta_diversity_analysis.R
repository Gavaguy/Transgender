# Load necessary libraries
require(vegan)
require(OTUtable)
require(magrittr)
require(pals)
require(pairwiseAdonis)

# Read metadata
Data <- read.csv(file = "~/path/to/meta_data.csv")

# Convert column names to lowercase
colnames(Data) <- tolower(colnames(Data))

# Read species abundance table
abund_table <- read.csv(file = "~/path/to/species.csv", row.names = 1)

# Check if samples in metadata match columns in abundance table
all(Data$sample == colnames(abund_table))

# Attach metadata
attach(what = Data, warn.conflicts = FALSE)

# Filter taxa based on abundance and persistence
abund_table <- filter_taxa(table = abund_table, abundance = 10^-2, persistence = 50)

# Transpose the abundance table
abund_table <- t(abund_table)

# Standardize abundance data (Hellinger transformation)
abund_table_hell <- decostand(x = abund_table, method = "hell")

# Calculate Bray-Curtis distance
bray_dist <- dist(x = abund_table_hell)

# Perform betadisper analysis
mod <- betadisper(bray_dist, paste(time, transition), type = "centroid")
anova(mod)

# Pairwise Adonis analysis
set.seed(123)
pairwise.adonis2(bray_dist ~ time * transition, data = Data, strata = "couple", nperm = 9999)

# Perform Principal Coordinate Analysis (PCoA)
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2) # k = 2 for 2D biplot
eigvals <- pcoa$eig
total_variance <- sum(eigvals)

# Calculate variance explained by PCoA 1 and PCoA 2
var_explained_pcoa1 <- (eigvals[1] / total_variance) * 100
var_explained_pcoa2 <- (eigvals[2] / total_variance) * 100

# Plotting
par(ps = 10)
plot(0, 0, xlim = range(pcoa$points[, 1]), ylim = range(pcoa$points[, 2]), type = "n",
     bty = "l",
     xlab = paste("PC1 (", round(var_explained_pcoa1, 2), "%)", sep = ""),
     ylab = paste("PC2 (", round(var_explained_pcoa2, 2), "%)", sep = ""))

abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

ordiellipse(ord = pcoa, groups = paste(time, transition), col = brewer.paired(n = 4), draw = "polygon", border = 0)
points(pcoa$points[, 1], pcoa$points[, 2], pch = 16, col = brewer.paired(n = 4), cex = 1.2)

legend("bottomright", legend = unique(paste(time, transition)), pch = 16, col = unique(brewer.paired(n = 4)), pt.cex = 1.2)
