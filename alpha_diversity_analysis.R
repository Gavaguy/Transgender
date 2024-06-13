# Load necessary libraries
require(vegan)
require(OTUtable)
require(beeswarm)
require(magrittr)
require(pals)
require(mgcv)

# Read the metadata
Data <- read.csv(file = "~/my_path/meta_data.csv")

# Convert column names to lowercase
colnames(Data) <- tolower(colnames(Data))

# Read the species abundance table
abund_table <- read.csv(file = "~/my_path/species.csv", row.names = 1)

# Attach metadata
attach(what = Data, warn.conflicts = FALSE)

# Filter taxa based on abundance and persistence
abund_table <- filter_taxa(table = abund_table, abundance = 10^-4, persistence = 2/72)

# Transpose abundance table
abund_table <- t(abund_table)

# Calculate species richness
Data$spec_number <- specnumber(x = abund_table)

# Standardize abundance data (Hellinger transformation)
abund_table_hell <- decostand(x = abund_table, method = "hell")

# Calculate Shannon diversity index
Data$shannon <- diversity(x = abund_table_hell)

# Calculate total reads per sample
Data$total_reads <- rowSums(abund_table)

# Attach updated Data
attach(Data, warn.conflicts = FALSE)

# Define color code for plots
colorCode <- brewer.paired(n = 4)[as.numeric(factor(paste(time, transition)))]

# Perform variance tests for Shannon index
var.test(shannon ~ time, alternative = "two.sided", data = subset(Data, transition == "FTM"))
var.test(shannon ~ time, alternative = "two.sided", data = subset(Data, transition == "MTF"))

# Fit generalized additive model (GAM) for Shannon index
mod <- gam(shannon ~ transition * time + s(factor(x = couple), bs = "re"), data = Data)

# Display summary of GAM
summary(mod)

# Set up plot layout
par(ps = 12, mfrow = c(3, 2))

# Create boxplots and beeswarm plots for Shannon index
boxplot(shannon ~ time * transition, outline = FALSE,
        col = brewer.paired(n = 4), xlab = NA, ylab = "Shannon Index")
beeswarm(shannon ~ time * transition, pch = 21, cex = 1.2, bg = brewer.paired(n = 4), add = TRUE)

dev.off()
