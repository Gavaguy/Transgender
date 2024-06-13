# Load necessary libraries
require(mgcv)
require(vegan)
require(beeswarm)
require(magrittr)
require(pals)
require(OTUtable)

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
abund_table <- as.data.frame(t(abund_table))

# Standardize abundance data (Hellinger transformation)
abund_table_hell <- decostand(x = abund_table, method = "hell")

# Calculate library size for each sample
library_size <- rowSums(x = abund_table)

# Fit generalized additive models (GAMs) for each species
mod_list <- lapply(abund_table, function(x) gam(x ~ transition * time + s(factor(x = couple), bs = "re"), family = nb(), offset = log(library_size)) %>% summary)

# Extract beta coefficients, standard errors, and p-values
beta <- lapply(mod_list, function(x) x$p.coeff[4] %>% unname) %>% unlist()
se <- lapply(mod_list, function(x) x$se[4] %>% unname) %>% unlist()
p_values <- lapply(mod_list, function(x) x$pTerms.table[3, 3]) %>% unlist()

# Print significant species based on p-value
cat(paste0(names(p_values)[p_values < 0.05], " (p = ", round(p_values[p_values < 0.05], 4), ")"))
cat(paste0(names(beta)[p_values < 0.05], " (beta = ", round(beta[p_values < 0.05], 2), "±", round(se[p_values < 0.05], 2), ")"))

# Adjust p-values using FDR correction
FDR <- p.adjust(p = p_values, method = "fdr")

# Print significant species based on FDR
cat(paste0(names(FDR)[FDR < 0.05], " (FDR = ", round(FDR[FDR < 0.05], 4), ")"))
cat(paste0(names(beta)[FDR < 0.05], " (beta = ", round(beta[FDR < 0.05], 2), "±", round(se[FDR < 0.05], 2), ")"))

# Loop over columns of the data frame and create boxplots
for (col_name in names(FDR)[FDR < 0.05]) {
  boxplot(abund_table_hell[, col_name] ~ time * transition, outline = FALSE,
          col = brewer.paired(n = 4), xlab = NA, ylab = col_name, notch = TRUE)
  beeswarm(abund_table_hell[, col_name] ~ time * transition, pch = 21, cex = 1.2, bg = brewer.paired(n = 4), add = TRUE)
}
