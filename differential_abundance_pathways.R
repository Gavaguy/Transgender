require(mgcv)
require(vegan)
require(beeswarm)
require(magrittr)
require(pals)
require(writexl)

# Load the metadata
Data <- read.csv(file = "path/to/meta_data.csv")

# Standardize column names
colnames(x = Data) <- tolower(x = colnames(Data))

# Load the pathways data
metab <- read.csv(file = "path/to/pathways.csv", row.names = 1)

# Attach the Data
attach(what = Data, warn.conflicts = FALSE)

# Transpose the pathways data and convert to data frame
metab <- metab %>%
  t() %>%
  as.data.frame()

# Extract pathways names
pathways <- colnames(x = metab)

# Clean column names
colnames(metab) <- sapply(strsplit(x=as.character(pathways), split=": "), `[[`, 2)

# Fit generalized additive models
mod_list <- lapply(metab, function(x) gam(x ~ transition * time + s(factor(x = couple), bs = "re")) %>% summary)

# Extract p-values
p_values <- lapply(mod_list, function(x) x$pTerms.table[3, 3])

# Calculate the proportion of significant p-values
length(which(p_values < 0.05)) / dim(metab)[2]

# Adjust p-values for false discovery rate
FDR <- p.adjust(p = p_values, method = "fdr")

# Extract betas and standard errors
beta <- lapply(mod_list, function(x) x$p.coeff[4] %>% unname) %>% unlist()
se <- lapply(mod_list, function(x) x$se[4] %>% unname) %>% unlist()

# Calculate the proportion of significant FDR values
length(pathways[which(FDR < 0.05)]) / dim(metab)[2]

# Identify highly significant pathways
pathways[which(p_values < 0.001)]
length(pathways[which(p_values < 0.001)])

# Prepare the export data frame
Export <- data.frame(cbind(colnames(metab), beta, se, unlist(p_values), FDR), row.names = NULL)
colnames(Export) <- c("HMP Unified Metabolic Analysis Network microbial metabolic pathways", "Estimate", "Standard error", "P-value", "FDR")

# Order the export data frame by beta values
Export <- Export[order(beta, decreasing = TRUE),]

# Plot boxplots for highly significant pathways
for (col_name in colnames(metab)[which(p_values < 0.001)]) {
  boxplot(metab[, col_name] ~ time * transition, outline = FALSE,
          col = brewer.paired(n = 4), xlab = NA, ylab = col_name, notch = T)  
  beeswarm(metab[, col_name] ~ time * transition, pch = 21, cex = 1.2, bg = brewer.paired(n = 4), add = TRUE)
}

# Save the results to an Excel file
write_xlsx(x = Export, path = "path/to/Supplementary_Table_1.xlsx")
