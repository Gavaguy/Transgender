require(OTUtable)
require(magrittr)
require(rms)
require(selbal)
require(vegan)

# Clear the workspace
rm(list = ls())

# Load necessary packages
# If 'selbal' is not installed, uncomment and run the line below:
# devtools::install_github(repo = "malucalle/selbal")

# Load metadata
Data <- read.csv(file = "path/to/meta_data.csv")

# Standardize column names
colnames(x = Data) <- tolower(x = colnames(Data))

# Load abundance table
abund_table <- read.csv(file = "path/to/species.csv", row.names = 1)

# Check if samples match
stopifnot(all(Data$sample == colnames(abund_table)))

# Attach metadata for easy access
attach(what = Data, warn.conflicts = FALSE)

# Subset abundance table to match samples
abund_table <- abund_table[, colnames(abund_table) %in% sample]

# Filter taxa based on abundance and persistence thresholds
abund_table <- filter_taxa(table = abund_table, abundance = 10^-2, persistence = 50)

# Standardize abundance table using Hellinger transformation
abund_hell <- decostand(x = t(abund_table), method = "hell")

# Define gender based on transition and time
Gender <- ifelse(transition == "MTF" & time == "M3", 0,
                 ifelse(transition == "FTM" & time == "BL", 0, 1)) %>% as.factor()

# Perform cross-validation for balance selection
Bal.cv <- selbal.cv(x = abund_hell, y = Gender, n.fold = 3, n.iter = 1, seed = 123)

# Save global plot to PDF
pdf(file = "path/to/figure2_B.pdf", width = 5.5, height = 5.0)
plot.new()
grid.draw(Bal.cv$global.plot2)
dev.off()

# Save ROC plot to PDF
pdf(file = "path/to/figure2_C.pdf", width = 4.0, height = 4.0)
plot.new()
grid.draw(Bal.cv$ROC.plot)
dev.off()
