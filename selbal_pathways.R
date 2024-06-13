# Load necessary libraries
library(magrittr)
library(selbal)
library(OTUtable)  # Assuming this is required based on your usage

# Read data
Data <- read.csv(file = "path/to/meta_data.csv", row.names = 1)
colnames(Data) <- tolower(colnames(Data))

metab <- read.csv(file = "path/to/pathways.csv", row.names = 1)

# Filter taxa
metab <- filter_taxa(table = metab, abundance = 10^-2, persistence = 50) %>%
  t() %>%
  as.data.frame()

# Attach data (consider avoiding attach due to potential namespace conflicts)
attach(Data, warn.conflicts = FALSE)

# Define Gender factor
Gender <- ifelse(transition == "MTF" & time == "M3", 0,
                 ifelse(transition == "FTM" & time == "BL", 0, 1)) %>% as.factor()

# Perform selbal cross-validation
Bal.cv <- selbal.cv(x = metab, y = Gender, n.fold = 3, n.iter = 1, seed = 123)

# Plotting results
pdf(file = "path/to/A.pdf",
    width = 5.5, height = 5.0)
plot.new()
grid.draw(Bal.cv$global.plot2)
dev.off()

pdf(file = "path/to/B.pdf",
    width = 4.0, height = 4.0)
plot.new()
grid.draw(Bal.cv$ROC.plot)
dev.off()
