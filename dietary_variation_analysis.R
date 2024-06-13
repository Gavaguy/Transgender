# Load necessary libraries
library(dplyr)

# Read the dataset
food <- read.csv(file = "~/my_path/food.csv", row.names = 1)

# Attach the dataset
attach(what = food, warn.conflicts = FALSE)

# Perform repeated measures ANOVA and calculate p-values for each column
# Reference: https://www.statology.org/repeated-measures-anova-in-r/
p_value <- lapply(food[,5:44], function(x) {
  summary(aov(rank(x) ~ factor(time) * 
                factor(transition) + Error(factor(ID))))$`Error: Within`[[1]][5][2,1]
})

# Filter p-values less than 0.05
significant_p_values <- p_value[p_value < 0.05]
print(significant_p_values)

# Adjust p-values using the false discovery rate (FDR) method
fdr <- p.adjust(p = p_value, method = "fdr")

# Filter adjusted p-values less than 0.1
significant_fdr <- fdr[fdr < 0.1]
print(significant_fdr)

# Display the session information
sessionInfo()
