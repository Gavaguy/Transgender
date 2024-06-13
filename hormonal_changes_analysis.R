# Load necessary libraries
require(vegan)
require(pals)
require(beeswarm)
require(tidyverse)

# Load the dataset
Dat <- read.csv(file = "~/my_path/Hormone_Sheet1.csv")

# Clean and prepare the data
Dat$Gruppe <- gsub(pattern = " ", replacement = "", x = Dat$Gruppe)
Dat <- subset(x = Dat, Gruppe %in% c("FTMT0", "FTMT2", "MTFT0", "MTFT2"))
Dat$Gruppe <- recode(.x = Dat$Gruppe, `FTMT0` = "BL TM", `FTMT2` = "M3 TM", 
                     `MTFT0` = "BL TW", `MTFT2` = "M3 TW")
Dat$Gruppe <- factor(Dat$Gruppe, levels = c("BL TM", "M3 TM", "BL TW", "M3 TW"))

# Separate metadata and values
Meta <- Dat[, 1:3]
Vals <- Dat[, 4:18]

# Replace "<" values with 0 and convert to numeric
Vals[] <- lapply(Vals, function(x) ifelse(grepl(pattern = "<", x = x), 0, x))
Vals[] <- lapply(Vals, as.numeric)

# Attach datasets
attach(Vals, warn.conflicts = FALSE)
attach(what = Meta, warn.conflicts = FALSE)

# Set up the plot layout and parameters
par(mfrow = c(2, 3), ps = 12)

# Create boxplots and beeswarm plots for each hormone
hormones <- c("Testosteron.nmol.L", "Androstendion.nmol.L", "DHEAS.nmol.L", 
              "DHT.nmol.L", "X17.OH.Progesteron.nmol.L", "Estradiol.nmol.L")

labels <- c("Testosteron [nmol/L]", "Androstendion [nmol/L]", 
            "Dehydroepiandrosterone [nmol/L]", "Dihydrotestosterone [nmol/L]", 
            "17-Hydroxyprogesterone [nmol/L]", "Estradiol [nmol/L]")

titles <- c("Testosteron", "Androstendion", "Dehydroepiandrosterone", 
            "Dihydrotestosterone", "17-Hydroxyprogesterone", "Estradiol")

for (i in 1:length(hormones)) {
  boxplot(as.formula(paste(hormones[i], "~Gruppe")), outline = FALSE, 
          col = brewer.paired(n = 4), xlab = NA, ylab = labels[i])
  beeswarm(as.formula(paste(hormones[i], "~Gruppe")), pch = 21, cex = 1.2, 
           bg = brewer.paired(n = 4), add = TRUE, method = "compactswarm", 
           corral = "gutter")
  title(main = titles[i])
}

# Display session information
sessionInfo()
