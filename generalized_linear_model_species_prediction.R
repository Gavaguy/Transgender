require(vegan)
require(OTUtable)
require(magrittr)
require(rms)

# Load metadata
Data <- read.csv(file = "path/to/meta_data.csv")

# Standardize column names
colnames(x = Data) <- tolower(x = colnames(Data))

# Load abundance table
abund_table <- read.csv(file = "path/to/species.csv", row.names = 1)

# Attach metadata for easy access
attach(what = Data, warn.conflicts = FALSE)

# Filter taxa based on abundance and persistence thresholds
abund_table <- filter_taxa(table = abund_table, abundance = 10^-2, persistence = 50)

# Standardize abundance table using Hellinger transformation
abund_hell <- decostand(x = abund_table, method = "hell") %>% t() %>% as.data.frame()

# Attach transformed abundance table
attach(what = abund_hell, warn.conflicts = FALSE)

# Define gender based on transition and time
Gender <- ifelse(transition == "MTF" & time == "M3", 0,
                 ifelse(transition == "FTM" & time == "BL", 0, 1))

# Create model dataframe
Model_df <- data.frame(Gender, Parabacteroides_goldsteinii, `Coprococcus_sp._ART55/1`, Coprococcus_eutactus, Escherichia_coli)

# Set datadist options for rms
dd <- datadist(Model_df)
options(datadist = "dd")

# Fit logistic regression model with restricted cubic splines
m1 <- rms::lrm(Gender ~ rcs(Parabacteroides_goldsteinii, 3) + rcs(`Coprococcus_sp._ART55/1`, 3) + rcs(Coprococcus_eutactus, 3) + rcs(Escherichia_coli, 3),
               data = Model_df, x = TRUE, y = TRUE)

# Evaluate model with penalty trace
pentrace(m1, list(simple = c(0, 1, 2), nonlinear = c(0, 100, 200)))

# Update model with specified penalty
pen_m1 <- update(m1, penalty = list(simple = 2, nonlinear = 200))

# Validate penalized model
val <- rms::validate(pen_m1)

# Calculate optimism-corrected C-index
c_opt_corr <- 0.5 * (val[1, 5] + 1)

# Calibrate penalized model with bootstrapping
cal <- rms::calibrate(pen_m1, B = 200)

# Plot calibration curve
par(ps = 11)
plot(cal)
