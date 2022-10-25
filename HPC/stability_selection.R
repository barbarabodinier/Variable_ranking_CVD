library(sharp)
library(openxlsx)
print(packageVersion("sharp"))

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
data_input <- as.character(args[1])
outcome <- as.character(args[2])
model_id <- as.numeric(args[3])
gender_id <- as.numeric(args[4])
gender <- c("Female", "Male")[gender_id + 1]
print(model_id)

# Loading the data
mydata <- readRDS(paste0("../Data/", outcome, "_imputed_MW_", data_input, ".rds"))
mydict <- read.xlsx("../Data/Data_dictionary.xlsx")

# Extracting selection set
if (gender == "Female") {
  eids <- readRDS("../Data/Split/selection_set_eids_0_1.rds")
} else {
  eids <- readRDS("../Data/Split/selection_set_eids_1_1.rds")
}
mydata <- mydata[as.character(eids), ]

# Using subset with Nightingale data
if (model_id == 7) {
  predictors <- mydict[which(mydict[, paste0("Model.", 3, ".(", gender, ")")] == "X"), 1]
  z <- mydata[, predictors]
  z <- na.exclude(z)
}

# Model construction
if (model_id == 7) {
  predictors <- mydict[which(mydict[, paste0("Model.", 1, ".(", gender, ")")] == "X"), 1]
} else {
  predictors <- mydict[which(mydict[, paste0("Model.", model_id, ".(", gender, ")")] == "X"), 1]
}
x <- mydata[, predictors]
x <- na.exclude(x)
if (model_id == 7) {
  x <- x[rownames(z), ]
}
mydata <- mydata[rownames(x), ]
y <- mydata[, c("time", "case")]
print(table(y[, 2]))

# Running stability selection
Lambda <- LambdaSequence(lmax = 0.05, lmin = 1e-4, cardinal = 100)
system.time({
  stab <- VariableSelection(
    xdata = x, ydata = y, family = "cox", K = 1000,
    pi_list = seq(0.55, 0.95, by = 0.01), Lambda = Lambda
  )
})

# Saving the output
dir.create("../Results", showWarnings = FALSE)
saveRDS(stab, paste0("../Results/stability_", tolower(outcome), "_m", model_id, "_", data_input, "_", tolower(gender), ".rds"))
