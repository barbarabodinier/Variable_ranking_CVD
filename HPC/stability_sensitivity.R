library(focus)
library(openxlsx)

# Reading arguments
args <- commandArgs(trailingOnly = TRUE)
data_input <- as.character(args[1])
outcome <- as.character(args[2])
model_id <- as.numeric(args[3])
gender_id <- as.numeric(args[4])
gender <- c("Female", "Male")[gender_id + 1]
seed_id <- as.numeric(args[5])

# Loading the data
mydata <- readRDS(paste0("../Data/", outcome, "_imputed_MW_", data_input, ".rds"))
mydict <- read.xlsx("../Data/Data_dictionary.xlsx")

# Extracting selection set
if (gender == "Female") {
  eids <- readRDS(paste0("../Data/Split/selection_set_eids_0_", seed_id, ".rds"))
} else {
  eids <- readRDS(paste0("../Data/Split/selection_set_eids_1_", seed_id, ".rds"))
}
mydata <- mydata[as.character(eids), ]

# Model construction
predictors <- mydict[which(mydict[, paste0("Model.", model_id, ".(", gender, ")")] == "X"), 1]
x <- mydata[, predictors]
x <- na.exclude(x)
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

# Extracting selection proportions and stable selection status
out <- cbind(
  selprop = SelectionProportions(stab),
  selected = SelectedVariables(stab)
)

# Saving the output
dir.create("../Results/Sensitivity", showWarnings = FALSE)
saveRDS(out, paste0(
  "../Results/Sensitivity/stability_", tolower(outcome), "_m", model_id, "_",
  data_input, "_", tolower(gender), "_", seed_id, ".rds"
))
