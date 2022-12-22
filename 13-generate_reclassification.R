rm(list = ls())

library(survival)
library(openxlsx)

source("Scripts/functions.R")

# Parameters
outcome <- "cvd"
genders <- c("male", "female")
data_input <- "updated"

# Variable annotation
mydict <- read.xlsx("Data/Data_dictionary.xlsx")

# Extracting ever-used variables
mydict <- mydict[which(apply(mydict[, 4:ncol(mydict)], 1, FUN = function(x) {
  any(x == "X")
})), ]
rownames(mydict) <- mydict[, 1]
variable_cat <- mydict[, 3]
names(variable_cat) <- mydict[, 2]

# Looping over models
for (model_id in 1:2) {
  for (gender in c("male", "female")) {
    print(gender)

    # Defining input/output name
    full_model_name <- paste0(outcome, "_m", model_id, "_", data_input, "_", gender)

    # PCE
    cox <- readRDS(paste0("Results/Refitted_models/cox_pce_", full_model_name, ".rds"))
    S_pce_prs <- cox$S_test
    observed_test <- cox$observed_test

    # LASSO stability selection
    cox <- readRDS(paste0("Results/Refitted_models/cox_lasso_", full_model_name, ".rds"))
    S_selected <- cox$S_test

    # Reclassification table
    tmp_original <- 1 - S_pce_prs[which(observed_test == 1)]
    tmp_incremented <- 1 - S_selected[which(observed_test == 1)]
    props <- prop.table(table(tmp_original >= 0.075, tmp_incremented >= 0.075), margin = 1)[rbind(c(1, 2), c(2, 1))]
    myreclassif <- cbind(
      table(tmp_original >= 0.075, tmp_incremented >= 0.075),
      formatC(props * 100, format = "f", digits = 1), c("+", "-")
    )

    tmp_original <- 1 - S_pce_prs[which(observed_test == 0)]
    tmp_incremented <- 1 - S_selected[which(observed_test == 0)]
    props <- prop.table(table(tmp_original >= 0.075, tmp_incremented >= 0.075), margin = 1)[rbind(c(1, 2), c(2, 1))]
    myreclassif <- rbind(
      myreclassif,
      cbind(
        table(tmp_original >= 0.075, tmp_incremented >= 0.075),
        formatC(props * 100, format = "f", digits = 1), c("-", "+")
      )
    )
    colnames(myreclassif) <- c("<7.5", ">=7.5", "% Reclassified", "Improved Classification")
    myreclassif <- cbind(
      c("PCE", rep(NA, 3)), c("Cases", "Cases", "Non-cases", "Non-cases"),
      rep(c("<7.5", ">=7.5"), 2), myreclassif
    )
    myreclassif <- rbind(
      c(NA, NA, NA, "Selected", rep(NA, ncol(myreclassif) - 4)),
      colnames(myreclassif), myreclassif
    )
    write.xlsx(as.data.frame(myreclassif), paste0("Tables/Reclassification_table_", outcome, "_m", model_id, "_", data_input, "_", gender, ".xlsx"),
      colNames = FALSE, rowNames = FALSE, overwrite = TRUE
    )
  }
}
