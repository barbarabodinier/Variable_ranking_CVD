rm(list = ls())

library(survival)
library(openxlsx)
library(Hmisc)
library(PredictABEL)

source("Scripts/functions.R")
source("Scripts/reclassification_PredictABEL_function.R")

# Parameters
outcome <- "cvd"
genders <- c("male", "female")
data_input <- "updated"
model_id <- 1

# Variable annotation
mydict <- read.xlsx("Data/Data_dictionary.xlsx")

# Extracting ever-used variables
mydict <- mydict[which(apply(mydict[, 4:ncol(mydict)], 1, FUN = function(x) {
  any(x == "X")
})), ]
rownames(mydict) <- mydict[, 1]
variable_cat <- mydict[, 3]
names(variable_cat) <- mydict[, 2]

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

    # Reclassification indices (7.5%)
    nri075 <- reclassification(cbind(observed_test),
      cOutcome = 1, cutoff = c(0, 0.075, 1),
      predrisk1 = 1 - S_pce_prs, predrisk2 = 1 - S_selected
    )
    probs <- improveProb(1 - S_pce_prs, 1 - S_selected, y = observed_test)
    mycontnri_10 <- c(
      paste0(
        formatC(probs$nri.ev, format = "f", digits = 3), " (",
        formatC(probs$nri.ev - 1.96 * probs$se.nri.ev, format = "f", digits = 3), ";",
        formatC(probs$nri.ev + 1.96 * probs$se.nri.ev, format = "f", digits = 3), ")"
      ),
      paste0(
        formatC(probs$nri.ne, format = "f", digits = 3), " (",
        formatC(probs$nri.ne - 1.96 * probs$se.nri.ne, format = "f", digits = 3), ";",
        formatC(probs$nri.ne + 1.96 * probs$se.nri.ne, format = "f", digits = 3), ")"
      ),
      paste0(
        formatC(probs$nri, format = "f", digits = 3), " (",
        formatC(probs$nri - 1.96 * probs$se.nri, format = "f", digits = 3), ";",
        formatC(probs$nri + 1.96 * probs$se.nri, format = "f", digits = 3), ")"
      )
    )
    c1 <- cut(1 - S_pce_prs, breaks = c(0, 0.075, 1), include.lowest = TRUE, right = FALSE)
    c2 <- cut(1 - S_selected, breaks = c(0, 0.075, 1), include.lowest = TRUE, right = FALSE)
    c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
    c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
    probs <- improveProb(
      x1 = as.numeric(c11) * (1 / (length(levels(c11)))),
      x2 = as.numeric(c22) * (1 / (length(levels(c22)))),
      y = observed_test
    )
    mycatnri_10 <- c(
      paste0(
        formatC(probs$nri.ev, format = "f", digits = 3), " (",
        formatC(probs$nri.ev - 1.96 * probs$se.nri.ev, format = "f", digits = 3), ";",
        formatC(probs$nri.ev + 1.96 * probs$se.nri.ev, format = "f", digits = 3), ")"
      ),
      paste0(
        formatC(probs$nri.ne, format = "f", digits = 3), " (",
        formatC(probs$nri.ne - 1.96 * probs$se.nri.ne, format = "f", digits = 3), ";",
        formatC(probs$nri.ne + 1.96 * probs$se.nri.ne, format = "f", digits = 3), ")"
      ),
      paste0(
        formatC(probs$nri, format = "f", digits = 3), " (",
        formatC(probs$nri - 1.96 * probs$se.nri, format = "f", digits = 3), ";",
        formatC(probs$nri + 1.96 * probs$se.nri, format = "f", digits = 3), ")"
      )
    )
    myidi_10 <- paste0(
      formatC(nri075$idi, format = "f", digits = 4), " (",
      formatC(nri075$idi - 1.96 * nri075$idi_se, format = "f", digits = 4), "-",
      formatC(nri075$idi + 1.96 * nri075$idi_se, format = "f", digits = 4), ")"
    )
    myreclasstable <- cbind(mycontnri_10, mycatnri_10, c(rep("", 2), myidi_10))
    colnames(myreclasstable) <- c("Continuous NRI", "Categorical NRI", "IDI")
    rownames(myreclasstable) <- c(
      paste0("Cases (N=", formatC(sum(observed_test == 1), big.mark = ","), ")"),
      paste0("Non-cases (N=", formatC(sum(observed_test == 0), big.mark = ","), ")"),
      paste0("Full population (N=", formatC(length(observed_test), big.mark = ","), ")")
    )
    write.xlsx(as.data.frame(myreclasstable), paste0("Tables/Reclassification_indices_", outcome, "_m", model_id, "_", data_input, "_", gender, ".xlsx"),
      colNames = TRUE, rowNames = TRUE, overwrite = TRUE
    )
  }
}
