rm(list = ls())

library(data.table)
library(glmnet)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)
library(c060)
library(plotrix)
library(PredictABEL)

source("Scripts/functions.R")
source("Scripts/reclassification_PredictABEL_function.R")

# outcomes=c("_cad", "_cvd", "_cad_pce", "_cvd_pce")
outcomes <- c("cvd")
genders <- c("male", "female")
mymethod <- "lasso"
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

for (i in 1:length(outcomes)) {
  outcome <- outcomes[i]
  print(outcome)
  for (model_id in 1:2) {
    for (gender in c("male", "female")) {

      # Loading the data
      mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_", data_input, ".rds")))

      # Extracting training set
      eids0 <- readRDS("Data/Split/performance_set_eids_0_1.rds")
      eids1 <- readRDS("Data/Split/performance_set_eids_1_1.rds")
      mydata_test <- mydata[c(eids0, eids1), ]

      # Extracting selection set
      eids0 <- readRDS("Data/Split/estimation_set_eids_0_1.rds")
      eids1 <- readRDS("Data/Split/estimation_set_eids_1_1.rds")
      mydata <- mydata[c(eids0, eids1), ]


      ### PCE

      cox0 <- readRDS(paste0("Results/Cox_models/cox_pce_female_cvd.rds"))
      cox1 <- readRDS(paste0("Results/Cox_models/cox_pce_male_cvd.rds"))

      # Get probabilities
      if (gender == "merged") {
        myS <- c(cox0$S_test, cox1$S_test)
      }
      if (gender == "female") {
        myS <- c(cox0$S_test)
      }
      if (gender == "male") {
        myS <- c(cox1$S_test)
      }
      mydata_test <- mydata_test[names(myS), ]
      print(all(names(myS) == rownames(mydata_test)))
      S_pce_prs <- myS


      ### LASSO stability selection

      cox0 <- readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_female.rds"))
      cox1 <- readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_male.rds"))

      # Get probabilities
      if (gender == "merged") {
        myS <- c(cox0$S_test, cox1$S_test)
      }
      if (gender == "female") {
        myS <- c(cox0$S_test)
      }
      if (gender == "male") {
        myS <- c(cox1$S_test)
      }
      mydata_test <- mydata_test[names(myS), ]
      print(all(names(myS) == rownames(mydata_test)))
      S_selected <- myS

      # Reclassification indices (7.5%)
      nri075 <- reclassification(mydata_test[, "case", drop = FALSE],
        cOutcome = 1, cutoff = c(0, 0.075, 1),
        predrisk1 = 1 - S_pce_prs, predrisk2 = 1 - S_selected
      )
      probs <- improveProb(1 - S_pce_prs, 1 - S_selected, y = mydata_test$case)
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
        y = mydata_test$case
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
        paste0("Cases (N=", formatC(sum(mydata_test$case == 1), big.mark = ","), ")"),
        paste0("Non-cases (N=", formatC(sum(mydata_test$case == 0), big.mark = ","), ")"),
        paste0("Full population (N=", formatC(nrow(mydata_test), big.mark = ","), ")")
      )
      write.xlsx(as.data.frame(myreclasstable), paste0("Tables/Reclassification_indices_", outcome, "_m", model_id, "_", data_input, "_", gender, ".xlsx"),
        colNames = TRUE, rowNames = TRUE, overwrite = TRUE
      )
    }
  }
}
