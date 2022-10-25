rm(list = ls())

library(data.table)
library(glmnet)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)
library(c060)
library(plotrix)

source("Scripts/functions.R")

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
    for (gender in c("male", "female", "merged")) {

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


      ### Reclassification table

      tmp_original <- 1 - S_pce_prs[mydata_test$case == 1]
      tmp_incremented <- 1 - S_selected[mydata_test$case == 1]
      props <- prop.table(table(tmp_original >= 0.075, tmp_incremented >= 0.075), margin = 1)[rbind(c(1, 2), c(2, 1))]
      myreclassif <- cbind(
        table(tmp_original >= 0.075, tmp_incremented >= 0.075),
        formatC(props * 100, format = "f", digits = 1), c("+", "-")
      )

      tmp_original <- 1 - S_pce_prs[mydata_test$case == 0]
      tmp_incremented <- 1 - S_selected[mydata_test$case == 0]
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
        col.names = FALSE, row.names = FALSE, overwrite = TRUE
      )
    }
  }
}
