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
  for (model_id in 1:4) {
    for (gender in c("male", "female", "merged")) {

      # Loading the data
      mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_", data_input, ".rds")))

      # Extracting training set
      eids0 <- readRDS("Data/Split/performance_set_eids_0.rds")
      eids1 <- readRDS("Data/Split/performance_set_eids_1.rds")
      mydata_test <- mydata[c(eids0, eids1), ]

      # Extracting selection set
      eids0 <- readRDS("Data/Split/estimation_set_eids_0.rds")
      eids1 <- readRDS("Data/Split/estimation_set_eids_1.rds")
      mydata <- mydata[c(eids0, eids1), ]


      ### PCE

      cox0 <- readRDS(paste0("Results/Cox_models/cox_pce_female_cvd.rds"))
      cox1 <- readRDS(paste0("Results/Cox_models/cox_pce_male_cvd.rds"))

      ## Test

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
      S_pce <- myS


      ### LASSO stability selection

      cox0 <- readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_female.rds"))
      cox1 <- readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_male.rds"))

      ## Test

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
      S_lasso <- myS


      ### Capacity plot

      censored_case <- mydata_test$case
      censored_case[which(mydata_test$time >= 3653)] <- 0

      cap_pce <- GetCapacity(proba_of_event = 1 - S_pce, case = mydata_test$case, N = 100)
      cap_pce_thr <- GetCapacity(proba_of_event = ifelse((1 - S_pce) >= 0.075, yes = 1, no = 0), case = mydata_test$case, N = 100)
      cap_lasso <- GetCapacity(proba_of_event = 1 - S_lasso, case = mydata_test$case, N = 100)

      {
        pdf(paste0("Figures/Capacity_plot_", outcome, "_m", model_id, "_", data_input, "_", gender, ".pdf"),
          width = 8, height = 8
        )
        par(mar = c(5, 5, 1, 1))
        plot(cap_pce$proportion, cap_pce$tpr,
          type = "l", las = 1, cex.lab = 1.5,
          col = "navy", lwd = 2,
          xlab = "Proportion classified 'at-risk'",
          ylab = "Proportion of detected incident cases",
          panel.first = c(
            abline(h = seq(0, 1, by = 0.2), lty = 3, col = "grey"),
            abline(v = seq(0, 1, by = 0.2), lty = 3, col = "grey")
          )
        )
        lines(cap_lasso$proportion, cap_lasso$tpr, type = "l", col = "red", lwd = 2)
        abline(v = cap_pce_thr$proportion[2], col = "navy", lty = 2, lwd = 2)
        abline(h = cap_pce_thr$tpr[2], col = "navy", lty = 2, lwd = 2)
        abline(
          h = cap_lasso$tpr[which.min(abs(cap_lasso$proportion - cap_pce_thr$proportion[2]))],
          col = "red", lty = 2, lwd = 2
        )
        text(
          x = 0.9, y = cap_lasso$tpr[which.min(abs(cap_lasso$proportion - cap_pce_thr$proportion[2]))], col = "red", cex = 1.4,
          pos = 3, labels = paste0("TPR = ", formatC(cap_lasso$tpr[which.min(abs(cap_lasso$proportion - cap_pce_thr$proportion[2]))],
            format = "f", digits = 2
          ))
        )
        text(
          x = 0.9, y = cap_pce_thr$tpr[2], col = "navy", cex = 1.4,
          pos = 1, labels = paste0("TPR = ", formatC(cap_pce_thr$tpr[2],
            format = "f", digits = 2
          ))
        )
        legend("bottomright",
          lty = 1, lwd = 2, col = c("navy", "red"),
          legend = c("PCE", "LASSO"), cex = 1.5, bty = "n"
        )
        dev.off()
      }
    }
  }
}
