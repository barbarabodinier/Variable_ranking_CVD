rm(list = ls())

library(survival)
library(openxlsx)
library(colorspace)
library(pROC)

source("Scripts/functions.R")
source("Scripts/roc_functions.R")

# Parameters
outcome <- "cvd"
model_list <- c(1, 3, 7)
data_input <- "updated"

# Looping over models
for (model_id in model_list) {
  print(model_id)

  if (model_id == 3) {
    legend_letters <- c("C", "D")
  } else {
    legend_letters <- c("A", "B")
  }
  names(legend_letters) <- c("male", "female")

  pdf(paste0("Figures/ROC_curve_", outcome, "_m", model_id, "_", data_input, ".pdf"),
    width = 12, height = 6
  )
  par(mar = c(5, 5, 1, 1), mfrow = c(1, 2))
  for (gender in c("male", "female")) {
    print(gender)

    # Defining input/output name
    full_model_name <- paste0(outcome, "_m", model_id, "_", data_input)

    # Loading the data
    mydata <- data.frame(readRDS("Data/CVD_imputed_MW_updated.rds"))

    eids_test <- c(
      as.character(readRDS(paste0("Data/Split/performance_set_eids_", 0, "_1.rds"))),
      as.character(readRDS(paste0("Data/Split/performance_set_eids_", 1, "_1.rds")))
    )
    mydata_test <- mydata[eids_test, ]

    # Loading the models (stability)
    logit_stable0 <- readRDS(paste0("Results/Refitted_models/logit_lasso_", full_model_name, "_female.rds"))
    logit_stable1 <- readRDS(paste0("Results/Refitted_models/logit_lasso_", full_model_name, "_male.rds"))

    if (gender == "merged") {
      S_lasso <- c(logit_stable0$predicted_test, logit_stable1$predicted_test)
    }
    if (gender == "female") {
      S_lasso <- c(logit_stable0$predicted_test)
    }
    if (gender == "male") {
      S_lasso <- c(logit_stable1$predicted_test)
    }
    mydata_test <- mydata_test[names(S_lasso), ]

    # Loading the models (PCE)
    logit_pce0 <- readRDS(paste0("Results/Refitted_models/logit_pce_", full_model_name, "_female.rds"))
    logit_pce1 <- readRDS(paste0("Results/Refitted_models/logit_pce_", full_model_name, "_male.rds"))

    if (gender == "merged") {
      S_pce <- c(logit_pce0$predicted_test, logit_pce1$predicted_test)
    }
    if (gender == "female") {
      S_pce <- c(logit_pce0$predicted_test)
    }
    if (gender == "male") {
      S_pce <- c(logit_pce1$predicted_test)
    }
    print(all(names(S_pce) == rownames(mydata_test)))
    print(all(names(S_lasso) == rownames(mydata_test)))

    # Computing sensitivity, specificity
    roc_pce <- pROC::roc(predictor = S_pce, response = mydata_test$case, ci = TRUE)
    roc_lasso <- pROC::roc(predictor = S_lasso, response = mydata_test$case, ci = TRUE)

    # Computing AUC
    auc_pce <- ReformatAUC(as.numeric(roc_pce$ci), digits = 2)
    auc_lasso <- ReformatAUC(as.numeric(roc_lasso$ci), digits = 2)

    # ROC curve
    plot(1 - roc_pce$specificities, roc_pce$sensitivities,
      type = "l", las = 1, cex.lab = 1.5,
      xlab = "False Positive Rate", ylab = "True Positive Rate",
      col = "navy", lwd = 1.5
    )
    mtext(legend_letters[gender], side = 2, line = 2.5, at = 1, cex = 4, las = 1)
    lines(1 - roc_lasso$specificities, roc_lasso$sensitivities,
      type = "l", las = 1, cex.lab = 1.5,
      col = "tomato", lwd = 1.5, lty = 1
    )
    abline(0, 1, lty = 3)
    legend("bottomright",
      lty = 1, col = c("navy", "tomato", "forestgreen"),
      lwd = 1.5, bty = "n", cex = 1,
      legend = c(
        paste0("PCE: AUC = ", auc_pce),
        paste0("LASSO Stability Selection: AUC = ", auc_lasso)
      )
    )
  }
  dev.off()
}
