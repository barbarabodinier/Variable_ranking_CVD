rm(list = ls())

library(sharp)
library(survival)
library(openxlsx)

source("Scripts/functions.R")

dir.create("Results/Refitted_models", showWarnings = FALSE)

# Parameters
outcome <- "cvd"
data_input <- "updated"
genders <- c("male", "female")
model_list <- c(1:4, 7)
n_digits <- 2

# Variable annotation
mydict <- read.xlsx("Data/Data_dictionary.xlsx")

# Extracting ever-used variables
mydict <- mydict[which(apply(mydict[, 4:ncol(mydict)], 1, FUN = function(x) {
  any(x == "X")
})), ]
rownames(mydict) <- mydict[, 1]

# Looping over models
for (j in 1:length(genders)) {
  gender <- genders[j]
  print(gender)
  cat("\n")

  for (model_id in model_list) {
    print(paste0("Model ", model_id))

    # Defining input/output name
    full_model_name <- paste0(outcome, "_m", model_id, "_", data_input, "_", gender)


    ## Data preparation

    # Loading the data
    mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_updated.rds")))
    mydata <- cbind(scale(mydata[, rownames(mydict)]), mydata[, c("case", "time")])

    # Using subset with Nightingale data
    if (model_id %in% c(3, 4, 7)) {
      id_column <- which(tolower(colnames(mydict)) == paste0("model.", 3, ".(", gender, ")"))
      predictors <- mydict[which(mydict[, id_column] == "X"), 1]
      z <- mydata[, predictors]
      z <- na.exclude(z)
      mydata <- mydata[rownames(z), ]
    }

    # Extracting training set
    if (gender == "female") {
      eids <- readRDS("Data/Split/estimation_set_eids_0_1.rds")
    } else {
      eids <- readRDS("Data/Split/estimation_set_eids_1_1.rds")
    }
    eids <- intersect(eids, rownames(mydata))
    mydata_training <- mydata[eids, ]

    # Extracting test set
    if (gender == "female") {
      eids <- readRDS("Data/Split/performance_set_eids_0_1.rds")
    } else {
      eids <- readRDS("Data/Split/performance_set_eids_1_1.rds")
    }
    eids <- intersect(eids, rownames(mydata))
    mydata_test <- mydata[eids, ]
    rm(mydata)


    ## Survival models

    # Survival model (PCE only)
    cindex_pce <- RunCoxSelected(
      x = mydata_training,
      x_test = mydata_test,
      selected = "logHR",
      digits = n_digits
    )
    cindex_pce <- c(cindex_pce, list(observed_test = mydata_test$case))
    print(cindex_pce$c_test)
    saveRDS(cindex_pce, paste0(
      "Results/Refitted_models/cox_pce_",
      full_model_name, ".rds"
    ))

    # Extracting stably selected variables
    stab <- readRDS(paste0(
      "Results/HPC_results/stability_",
      full_model_name, ".rds"
    ))
    class(stab) <- "variable_selection"
    selected <- SelectedVariables(stab)
    print(sum(selected))
    selected <- names(selected)[selected == 1]

    # Survival model (stably selected)
    cindex_stable <- RunCoxSelected(
      x = mydata_training,
      x_test = mydata_test,
      selected = selected,
      digits = n_digits
    )
    cindex_stable <- c(cindex_stable, list(observed_test = mydata_test$case))
    print(cindex_stable$c_test)
    saveRDS(cindex_stable, paste0(
      "Results/Refitted_models/cox_lasso_",
      full_model_name, ".rds"
    ))


    ## Logistic models

    # Logistic model (PCE only)
    refitted <- Refit(
      xdata = mydata_training[, "logHR", drop = FALSE],
      ydata = mydata_training$case,
      family = "binomial"
    )
    predicted <- predict.glm(refitted, newdata = mydata_test)
    myroc <- pROC::roc(predictor = predicted, response = mydata_test$case, ci = TRUE)
    logit_pce <- list(
      glm = refitted,
      predicted_test = predicted,
      observed_test = mydata_test$case,
      c_test = ReformatAUC(as.numeric(myroc$ci), digits = n_digits)
    )
    saveRDS(logit_pce, paste0(
      "Results/Refitted_models/logit_pce_",
      full_model_name, ".rds"
    ))

    # Extracting stably selected variables
    stab <- readRDS(paste0(
      "Results/HPC_results/stability_",
      full_model_name, ".rds"
    ))
    class(stab) <- "variable_selection"
    selected <- SelectedVariables(stab)
    print(sum(selected))
    selected <- names(selected)[selected == 1]

    # Logistic model (stably selected)
    refitted <- Refit(
      xdata = mydata_training[, selected, drop = FALSE],
      ydata = mydata_training$case,
      family = "binomial"
    )
    predicted <- predict.glm(refitted, newdata = mydata_test)
    myroc <- pROC::roc(predictor = predicted, response = mydata_test$case, ci = TRUE)
    logit_stable <- list(
      glm = refitted,
      predicted_test = predicted,
      observed_test = mydata_test$case,
      c_test = ReformatAUC(as.numeric(myroc$ci), digits = n_digits)
    )
    saveRDS(logit_stable, paste0(
      "Results/Refitted_models/logit_lasso_",
      full_model_name, ".rds"
    ))
    cat("\n")
  }
}
