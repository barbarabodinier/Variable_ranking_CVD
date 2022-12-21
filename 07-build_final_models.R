rm(list = ls())

library(fake)
library(sharp)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)
library(plotrix)
library(colorspace)
source("Scripts/functions.R")

dir.create("Results/Cox_models", showWarnings = FALSE)

genders <- c("male", "female")
outcome_names <- c("cvd")
data_input <- "updated"
penalised <- FALSE

# Variable annotation
mydict <- read.xlsx("Data/Data_dictionary.xlsx")

# Extracting ever-used variables
mydict <- mydict[which(apply(mydict[, 4:ncol(mydict)], 1, FUN = function(x) {
  any(x == "X")
})), ]
rownames(mydict) <- mydict[, 1]

for (i in 1:length(outcome_names)) {
  outcome <- outcome_names[i]
  print(outcome)
  for (j in 1:length(genders)) {
    gender <- genders[j]
    print(gender)
    
    # Loading the data
    mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_updated.rds")))
    mydata <- cbind(scale(mydata[, rownames(mydict)]), mydata[, c("case", "time")])
    
    # Extracting test set
    if (gender == "female") {
      eids <- readRDS("Data/Split/performance_set_eids_0_1.rds")
    } else {
      eids <- readRDS("Data/Split/performance_set_eids_1_1.rds")
    }
    mydata_test <- mydata[eids, ]
    
    # Extracting training set
    if (gender == "female") {
      eids <- readRDS("Data/Split/estimation_set_eids_0_1.rds")
    } else {
      eids <- readRDS("Data/Split/estimation_set_eids_1_1.rds")
    }
    mydata <- mydata[eids, ]
    
    # Survival model (PCE only)
    cindex_stable <- RunCoxSelected(mydata = mydata, mydata_test = mydata_test, selected = "logHR", penalisation = penalised, digits = 5)
    saveRDS(cindex_stable, paste0("Results/Cox_models/cox_pce_", gender, "_", outcome, ".rds"))
    
    # Survival model (stably selected)
    for (model_id in c(1:4, 7)) {
      print(model_id)
      
      # Extracting stably selected variables
      stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
      class(stab) <- "variable_selection"
      selected <- SelectedVariables(stab)
      print(sum(selected))
      selected <- names(selected)[selected == 1]
      
      # Refitting the Cox model
      cindex_stable <- RunCoxSelected(mydata = mydata, mydata_test = mydata_test, selected = selected, penalisation = penalised, digits = 5)
      saveRDS(cindex_stable, paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
    }
    
    # Logistic model (PCE only)
    refitted=Refit(xdata=mydata[,"logHR",drop=FALSE], ydata = mydata$case, family = "binomial")
    predicted=predict.glm(refitted, newdata = mydata_test)
    myroc=pROC::roc(predictor=predicted, response=mydata_test$case, ci=TRUE)
    logit_pce=list(glm=refitted,
                   predicted_test=predicted, 
                   observed_test=mydata_test$case,
                   c_test=ReformatAUC(as.numeric(myroc$ci)))
    saveRDS(logit_pce, paste0("Results/Cox_models/logit_pce_", gender, "_", outcome, ".rds"))
    
    # Logistic model (stably selected)
    for (model_id in c(1:4, 7)) {
      print(model_id)
      
      # Extracting stably selected variables
      stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
      class(stab) <- "variable_selection"
      selected <- SelectedVariables(stab)
      print(sum(selected))
      selected <- names(selected)[selected == 1]
      
      # Refitting the logistic model
      refitted=Refit(xdata=mydata[,selected,drop=FALSE], ydata = mydata$case, family = "binomial")
      predicted=predict.glm(refitted, newdata = mydata_test)
      myroc=pROC::roc(predictor=predicted, response=mydata_test$case, ci=TRUE)
      logit_stable=list(glm=refitted,
                        predicted_test=predicted, 
                        observed_test=mydata_test$case,
                        c_test=ReformatAUC(as.numeric(myroc$ci)))
      saveRDS(logit_stable, paste0("Results/Cox_models/logit_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
    }
    cat("\n")
  }
}
