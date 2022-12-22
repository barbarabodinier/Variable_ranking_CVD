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
model_id <- 1

for (gender in genders) {
  print(gender)

  # Defining input/output name
  full_model_name <- paste0(outcome, "_m", model_id, "_", data_input, "_", gender)

  # Loading refitted model output
  cindex_lasso <- readRDS(paste0(
    "Results/Refitted_models/cox_lasso_",
    full_model_name, ".rds"
  ))

  # Showing model coefficients
  print(cindex_lasso$cox_model)
  cat("\n")
}
