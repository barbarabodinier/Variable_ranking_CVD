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

for (model_id in model_list) {
  print(model_id)

  mytable <- c("N", "PCE", "LASSO")
  for (gender in genders) {
    # Defining input/output name
    full_model_name <- paste0(outcome, "_m", model_id, "_", data_input, "_", gender)

    # Loading refitted model outputs
    cindex_pce <- readRDS(paste0(
      "Results/Refitted_models/cox_pce_",
      full_model_name, ".rds"
    ))
    cindex_lasso <- readRDS(paste0(
      "Results/Refitted_models/cox_lasso_",
      full_model_name, ".rds"
    ))

    tmp <- c(
      paste0(
        formatC(sum(cindex_pce$observed_test == 1), big.mark = ","),
        "/",
        formatC(sum(cindex_pce$observed_test == 0), big.mark = ",")
      ),
      cindex_pce$c_test,
      cindex_lasso$c_test
    )
    mytable <- cbind(mytable, tmp)
  }
  colnames(mytable) <- c("", genders)
  write.xlsx(
    data.frame(mytable),
    paste0("Tables/C_statistics_", outcome, "_m", model_id, "_", data_input, ".xlsx")
  )
}
