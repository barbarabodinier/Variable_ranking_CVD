rm(list = ls())

library(sharp)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)

source("Scripts/functions.R")

dir.create("Results/Incremental", showWarnings = FALSE)

# Parameters
outcome <- "cvd"
data_input <- "updated"
genders <- c("male", "female")

for (model_id in c(1, 3, 7)) {
  print(model_id)
  for (gender in genders) {
    print(gender)

    full_data <- TRUE

    # Variable annotation
    mydict <- read.xlsx("Data/Data_dictionary.xlsx")

    # Loading the data
    mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_", data_input, ".rds")))

    # Using subset with Nightingale data
    if (model_id %in% c(3, 4, 7)) {
      id_column <- which(tolower(colnames(mydict)) == paste0("model.", 3, ".(", gender, ")"))
      predictors <- mydict[which(mydict[, id_column] == "X"), 1]
      z <- mydata[, predictors]
      z <- na.exclude(z)
      mydata <- mydata[rownames(z), ]
    }

    # Extracting ever-used variables
    mydict <- mydict[which(apply(mydict[, 4:ncol(mydict)], 1, FUN = function(x) {
      any(x == "X")
    })), ]
    rownames(mydict) <- mydict[, 1]
    mydata <- mydata[, c(rownames(mydict), "time", "case")]
    variable_cat <- mydict[, 3]
    names(variable_cat) <- mydict[, 2]

    # Extracting test set
    if (gender == "female") {
      eids <- readRDS("Data/Split/performance_set_eids_0_1.rds")
    } else {
      eids <- readRDS("Data/Split/performance_set_eids_1_1.rds")
    }
    eids <- intersect(eids, rownames(mydata))
    mydata_test <- mydata[eids, ]

    # Extracting training set
    if (gender == "female") {
      eids <- readRDS("Data/Split/estimation_set_eids_0_1.rds")
    } else {
      eids <- readRDS("Data/Split/estimation_set_eids_1_1.rds")
    }
    eids <- intersect(eids, rownames(mydata))
    mydata <- mydata[eids, ]

    # Printing dimensions
    print(dim(mydata))
    print(dim(mydata_test))

    # Reading LASSO results
    stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
    class(stab) <- "variable_selection"
    selprop <- SelectionProportions(stab)
    myorder <- data.frame(selprop = selprop, cat = mydict[names(selprop), 3])
    myorder$cat <- factor(myorder$cat, levels = c("Established Risk Factors", "Genetic", "Biochemistry", "Haematology", "Nightingale"))

    # Defining order of inclusion
    myorder$cat <- -as.numeric(myorder$cat)
    if ("logHR" %in% rownames(myorder)) {
      myorder["logHR", 2] <- 0
    }
    myorder <- myorder[with(myorder, order(selprop, cat, decreasing = TRUE)), ]
    myorder <- rownames(myorder)
    names(myorder) <- mydict[myorder, 2]

    # Keep those selected at least once / nominally significant variable importance
    variable_cat <- variable_cat[names(variable_cat)[names(variable_cat) %in% names(myorder)]]
    myorder <- myorder[myorder > 0]
    all_variables <- names(myorder)[sort.list(myorder, decreasing = TRUE)]
    all_variables <- unname(myorder)
    mycat <- variable_cat

    # Run Cox models with variables added one by one
    c_test <- NULL
    for (k in 1:length(all_variables)) {
      print(k)
      selected <- all_variables[1:k]
      cindex_tmp <- RunCoxSelected(
        x = mydata,
        x_test = mydata_test,
        selected = selected,
        verbose = FALSE,
        digits = 5
      )
      c_test <- c(c_test, cindex_tmp$c_test)
    }
    names(c_test) <- all_variables

    if (!full_data) {
      model_id <- 0
    }
    saveRDS(c_test, paste0("Results/Incremental/C_statistic_incremental_test_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))

    # Run Cox model with PCE as predictor
    cindex_tmp <- RunCoxSelected(
      x = mydata,
      x_test = mydata_test,
      selected = "logHR",
      verbose = FALSE,
      digits = 5
    )
    c_test <- cindex_tmp$c_test
    saveRDS(c_test, paste0("Results/Incremental/C_statistic_incremental_test_", outcome, "_m", model_id, "_PCE_", data_input, "_", gender, ".rds"))
  }
}
