rm(list = ls())

library(focus)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)
library(plotrix)
library(colorspace)
source("Scripts/functions.R")

genders <- c("male", "female")
outcome_names <- c("cvd")
mymethod <- "lasso"
data_input <- "updated"

# for (model_id in c(0:4,7)) {
for (model_id in c(1,3)) {
  for (i in 1:length(outcome_names)) {
    outcome <- outcome_names[i]
    print(outcome)
    for (j in 1:length(genders)) {
      gender <- genders[j]
      print(gender)

      full_data <- TRUE

      # Loading the data
      mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_", data_input, ".rds")))

      # Variable annotation
      mydict <- read.xlsx("Data/Data_dictionary.xlsx")

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
      mydata_test <- mydata[eids, ]

      # Extracting training set
      if (gender == "female") {
        eids <- readRDS("Data/Split/estimation_set_eids_0_1.rds")
      } else {
        eids <- readRDS("Data/Split/estimation_set_eids_1_1.rds")
      }
      mydata <- mydata[eids, ]

      print(dim(mydata))
      print(dim(mydata_test))
      if (model_id %in% c(0, 3, 4, 7)) {
        mydata <- na.exclude(mydata)
        mydata_test <- na.exclude(mydata_test)
        if (model_id == 0) {
          full_data <- FALSE
          model_id <- 1
        }
      }
      print(dim(mydata))
      print(dim(mydata_test))

      # Reading LASSO results
      stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
      selprop <- SelectionProportions(stab)

      if (model_id %in% c(1, 3)) {
        # Reading Random Forest results
        if (gender == "male") { # men and women inverted
          rf <- readRDS(paste0("Results/Random_forest/mod", model_id, "men_VarImp_results_pvals.rds"))
        } else {
          rf <- readRDS(paste0("Results/Random_forest/mod", model_id, "women_VarImp_results_pvals.rds"))
        }
        rownames(rf) <- rf$var
        print(all(rownames(rf) %in% names(selprop)))
        print(all(names(selprop) %in% rownames(rf)))
        rf <- rf[names(selprop), ]

        myorder <- data.frame(selprop = selprop, cat = mydict[names(selprop), 3], varimp = rf$mean)
      } else {
        myorder <- data.frame(selprop = selprop, cat = mydict[names(selprop), 3])
      }
      myorder$cat <- factor(myorder$cat, levels = c("Established Risk Factors", "Genetic", "Biochemistry", "Haematology", "Nightingale"))

      if (model_id %in% c(1, 3)) {
        # Saving selection counts
        mytable <- myorder
        mytable <- mytable[with(myorder, order(-as.numeric(cat), selprop, varimp, decreasing = TRUE)), ]
        mytable[, 1] <- mytable[, 1] * 1000
        mytable[, 3] <- paste0(
          formatC(rf[rownames(mytable), "mean"], format = "e", digits = 2), " (",
          formatC(rf[rownames(mytable), "lower_CI"], format = "e", digits = 2), ", ",
          formatC(rf[rownames(mytable), "upper_CI"], format = "e", digits = 2), ")"
        )
        rownames(mytable) <- mydict[rownames(mytable), 2]
        mytable <- mytable[, c(1, 3, 2)]
        write.xlsx(mytable, paste0("Tables/Selection_counts_", outcome, "_m", model_id, "_", data_input, "_", gender, ".xlsx"), 
                   row.names = TRUE, overwrite = TRUE)
      }

      # Defining order of inclusion
      myorder$cat <- -as.numeric(myorder$cat)
      if ("logHR" %in% rownames(myorder)) {
        myorder["logHR", 2] <- 0
      }
      myorder <- myorder[with(myorder, order(selprop, cat, decreasing = TRUE)), ]
      myorder <- rownames(myorder)
      names(myorder) <- mydict[myorder, 2]

      if (mymethod == "_rf") {
        myorder <- -log10(readRDS(paste0("Results/Random_Forest_variable_importance_non_corrected_pvalue_", gender, outcome, ".rds")))
      }

      # Keep those selected at least once / nominally significant variable importance
      variable_cat <- variable_cat[names(variable_cat)[names(variable_cat) %in% names(myorder)]]
      if (mymethod == "lasso") {
        myorder <- myorder[myorder > 0]
        all_variables <- names(myorder)[sort.list(myorder, decreasing = TRUE)]
      }
      if (mymethod == "_rf") {
        vi <- readRDS(paste0("Results/Random_Forest_variable_importance_", gender, outcome, ".rds"))
        print(all(names(myorder) == names(vi)))
        myorder <- myorder[order(-myorder, -vi)]
        all_variables <- names(myorder)[myorder > -log10(0.3)]
        pvalue <- readRDS(paste0("Results/Random_Forest_variable_importance_pvalue_", gender, outcome, ".rds"))
        print(sum(myorder >= min(myorder[names(pvalue)[pvalue < 0.05]])) == sum(pvalue < 0.05))
        print(sum(names(which(myorder >= min(myorder[names(pvalue)[pvalue < 0.05]]))) %in% names(which(pvalue < 0.05))) == sum(pvalue < 0.05))
      }
      all_variables <- unname(myorder)
      mycat <- variable_cat

      # Run Cox models with variables added one by one
      penalised <- FALSE
      c_training <- c_test <- NULL
      for (k in 1:length(all_variables)) {
        print(k)
        selected <- all_variables[1:k]
        cindex_tmp <- RunCoxSelected(mydata = mydata, mydata_test = mydata_test, selected = selected, penalisation = penalised, digits = 5)
        c_training <- c(c_training, cindex_tmp$c_training)
        c_test <- c(c_test, cindex_tmp$c_test)
      }
      names(c_training) <- all_variables
      names(c_test) <- all_variables

      if (!full_data) {
        model_id <- 0
      }

      saveRDS(c_training, paste0("Results/C_statistic_incremental_training_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
      saveRDS(c_test, paste0("Results/C_statistic_incremental_test_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))

      # Run Cox model with PCE as predictor
      penalised <- FALSE
      cindex_tmp <- RunCoxSelected(mydata = mydata, mydata_test = mydata_test, selected = "logHR", penalisation = penalised, digits = 5)
      c_training <- cindex_tmp$c_training
      c_test <- cindex_tmp$c_test
      saveRDS(c_training, paste0("Results/C_statistic_incremental_training_", outcome, "_m", model_id, "_PCE_", data_input, "_", gender, ".rds"))
      saveRDS(c_test, paste0("Results/C_statistic_incremental_test_", outcome, "_m", model_id, "_PCE_", data_input, "_", gender, ".rds"))
    }
  }
}
