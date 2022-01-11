rm(list = ls())

library(openxlsx)
library(focus)

# Parameters
data_input <- "updated"
outcome <- "cvd"

# Variable annotation
mydict <- read.xlsx("Data/Data_dictionary.xlsx")

# Extracting ever-used variables
mydict <- mydict[which(apply(mydict[, 4:ncol(mydict)], 1, FUN = function(x) {
  any(x == "X")
})), ]
rownames(mydict) <- mydict[, 1]

# Loading the data
mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_updated.rds")))
mydata <- cbind(scale(mydata[, rownames(mydict)]), mydata[, c("case", "time", "sex", "logHR")])

for (model_id in c(1, 3)) {
  for (gender in c("Male", "Female")) {

    # Selection proportions
    stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id, "_", data_input, "_", tolower(gender), ".rds"))
    class(stab) <- "variable_selection"

    # Preparing sex-specific data
    if (gender == "male") {
      tmpdata <- mydata[which(mydata$sex == 1), ]
    } else {
      tmpdata <- mydata[which(mydata$sex == 0), ]
    }

    # Identifying stably selected predictors
    predictors <- names(SelectedVariables(stab))[which(SelectedVariables(stab) == 1)]

    # Defining order of inclusion for incremental analyses
    selprop <- SelectionProportions(stab)
    myorder <- data.frame(selprop = selprop, cat = mydict[names(selprop), 3])
    myorder$cat <- factor(myorder$cat, levels = c("Established Risk Factors", "Genetic", "Biochemistry", "Haematology", "Nightingale"))
    myorder <- myorder[with(myorder, order(selprop, -as.numeric(cat), decreasing = TRUE)), ]

    # Extracting stably selected predictors only
    tmporder <- myorder[which(rownames(myorder) %in% predictors), ]

    # Data preparation for performance of calibrated model
    x <- tmpdata[, rownames(tmporder)]
    x <- na.exclude(x)
    tmpdata <- tmpdata[rownames(x), ]
    y <- tmpdata[, c("time", "case")]

    cstat <- ExplanatoryPerformance(
      xdata = x, ydata = y, tau = 0.5,
      family = "cox", ij_method = TRUE
    )
    saveRDS(cstat, paste0("Outputs/expl_perf_", outcome, "_m", model_id, "_", tolower(gender), ".rds"))

    # Calculating C statistic for PCE only (on same data)
    x <- tmpdata[rownames(x), "logHR", drop = FALSE]
    y <- tmpdata[rownames(x), c("time", "case")]
    cstat_pce <- ExplanatoryPerformance(
      xdata = x, ydata = y, tau = 0.5,
      family = "cox", ij_method = TRUE
    )
    saveRDS(cstat_pce, paste0("Outputs/expl_perf_pce_", outcome, "_m", model_id, "_", tolower(gender), ".rds"))

    # Data preparation for incremental analyses
    x <- tmpdata[, rownames(myorder)]
    x <- na.exclude(x)
    tmpdata <- tmpdata[rownames(x), ]
    y <- tmpdata[, c("time", "case")]
    incremental <- Incremental(
      xdata = x, ydata = y, tau = 0.5,
      family = "cox", ij_method = TRUE
    )
    saveRDS(incremental, paste0("Outputs/incr_perf_", outcome, "_m", model_id, "_", tolower(gender), ".rds"))
  }
}
