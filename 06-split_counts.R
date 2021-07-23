rm(list = ls())

library(focus)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)
library(plotrix)
library(colorspace)
source("Scripts/functions.R")

dir.create("Results/Cox_models", showWarnings = FALSE)

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
mydata <- cbind(scale(mydata[, rownames(mydict)]), mydata[, c("case", "time")])

print(nrow(mydata))
print(nrow(na.exclude(mydata)))

for (gender in c("male", "female")) {
  print(gender)
  # Extracting selection set
  if (gender == "female") {
    eids <- readRDS("Data/Split/selection_set_eids_0_0.rds")
  } else {
    eids <- readRDS("Data/Split/selection_set_eids_1_0.rds")
  }
  mydata_test <- mydata[eids, ]
  print("Selection:")
  print(table(mydata_test$case))
  print(table(na.exclude(mydata_test)$case))
  cat("\n")

  # Extracting training set
  if (gender == "female") {
    eids <- readRDS("Data/Split/performance_set_eids_0_0.rds")
  } else {
    eids <- readRDS("Data/Split/performance_set_eids_1_0.rds")
  }
  mydata_test <- mydata[eids, ]
  print("Training:")
  print(table(mydata_test$case))
  print(table(na.exclude(mydata_test)$case))
  cat("\n")

  # Extracting selection set
  if (gender == "female") {
    eids <- readRDS("Data/Split/estimation_set_eids_0_0.rds")
  } else {
    eids <- readRDS("Data/Split/estimation_set_eids_1_0.rds")
  }
  mydata_test <- mydata[eids, ]
  print("Test:")
  print(table(mydata_test$case))
  print(table(na.exclude(mydata_test)$case))
  cat("\n")
}
