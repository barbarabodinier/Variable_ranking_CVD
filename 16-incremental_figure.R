rm(list = ls())

library(sharp)
library(survival)
library(openxlsx)
library(plotrix)
library(colorspace)

source("Scripts/functions.R")

# Parameters
outcome <- "cvd"
genders <- c("male", "female")
mymethod <- "lasso"
data_input <- "updated"
cex.axis <- 0.8

for (model_id in c(1, 3, 7)) {
  i <- 0

  # Incremental figure
  background <- TRUE
  {
    pdf(paste0("Figures/incremented_c_statistic_", outcome, "_m", model_id, "_", data_input, ".pdf"), useDingbats = FALSE, width = 14, height = 16)
    par(mar = c(18, 5.5, 1, 1), mfrow = c(2, 1))

    for (j in 1:length(genders)) {
      gender <- genders[j]
      print(gender)

      # Reading the output
      c_test_lasso <- readRDS(paste0("Results/Incremental/C_statistic_incremental_test_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
      if (model_id == 0) {
        stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m1_", data_input, "_", gender, ".rds"))
      } else {
        stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
      }
      class(stab) <- "variable_selection"
      c_test_rf <- c_test_lasso

      # Performance of PCE only
      c_test_pce <- readRDS(paste0("Results/Incremental/C_statistic_incremental_test_", outcome, "_m", model_id, "_PCE_", data_input, "_", gender, ".rds"))
      i <- i + 1
      c_test <- eval(parse(text = paste0("c_test_", mymethod)))
      myorder <- names(c_test)

      # Variable annotation
      data_dico <- read.xlsx("Data/Data_dictionary.xlsx")
      variable_names <- data_dico[, 2]
      names(variable_names) <- data_dico[, 1]
      variable_cat <- data_dico[, 3]
      names(variable_cat) <- data_dico[, 2]

      all_variables <- myorder
      all_variables <- variable_names[all_variables]
      mycat <- variable_cat

      # Definition of the colours
      mycolours <- c(
        lighten("darkolivegreen", amount = 0.1),
        darken("royalblue", amount = 0.2),
        darken("chocolate", amount = 0.05),
        lighten("red", amount = 0.1),
        darken("darkorchid4", amount = 0.05),
        darken("gold", amount = 0.2)
      )
      names(mycolours) <- unique(mycat)

      # Reformat units
      mynames <- all_variables
      mynames <- gsub("109", "'*10^9*'", mynames)
      mynames <- gsub("1012", "'*10^12*'", mynames)
      mynames <- gsub("m2", "'*m^2*'", mynames)
      mynames <- gsub("micro", "'*mu*'", mynames)

      # English/American spelling
      id <- grep("Haemoglobin", mynames)
      mynames[id] <- "Hemoglobin"
      id <- grep("Glycated haemoglobin", mynames)
      mynames[id] <- "Glycated hemoglobin"
      id <- grep("Haematocrit", mynames)
      mynames[id] <- "Hematocrit"
      id <- grep("Low density cholesterol", mynames)
      mynames[id] <- "Low-density lipoprotein cholesterol"
      id <- grep("High density lipoprotein cholesterol", mynames)
      mynames[id] <- "High-density lipoprotein cholesterol"

      yrange <- range(c(
        as.numeric(gsub("-.*", "", gsub(".* \\[", "", c_test_lasso))),
        as.numeric(gsub("-.*", "", gsub(".* \\[", "", c_test_lasso))),
        as.numeric(gsub("\\]", "", gsub(".*-", "", gsub(".* \\[", "", c_test_rf)))),
        as.numeric(gsub("\\]", "", gsub(".*-", "", gsub(".* \\[", "", c_test_rf))))
      )) + c(0, 0.005)
      plotCI(
        x = 1:length(c_test), y = as.numeric(gsub(" .*", "", c_test)),
        li = as.numeric(gsub("-.*", "", gsub(".* \\[", "", c_test))),
        ui = as.numeric(gsub("\\]", "", gsub(".*-", "", gsub(".* \\[", "", c_test)))),
        pch = 18, cex = 1.5, las = 1, xlab = "", ylab = "C statistic", cex.lab = 1.5, xaxt = "n",
        col = NA, ylim = yrange
      )

      par(xpd = TRUE)
      if (model_id == 3) {
        mtext(side = 2, text = LETTERS[i + 2], at = max(yrange), las = 1, line = 3, cex = 3)
      } else {
        mtext(side = 2, text = LETTERS[i], at = max(yrange), las = 1, line = 3, cex = 3)
      }
      par(xpd = FALSE)

      xseqgreysep <- c(1:(length(c_test) + 1)) - 0.5
      if (background) {
        ymax <- max(as.numeric(gsub("\\]", "", gsub(".*-", "", gsub(".* \\[", "", c_test)))), na.rm = TRUE) + 10
        for (k in seq(1, length(xseqgreysep), by = 2)) {
          polygon(
            x = c(xseqgreysep[k], xseqgreysep[k + 1], xseqgreysep[k + 1], xseqgreysep[k]),
            y = c(-ymax, -ymax, ymax, ymax), col = lighten("black", amount = 0.97), border = NA
          )
        }
        for (k in seq(2, length(xseqgreysep), by = 2)) {
          polygon(
            x = c(xseqgreysep[k], xseqgreysep[k + 1], xseqgreysep[k + 1], xseqgreysep[k]),
            y = c(-ymax, -ymax, ymax, ymax), col = lighten("black", amount = 0.99), border = NA
          )
        }
      }
      abline(h = axTicks(2), lty = 3, col = "grey")
      abline(v = xseqgreysep, lty = 3, col = "grey95")
      myc_pce <- as.numeric(gsub(" \\[.*", "", c_test_pce))
      abline(h = myc_pce, lty = 2, col = "black")
      text(x = max(xseqgreysep), y = myc_pce, labels = "PCE", pos = 1)
      box()
      if (mymethod == "lasso") {
        abline(v = sum(SelectedVariables(stab)) + 0.5, lty = 2, col = "red")
      }
      plotCI(
        x = 1:length(c_test), y = as.numeric(gsub(" .*", "", c_test)),
        li = as.numeric(gsub("-.*", "", gsub(".* \\[", "", c_test))),
        ui = as.numeric(gsub("\\]", "", gsub(".*-", "", gsub(".* \\[", "", c_test)))),
        pch = 18, cex = 1.2, las = 1, xlab = "", ylab = "C statistic", cex.lab = 1.5, xaxt = "n",
        col = darken(mycolours[mycat[all_variables]], amount = 0.3), sfrac = 0.003, add = TRUE
      )
      for (k in 1:length(c_test)) {
        if (k == 1) {
          axis(
            side = 1, at = k, labels = eval(parse(text = paste0("expression('", mynames[k], "')"))),
            col.axis = darken(mycolours[mycat[(all_variables[k])]], amount = 0.4), las = 2, cex.axis = cex.axis
          )
        } else {
          axis(
            side = 1, at = k, labels = eval(parse(text = paste0("expression('+ ", mynames[k], "')"))),
            col.axis = darken(mycolours[mycat[(all_variables[k])]], amount = 0.4), las = 2, cex.axis = cex.axis
          )
        }
      }
    }
    dev.off()
  }
}
