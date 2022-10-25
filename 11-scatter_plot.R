rm(list = ls())

library(sharp)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)
library(plotrix)
library(colorspace)
library(basicPlotteR)
source("Scripts/functions.R")

dir.create("Tables", showWarnings = FALSE)

outcome_names <- c("cvd")
model_id <- c(3, 7)
xlab <- "Base model"
ylab <- "Base model + Nightingale"
# model_id <- c(2, 1)
# xlab="Base model"
# ylab="PCE in place of constituent variables"

for (data_input in c("updated")) {
  for (i in 1:length(outcome_names)) {
    outcome <- outcome_names[i]
    print(outcome)

    # Variable annotation
    dict <- read.xlsx("Data/Data_dictionary.xlsx")
    rownames(dict) <- dict[, 1]
    variable_cat <- dict[, 3]
    names(variable_cat) <- dict[, 2]
    nightingale <- names(which(variable_cat == "Nightingale"))

    # Extracting selection proportions from both models
    for (j in 1:length(model_id)) {
      # Selection proportions (Men)
      stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id[j], "_", data_input, "_male.rds"))
      class(stab) <- "variable_selection"
      selprop_male <- SelectionProportions(stab)
      names(selprop_male) <- dict[names(selprop_male), 2]
      hat_pi_male <- Argmax(stab)[2]

      # Selection proportions (Women)
      stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id[j], "_", data_input, "_female.rds"))
      class(stab) <- "variable_selection"
      selprop_female <- SelectionProportions(stab)
      selprop_female <- c(selprop_female, ED = 0)
      names(selprop_female) <- dict[names(selprop_female), 2]
      hat_pi_female <- Argmax(stab)[2]

      if (model_id[j] == 3) {
        selprop_nightingale_male <- selprop_male[which(names(selprop_male) %in% nightingale)]
        selprop_nightingale_female <- selprop_female[which(names(selprop_female) %in% nightingale)]
      } else {
        if (!3 %in% model_id) {
          selprop_nightingale_male <- selprop_nightingale_female <- 2
        }
      }

      assign(paste0("selprop_male_", j), selprop_male)
      assign(paste0("selprop_female_", j), selprop_female)
      assign(paste0("hat_pi_male_", j), hat_pi_male)
      assign(paste0("hat_pi_female_", j), hat_pi_female)
    }

    variable_cat <- variable_cat[names(variable_cat)[names(variable_cat) %in% intersect(names(selprop_male_1), names(selprop_male_2))]]
    selprop_male <- selprop_male[names(variable_cat)]
    selprop_female <- selprop_female[names(variable_cat)]
    selprop_female <- selprop_female[names(selprop_male)]
    print(all(names(selprop_male) == names(selprop_female)))
    mycat <- variable_cat

    for (j in 1:length(model_id)) {
      selprop_male <- eval(parse(text = paste0("selprop_male_", j)))
      selprop_female <- eval(parse(text = paste0("selprop_female_", j)))

      selprop_male <- selprop_male[names(mycat)]
      selprop_female <- selprop_female[names(mycat)]

      assign(paste0("selprop_male_", j), selprop_male)
      assign(paste0("selprop_female_", j), selprop_female)
    }

    # Change group name
    if (outcome == "cvd") {
      mycat[mycat == "Established Risk Factors"] <- "PCE/QRISK3"
    } else {
      mycat[mycat == "Established Risk Factors"] <- "PCE/QRISK3"
    }

    # Define colours by category
    mycolours <- c(
      lighten("darkolivegreen", amount = 0.1),
      darken("royalblue", amount = 0.2),
      darken("chocolate", amount = 0.05),
      lighten("red", amount = 0.1),
      darken("darkorchid4", amount = 0.05),
      darken("gold", amount = 0.2)
    )
    names(mycolours) <- c("Biochemistry", "PCE/QRISK3", "Genetic", "Haematology", "Nightingale", "Diet")
    mycolours <- mycolours[mycat]

    # Reformat units
    names(selprop_male) <- gsub("109", "'*10^9*'", names(selprop_male))
    names(selprop_male) <- gsub("1012", "'*10^12*'", names(selprop_male))
    names(selprop_male) <- gsub("m2", "'*m^2*'", names(selprop_male))
    names(selprop_male) <- gsub("micro", "'*mu*'", names(selprop_male))

    # English/American spelling
    id <- grep("Haemoglobin", names(selprop_male))
    names(selprop_male)[id] <- "Hemoglobin"
    names(selprop_female)[id] <- "Hemoglobin"
    id <- grep("Haematocrit", names(selprop_male))
    names(selprop_male)[id] <- "Hematocrit"
    names(selprop_female)[id] <- "Hematocrit"

    plotname <- paste0("Figures/Scatter_m", model_id[1], "_m", model_id[2], "_", outcome, "_", data_input, ".pdf")
    {
      pdf(plotname, width = 14, height = 7)
      par(mfrow = c(1, 2), mar = c(5, 5, 1, 1))
      # Male
      plot(NA,
        xlab = xlab, ylab = ylab,
        xlim = c(0, 1), ylim = c(0, 1),
        cex.lab = 1.5,
        xaxt = "n", yaxt = "n"
      )
      abline(h = hat_pi_male_2, lty = 2, col = "darkred")
      abline(v = hat_pi_male_1, lty = 2, col = "darkred")
      abline(0, 1, lty = 2)
      for (l in 1:length(selprop_nightingale_male)) {
        abline(h = selprop_nightingale_male[l], col = darken("darkorchid4", amount = 0.05), lty = 3)
      }
      text(
        x = 0, y = selprop_nightingale_male,
        pos = 3, offset = 0.3, cex = 1.2,
        col = darken("darkorchid4", amount = 0.4),
        labels = ifelse(selprop_nightingale_male >= 0.3,
          yes = seq(
            length(selprop_male_2) + 1,
            length(selprop_male_2) + length(selprop_nightingale_male)
          ),
          no = ""
        )
      )
      points(selprop_male_2, selprop_male_1,
        pch = 19, col = mycolours
      )
      axis(side = 1, at = seq(0, 1, by = 0.1))
      axis(side = 2, at = seq(0, 1, by = 0.1), las = 1)
      mtext(text = "A", side = 2, at = 1, las = 1, cex = 3, line = 2.7)
      ids <- which(selprop_male_1 >= 0.07)
      addTextLabels(selprop_male_2[ids],
        selprop_male_1[ids],
        keepLabelsInside = FALSE,
        col.label = darken(mycolours, amount = 0.4)[ids],
        col.line = mycolours[ids],
        col.background = lighten(mycolours[ids], amount = 0.8),
        border = lighten(mycolours[ids], amount = 0.7),
        cex.label = 1.2,
        labels = seq(1, length(selprop_male_2))[ids]
      )

      # Female
      plot(NA,
        xlab = xlab, ylab = ylab,
        xlim = c(0, 1), ylim = c(0, 1),
        cex.lab = 1.5,
        xaxt = "n", yaxt = "n"
      )
      abline(h = hat_pi_female_2, lty = 2, col = "darkred")
      abline(v = hat_pi_female_1, lty = 2, col = "darkred")
      abline(0, 1, lty = 2)
      for (l in 1:length(selprop_nightingale_female)) {
        abline(h = selprop_nightingale_female[l], col = darken("darkorchid4", amount = 0.05), lty = 3)
      }
      text(
        x = 0, y = selprop_nightingale_female,
        pos = 3, offset = 0.3, cex = 1.2,
        col = darken("darkorchid4", amount = 0.4),
        labels = ifelse(selprop_nightingale_female >= 0.3,
          yes = seq(
            length(selprop_female_2) + 1,
            length(selprop_female_2) + length(selprop_nightingale_female)
          ),
          no = ""
        )
      )
      points(selprop_female_2, selprop_female_1,
        pch = 19, col = mycolours
      )
      axis(side = 1, at = seq(0, 1, by = 0.1))
      axis(side = 2, at = seq(0, 1, by = 0.1), las = 1)
      mtext(text = "B", side = 2, at = 1, las = 1, cex = 3, line = 2.7)
      ids <- which(selprop_female_1 >= 0.07)
      addTextLabels(selprop_female_2[ids],
        selprop_female_1[ids],
        keepLabelsInside = FALSE,
        col.label = darken(mycolours, amount = 0.4)[ids],
        col.line = mycolours[ids],
        col.background = lighten(mycolours[ids], amount = 0.8),
        border = lighten(mycolours[ids], amount = 0.7),
        cex.label = 1.2,
        labels = seq(1, length(selprop_female_2))[ids]
      )
      dev.off()
    }


    plotname <- paste0("Figures/Scatter_legend_m", model_id[1], "_m", model_id[2], "_", outcome, "_", data_input, ".pdf")
    {
      pdf(plotname, width = 20, height = 7)
      par(mfrow = c(1, 1), xpd = NA)
      plot.new()
      legend("top",
        bty = "n", ncol = 4,
        text.col = darken(c(mycolours, rep(darken("darkorchid4", amount = 0.05), length(selprop_nightingale_female))), amount = 0.4),
        legend = paste0(
          seq(1, length(selprop_female_2) + length(selprop_nightingale_female)),
          " - ",
          c(names(selprop_female_2), names(selprop_nightingale_female))
        )
      )
      dev.off()
    }
    system(paste("pdfcrop --margin 10", plotname, plotname))
  }
}
