rm(list = ls())

library(sharp)
library(survival)
library(openxlsx)
library(colorspace)

# source("Scripts/functions.R")

dir.create("Tables", showWarnings = FALSE)

# Parameters
outcome_names <- "cvd"
data_input <- "updated"
model_list <- c(1:4, 7)

# Looping over models
for (i in 1:length(outcome_names)) {
  outcome <- outcome_names[i]
  print(outcome)
  
  for (model_id in model_list) {
    print(model_id)
    
    # Variable annotation
    dict <- read.xlsx("Data/Data_dictionary.xlsx")
    rownames(dict) <- dict[, 1]
    variable_cat <- dict[, 3]
    names(variable_cat) <- dict[, 2]
    
    # Selection proportions (Men)
    stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id, "_", data_input, "_male.rds"))
    class(stab) <- "variable_selection"
    {
      pdf(paste0("Figures/Calibration_m", model_id, "_", outcome, "_", data_input, "_male.pdf"),
          width = 12, height = 7, useDingbats = FALSE
      )
      par(mar = c(7, 5, 7, 6))
      CalibrationPlot(stab)
      dev.off()
    }
    selprop_male <- SelectionProportions(stab)
    names(selprop_male) <- dict[names(selprop_male), 2]
    hat_pi1 <- Argmax(stab)[2]
    
    # Selection proportions (Women)
    stab <- readRDS(paste0("Results/HPC_results/stability_", outcome, "_m", model_id, "_", data_input, "_female.rds"))
    class(stab) <- "variable_selection"
    {
      pdf(paste0("Figures/Calibration_m", model_id, "_", outcome, "_", data_input, "_female.pdf"),
          width = 12, height = 7, useDingbats = FALSE
      )
      par(mar = c(7, 5, 7, 6))
      CalibrationPlot(stab)
      dev.off()
    }
    selprop_female <- SelectionProportions(stab)
    selprop_female <- c(selprop_female, ED = 0)
    names(selprop_female) <- dict[names(selprop_female), 2]
    hat_pi2 <- Argmax(stab)[2]
    
    variable_cat <- variable_cat[names(variable_cat)[names(variable_cat) %in% names(selprop_male)]]
    selprop_male <- selprop_male[names(variable_cat)]
    selprop_female <- selprop_female[names(variable_cat)]
    selprop_female <- selprop_female[names(selprop_male)]
    print(all(names(selprop_male) == names(selprop_female)))
    mycat <- variable_cat
    
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
    
    # Figure
    background <- TRUE
    plotname <- paste0("Figures/Selected_variables_m", model_id, "_", outcome, "_", data_input, ".pdf")
    {
      pdf(plotname, width = 12, height = 10.5, useDingbats = FALSE)
      mylwd <- 3
      par(mar = c(0.1, 7, 0, 1), mfrow = c(4, 1), lheight = 0.2)
      
      plot.new()
      
      plot(1:length(selprop_female), selprop_female,
           ylim = c(0, 1),
           col = ifelse(selprop_female >= hat_pi2, yes = darken(mycolours, amount = 0.3), no = lighten(mycolours, amount = 0.1)),
           type = "h", lwd = mylwd, xlab = "", ylab = expression(atop(NA, atop(
             textstyle("Selection Proportion"),
             textstyle("(Women)")
           ))),
           xaxt = "n", cex.lab = 1.3, panel.first = abline(v = 1:length(selprop_female), lty = 3, col = "grey")
      )
      xseqgreysep <- c(1:(length(selprop_female) + 1)) - 0.5
      if (background) {
        ymax <- max(selprop_female, na.rm = TRUE) + 10
        for (k in seq(1, length(xseqgreysep), by = 2)) {
          polygon(
            x = c(xseqgreysep[k], xseqgreysep[k + 1], xseqgreysep[k + 1], xseqgreysep[k]),
            y = c(-ymax, -ymax, ymax, ymax), col = lighten(mycolours[mycat[k]], amount = 0.97), border = NA
          )
        }
        for (k in seq(2, length(xseqgreysep), by = 2)) {
          polygon(
            x = c(xseqgreysep[k], xseqgreysep[k + 1], xseqgreysep[k + 1], xseqgreysep[k]),
            y = c(-ymax, -ymax, ymax, ymax), col = lighten(mycolours[mycat[k]], amount = 0.99), border = NA
          )
        }
        box()
      }
      points(1:length(selprop_female), selprop_female,
             col = ifelse(selprop_female >= hat_pi2, yes = darken(mycolours, amount = 0.3), no = lighten(mycolours, amount = 0.1)),
             type = "h", lwd = mylwd
      )
      abline(h = hat_pi2, lty = 2, col = "darkred")
      abline(v = xseqgreysep, lty = 1, col = "grey90", lwd = 0.5)
      box()
      abline(v = c(seq(1, length(selprop_male))[!duplicated(mycat)] - 0.5, length(selprop_male) + 0.5), lty = 2)
      par(xpd = TRUE)
      
      catseq <- c(seq(1, length(selprop_male))[!duplicated(mycat)] - 0.5, length(selprop_male) + 0.5)
      for (k in 1:(length(catseq) - 1)) {
        axis(side = 3, at = catseq[c(k, k + 1)] + c(0.2, -0.2), line = 0.7, labels = NA, col = unique(darken(mycolours, amount = 0.4))[k])
      }
      for (k in 1:(length(catseq) - 1)) {
        if (model_id %in% c(1, 3)) {
          axis(
            side = 3, at = mean(catseq[c(k, k + 1)]), line = 0.5, labels = unique(mycat)[k], tick = FALSE, cex.axis = 1.4,
            col.axis = unique(darken(mycolours, amount = 0.4))[k]
          )
        } else {
          tmp_cat_name <- unique(mycat)
          # tmp_cat_name[tmp_cat_name %in% c("PCE/QRISK3", "Genetic")] <- NA
          axis(
            side = 3, at = mean(catseq[c(k, k + 1)]), line = 0.5, labels = tmp_cat_name[k], tick = FALSE, cex.axis = 1.4,
            col.axis = unique(darken(mycolours, amount = 0.4))[k]
          )
        }
      }
      par(xpd = FALSE)
      
      plot(1:length(selprop_female), rep(0, length(selprop_female)),
           xlab = "", ylab = "",
           ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n", pch = NA
      )
      abline(v = xseqgreysep, lty = 1, col = "grey90", lwd = 0.5)
      
      plot(1:length(selprop_male), -selprop_male,
           col = ifelse(selprop_male >= hat_pi1, yes = darken(mycolours, amount = 0.3), no = lighten(mycolours, amount = 0.1)),
           type = "h", lwd = mylwd, xlab = "", ylab = expression(atop(NA, atop(
             textstyle("Selection Proportion"),
             textstyle("(Men)")
           ))),
           xaxt = "n", cex.lab = 1.3, yaxt = "n", panel.first = c(abline(v = 1:length(selprop_female), lty = 3, col = "grey"))
      )
      par(xpd = TRUE)
      for (k in 1:length(selprop_male)) {
        axis(
          side = 3, line = 9, at = k, labels = eval(parse(text = paste0("expression('", names(selprop_male)[k], "')"))), las = 2,
          col.axis = darken(mycolours, amount = 0.3)[k],
          hadj = 0.5, tick = FALSE
        )
      }
      par(xpd = FALSE)
      xseqgreysep <- c(1:(length(selprop_female) + 1)) - 0.5
      if (background) {
        ymax <- max(selprop_female, na.rm = TRUE) + 10
        for (k in seq(1, length(xseqgreysep), by = 2)) {
          polygon(
            x = c(xseqgreysep[k], xseqgreysep[k + 1], xseqgreysep[k + 1], xseqgreysep[k]),
            y = c(-ymax, -ymax, ymax, ymax), col = lighten(mycolours[mycat[k]], amount = 0.97), border = NA
          )
        }
        for (k in seq(2, length(xseqgreysep), by = 2)) {
          polygon(
            x = c(xseqgreysep[k], xseqgreysep[k + 1], xseqgreysep[k + 1], xseqgreysep[k]),
            y = c(-ymax, -ymax, ymax, ymax), col = lighten(mycolours[mycat[k]], amount = 0.99), border = NA
          )
        }
        box()
      }
      points(1:length(selprop_male), -selprop_male,
             col = ifelse(selprop_male >= hat_pi1, yes = darken(mycolours, amount = 0.3), no = lighten(mycolours, amount = 0.1)),
             type = "h", lwd = mylwd, xlab = "", ylab = "Selection Proportion (LASSO)", xaxt = "n", cex.lab = 1.2
      )
      abline(h = -hat_pi1, lty = 2, col = "darkred")
      abline(v = xseqgreysep, lty = 1, col = "grey90", lwd = 0.5)
      box()
      axis(2, at = axTicks(2), labels = -axTicks(2))
      abline(h = hat_pi1, col = "darkred", lty = 3)
      abline(v = c(seq(1, length(selprop_male))[!duplicated(mycat)] - 0.5, length(selprop_male) + 0.5), lty = 2)
      dev.off()
    }
    system(paste("pdfcrop --margin 10", plotname, plotname))
  }
}
