rm(list = ls())

outcome <- "cvd"
data_input <- "updated"

for (model_id in c(1, 3)) {
  # Variable annotation
  dict <- read.xlsx("Data/Data_dictionary.xlsx")
  rownames(dict) <- dict[, 1]
  variable_cat <- dict[, 3]
  names(variable_cat) <- dict[, 2]
  mycat <- variable_cat

  {
    pdf(paste0("Figures/Stability_sensitivity_m", model_id, "_", outcome, "_", data_input, ".pdf"),
      useDingbats = FALSE, width = 14, height = 14
    )
    par(mar = c(17, 5.5, 2, 1), mfrow = c(2, 1))

    for (gender in c("male", "female")) {
      bigselprop <- readRDS(paste0("Results/HPC_results/Sensitivity/stability_proportions_cvd_m", model_id, "_updated_", gender, "_merged.rds"))
      bigselected <- readRDS(paste0("Results/HPC_results/Sensitivity/stability_selected_cvd_m", model_id, "_updated_", gender, "_merged.rds"))
      selected <- formatC(apply(bigselected, 1, FUN = function(x) {
        sum(x > 0) / length(x)
      }), format = "f", digits = 2)
      rownames(bigselprop) <- dict[rownames(bigselprop), 2]
      tmp <- as.vector(t(cbind(bigselprop, matrix(NA, nrow = nrow(bigselprop), ncol = 3))))

      # Change group name
      if (outcome == "cvd") {
        mycat[mycat == "Established Risk Factors"] <- "PCE/QRISK3"
      } else {
        mycat[mycat == "Established Risk Factors"] <- "PCE/QRISK3"
      }
      mycat <- mycat[rownames(bigselprop)]

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

      # English/American spelling
      id <- grep("Haemoglobin", rownames(bigselprop))
      rownames(bigselprop)[id] <- "Hemoglobin"
      id <- grep("Haematocrit", rownames(bigselprop))
      rownames(bigselprop)[id] <- "Hematocrit"

      plot(tmp,
        type = "h", xaxt = "n", xlab = "", ylab = "Selection Proportion", cex.lab = 1.5, lwd = 0.5,
        las = 1, col = rep(mycolours, each = ncol(bigselprop) + 3)
      )
      if (model_id == 1) {
        mtext(text = ifelse(gender == "male", yes = "A", no = "B"), side = 2, at = 1.1, las = 1, cex = 3, line = 2.7)
      } else {
        mtext(text = ifelse(gender == "male", yes = "C", no = "D"), side = 2, at = 1.1, las = 1, cex = 3, line = 2.7)
      }
      abline(v = seq(-1, length(tmp), by = ncol(bigselprop) + 3), lty = 3, col = "grey")
      axis(
        side = 1, at = seq(mean(c(-1, ncol(bigselprop) + 3)), length(tmp), by = ncol(bigselprop) + 3),
        labels = rownames(bigselprop), las = 2, cex.axis = 0.8
      )
      for (i in 1:length(selected)) {
        axis(
          side = 3, at = seq(mean(c(-1, ncol(bigselprop) + 3)), length(tmp), by = ncol(bigselprop) + 3)[i],
          las = 2, cex.axis = 0.8, tick = FALSE, line = -0.5,
          labels = ifelse(as.numeric(selected[i]) >= 0.1,
            yes = eval(parse(text = paste0("expression(bold('", selected[i], "'))"))),
            no = eval(parse(text = paste0("expression('", selected[i], "')")))
          )
        )
      }
    }
    dev.off()
  }
}
