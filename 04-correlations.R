rm(list = ls())

library(openxlsx)
library(sharp)
library(colorspace)

outcome <- "cvd"

# Loading the data
mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_updated.rds")))

# Variable annotation
mydict <- read.xlsx("Data/Data_dictionary.xlsx")

# Extracting ever-used variables
mydict <- mydict[which(apply(mydict[, 4:11], 1, FUN = function(x) {
  any(x == "X")
})), ]
rownames(mydict) <- mydict[, 1]
x <- mydata[, rownames(mydict)]
colnames(x) <- mydict[colnames(x), "Description"]

# Defining categories
mycat <- mydict$Category
mycat[mycat == "Established Risk Factors"] <- "PCE/QRISK3"

# Defining colours by category
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

# Figure
{
  pdf("Figures/Correlation_heatmap_predictors.pdf",
    width = 9, height = 9
  )
  par(mar = c(11, 11, 1, 1))

  # Controls
  Heatmap(cor(na.exclude(x)),
    legend_range = c(-1, 1),
    axes = FALSE, legend = FALSE,
    col = c("royalblue", "white", "red")
  )
  axis(side = 1, at = c(0.5, ncol(x) - 0.5), labels = NA)
  axis(side = 2, at = c(0.5, ncol(x) - 0.5), labels = NA)
  for (k in 1:ncol(x)) {
    axis(
      side = 1, at = k - 0.5, labels = colnames(x)[k], cex.axis = 0.5,
      col.axis = darken(mycolours[k], amount = 0.4),
      col = darken(mycolours[k], amount = 0.4), las = 2
    )
    axis(
      side = 2, at = k - 0.5, labels = colnames(x)[ncol(x) - k + 1], cex.axis = 0.5,
      col.axis = darken(mycolours[ncol(x) - k + 1], amount = 0.4),
      col = darken(mycolours[ncol(x) - k + 1], amount = 0.4), las = 1
    )
  }
  abline(v = c(which(!duplicated(mydict$Category)) - 1, ncol(x)), lty = 2)
  abline(h = ncol(x) - c(which(!duplicated(mydict$Category)) - 1, ncol(x)), lty = 2)
  dev.off()
}

# # Stratified by future disease status
# {
#   pdf("Figures/Correlation_heatmap_predictors_stratified.pdf",
#     width = 18, height = 9
#   )
#   par(mar = c(11, 11, 1, 4), mfrow = c(1, 2))
#
#   # Controls
#   Heatmap(cor(na.exclude(x[which(mydata$case == 0), ])),
#     legend_range = c(-1, 1),
#     axes = FALSE, legend = FALSE,
#     col = c("royalblue", "white", "red")
#   )
#   axis(side = 1, at = c(0.5, ncol(x) - 0.5), labels = NA)
#   axis(side = 2, at = c(0.5, ncol(x) - 0.5), labels = NA)
#   for (k in 1:ncol(x)) {
#     axis(
#       side = 1, at = k - 0.5, labels = colnames(x)[k], cex.axis = 0.5,
#       col.axis = darken(mycolours[k], amount = 0.4),
#       col = darken(mycolours[k], amount = 0.4), las = 2
#     )
#     axis(
#       side = 2, at = k - 0.5, labels = colnames(x)[ncol(x) - k + 1], cex.axis = 0.5,
#       col.axis = darken(mycolours[ncol(x) - k + 1], amount = 0.4),
#       col = darken(mycolours[ncol(x) - k + 1], amount = 0.4), las = 1
#     )
#   }
#   abline(v = c(which(!duplicated(mydict$Category)) - 1, ncol(x)), lty = 2)
#   abline(h = ncol(x) - c(which(!duplicated(mydict$Category)) - 1, ncol(x)), lty = 2)
#   mtext(text = "A", side = 2, line = 9, at = ncol(x) + 3, cex = 3, las = 1)
#
#   # Cases
#   Heatmap(cor(na.exclude(x[which(mydata$case == 1), ])),
#     legend_range = c(-1, 1), legend_length = 15,
#     axes = FALSE,
#     col = c("royalblue", "white", "red")
#   )
#   axis(side = 1, at = c(0.5, ncol(x) - 0.5), labels = NA)
#   axis(side = 2, at = c(0.5, ncol(x) - 0.5), labels = NA)
#   for (k in 1:ncol(x)) {
#     axis(
#       side = 1, at = k - 0.5, labels = colnames(x)[k], cex.axis = 0.5,
#       col.axis = darken(mycolours[k], amount = 0.4),
#       col = darken(mycolours[k], amount = 0.4), las = 2
#     )
#     axis(
#       side = 2, at = k - 0.5, labels = colnames(x)[ncol(x) - k + 1], cex.axis = 0.5,
#       col.axis = darken(mycolours[ncol(x) - k + 1], amount = 0.4),
#       col = darken(mycolours[ncol(x) - k + 1], amount = 0.4), las = 1
#     )
#   }
#   abline(v = c(which(!duplicated(mydict$Category)) - 1, ncol(x)), lty = 2)
#   abline(h = ncol(x) - c(which(!duplicated(mydict$Category)) - 1, ncol(x)), lty = 2)
#   mtext(text = "B", side = 2, line = 9, at = ncol(x) + 3, cex = 3, las = 1)
#   dev.off()
# }
