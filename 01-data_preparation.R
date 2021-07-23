rm(list = ls())

library(openxlsx)
library(pheatmap)

# Loading data dictionary
dict <- read.xlsx("Data/Data_dictionary.xlsx")
rownames(dict) <- dict[, 1]

# Preparing and saving data
for (outcome in c("CAD", "CVD")) {
  # Loading data
  mydata <- readRDS(paste0("Data/", outcome, "_imputed_MW.rds"))
  nmr <- readRDS(paste0("Data/Nightingale_log_transformed.rds"))
  rownames(mydata) <- mydata$eid

  # Final data preparation
  mydata$bmi <- mydata$weight / (mydata$height / 100)^2
  mydata$smoking <- as.factor(mydata$smoking)
  mydata$ethnicity_basic <- factor(mydata$ethnicity_basic, levels = c("White", "Black", "Other"))

  # Log-transformation of bmk and haem variables
  var_list <- dict[which(dict[, 3] %in% c("Biochemistry", "Haematology")), 1]
  for (i in var_list) { # setting small value instead of zero
    mydata[which(mydata[, i] == 0), i] <- NA
    mydata[which(is.na(mydata[, i])), i] <- mean(c(0, min(mydata[, i], na.rm = TRUE)))
  }
  mydata[, c(var_list)] <- log(mydata[, c(var_list)])

  # Merging with NMR data
  ids <- intersect(rownames(nmr), rownames(mydata))
  tmp <- nmr[ids, ]
  add <- rownames(mydata)[!rownames(mydata) %in% ids]
  empty_mat <- matrix(NA, nrow = length(add), ncol = ncol(tmp))
  colnames(empty_mat) <- colnames(tmp)
  rownames(empty_mat) <- add
  tmp <- rbind(tmp, empty_mat)
  tmp <- tmp[rownames(mydata), ]
  print(all(rownames(tmp) == rownames(mydata)))
  mydata <- cbind(mydata, tmp)

  # Identifying columns to keep
  predictors <- dict[, 1]
  predictors[predictors == "ethnicity_basicBlack"] <- "ethnicity_basic"
  predictors <- predictors[!predictors %in% "ethnicity_basicOther"]
  predictors[predictors == "smoking1"] <- "smoking"
  predictors <- predictors[!predictors %in% "smoking2"]
  print(all(predictors %in% colnames(mydata)))

  # Recoding
  options(na.action = "na.pass") # ensuring NAs are kept in the resulting dataset
  mydata <- model.matrix(as.formula(paste0(
    "~", paste0(predictors, collapse = "+"), "+",
    paste0(c("sex", "time", "case"), collapse = "+")
  )), data = mydata)[, -1]
  mydata <- data.frame(mydata)
  options(na.action = "na.exclude") # re-setting default NA handling

  saveRDS(mydata, paste0("Data/", outcome, "_imputed_MW_full.rds"))
}

# Correlations between selected metabolites
selected_nmr <- na.exclude(mydata)
selected_nmr <- selected_nmr[, grepl("^X", colnames(selected_nmr))]
selected_nmr <- selected_nmr[, intersect(dict[which(dict$`Model.3.(Male)` == "X"), 1], colnames(selected_nmr))]
mycor <- cor(selected_nmr)
rownames(mycor) <- colnames(mycor) <- dict[rownames(mycor), 2]
pheatmap(mycor,
  breaks = seq(-1, 1, length.out = 100),
  border = NA, width = 12, height = 12,
  filename = paste0("NMR_descriptive/Figures/Correlations_selected.pdf")
)

# Correlations between selected metabolites and biochem/haem variables
mycor <- cor(selected_nmr, mydata[rownames(selected_nmr), dict[which(dict[, 3] %in% c("Biochemistry", "Haematology")), 1]])
rownames(mycor) <- dict[rownames(mycor), 2]
colnames(mycor) <- dict[colnames(mycor), 2]
pheatmap(mycor,
  breaks = seq(-1, 1, length.out = 100),
  border = NA, width = 12, height = 7,
  filename = paste0("NMR_descriptive/Figures/Correlations_selected_nmr_biochem_haem.pdf")
)

# Correlations between biochemistry
mycor <- cor(
  mydata[, dict[which(dict[, 3] %in% c("Biochemistry")), 1]],
  mydata[, dict[which(dict[, 3] %in% c("Biochemistry")), 1]]
)
rownames(mycor) <- colnames(mycor) <- dict[rownames(mycor), 2]
pheatmap(mycor,
  breaks = seq(-1, 1, length.out = 100),
  border = NA, width = 12, height = 12,
  filename = paste0("NMR_descriptive/Figures/Correlations_biochemistry.pdf")
)

# Correlations between biochemistry
mycor <- cor(
  mydata[, dict[which(dict[, 3] %in% c("Haematology")), 1]],
  mydata[, dict[which(dict[, 3] %in% c("Haematology")), 1]]
)
rownames(mycor) <- colnames(mycor) <- dict[rownames(mycor), 2]
pheatmap(mycor,
  breaks = seq(-1, 1, length.out = 100),
  border = NA, width = 12, height = 12,
  filename = paste0("NMR_descriptive/Figures/Correlations_haematology.pdf")
)

# Correlations between biochemistry
mycor <- cor(cbind(
  selected_nmr,
  mydata[rownames(selected_nmr), dict[which(dict[, 3] %in% c("Biochemistry")), 1]],
  mydata[rownames(selected_nmr), dict[which(dict[, 3] %in% c("Haematology")), 1]]
))
rownames(mycor) <- colnames(mycor) <- dict[rownames(mycor), 2]
pheatmap(mycor,
  breaks = seq(-1, 1, length.out = 100),
  border = NA, width = 15, height = 15,
  filename = paste0("NMR_descriptive/Figures/Correlations_all.pdf")
)


for (outcome_name in c("CAD", "CVD")) {
  mydata <- readRDS(paste0("Data/", outcome_name, "_imputed_MW_full.rds"))
  print(colnames(mydata))

  diet <- readRDS("Data/Raw_diet/first_diet_df.RDS")
  rownames(diet) <- diet$eid
  diet <- diet[, -1]
  table(rownames(diet) %in% rownames(mydata))
  mydata <- cbind(mydata, diet[rownames(mydata), ])
  print(colnames(mydata))

  saveRDS(mydata, paste0("Data/", outcome_name, "_imputed_MW_full.rds"))
}
