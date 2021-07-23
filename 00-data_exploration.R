rm(list = ls())

library(data.table)
library(openxlsx)
library(pheatmap)

# Creating directories
dir.create("NMR_descriptive", showWarnings = FALSE)
dir.create("NMR_descriptive/Figures", showWarnings = FALSE)
dir.create("NMR_descriptive/Figures/Distributions", showWarnings = FALSE)

# Loading existing data
mydata <- readRDS("Data/CVD_imputed_MW.rds")
rownames(mydata) <- mydata$eid

# Loading NMR data
nmr <- readRDS("Data/Raw_nightingale/ukb46678_19266.rds")

# Loading annotation
annot <- fread("Data/Raw_nightingale/Data_Dictionary_Showcase.csv", data.table = FALSE)
rownames(annot) <- annot$FieldID

# Excluding participants with no measurements
prop_missing <- apply(nmr, 1, FUN = function(x) {
  sum(is.na(x)) / length(x)
})
nmr <- nmr[prop_missing < 0.99, ] # only available data is eid
rm(prop_missing)

# Keeping baseline data
nmr <- nmr[, grep("-0", colnames(nmr))]
print(all(gsub(".*-0.", "", colnames(nmr)) == 0))
colnames(nmr) <- gsub("-0.*", "", colnames(nmr))

# Creating list of Nightingale metabolites
print(all(colnames(nmr) %in% rownames(annot)))
annot <- annot[colnames(nmr), ]
write.xlsx(annot, "NMR_descriptive/List_nmr_variables.xlsx")

# Proportions of missing
prop_missing_indiv <- apply(nmr, 1, FUN = function(x) {
  sum(is.na(x)) / length(x)
})
print("Proportions of missing by participant:")
print(range(prop_missing_indiv))
prop_missing <- apply(nmr, 2, FUN = function(x) {
  sum(is.na(x)) / length(x)
})
print("Proportions of missing by variable:")
print(range(prop_missing))

# Updating names
dict <- read.xlsx("Data/Data_dictionary.xlsx")
rownames(dict) <- dict[, 1]
colnames(nmr) <- paste0("X", colnames(nmr))

# Complete cases
nmr <- na.exclude(nmr)

# Distributions of raw levels
for (k in 1:ncol(nmr)) {{ pdf(paste0("NMR_descriptive/Figures/Distributions/Distribution_", colnames(nmr)[k], ".pdf"))
  par(mar = c(5, 3, 1, 1))
  plot(density(nmr[, k]),
    las = 1, lwd = 2,
    xlab = dict[colnames(nmr)[k], 2],
    main = "", ylab = "", cex.lab = 1.5
  )
  dev.off() }}

# Log-transformation of nmr variables
for (i in 1:ncol(nmr)) { # setting small value instead of zero
  nmr[which(nmr[, i] == 0), i] <- NA
  nmr[which(is.na(nmr[, i])), i] <- mean(c(0, min(nmr[, i], na.rm = TRUE)))
}
nmr <- log(nmr)
saveRDS(nmr, "Data/Nightingale_log_transformed.rds")

# Distributions of log-transformed levels
for (k in 1:ncol(nmr)) {{ pdf(paste0("NMR_descriptive/Figures/Distributions/Distribution_log_", colnames(nmr)[k], ".pdf"))
  par(mar = c(5, 3, 1, 1))
  plot(density(nmr[, k]),
    las = 1, lwd = 2,
    xlab = dict[colnames(nmr)[k], 2],
    main = "", ylab = "", cex.lab = 1.5
  )
  dev.off() }}

# Correlations
mycor <- cor(nmr)
rownames(mycor) <- colnames(mycor) <- dict[rownames(mycor), 2]
pheatmap(mycor,
  breaks = seq(-1, 1, length.out = 100),
  border = NA, width = 30, height = 30,
  filename = paste0("NMR_descriptive/Figures/Correlations.pdf")
)

# Log-transformation of biochemistry and haematology data
others_list <- c(dict[which(dict[, 3] %in% c("Biochemistry", "Haematology")), 1], "Cholesterol", "HDL_cholesterol")
others <- mydata[, others_list]
for (i in 1:ncol(others)) { # setting small value instead of zero
  others[which(others[, i] == 0), i] <- NA
  others[which(is.na(others[, i])), i] <- mean(c(0, min(others[, i], na.rm = TRUE)))
}
others <- log(others)
rownames(dict) <- dict[, 1]

# Correlations between NMR and biochem/haem data
ids <- intersect(rownames(others), rownames(nmr))
mycor <- cor(nmr[ids, ], others[ids, ])
rownames(mycor) <- dict[rownames(mycor), 2]
colnames(mycor) <- dict[colnames(mycor), 2]
pheatmap(mycor,
  breaks = seq(-1, 1, length.out = 100),
  border = NA, width = 15, height = 30,
  filename = paste0("NMR_descriptive/Figures/Correlations_nmr_biochem_haem.pdf")
)

# Correlations between lipid related biomarkers
mycor=cor(nmr[,grepl("DL", dict[colnames(nmr),2])])
rownames(mycor)=colnames(mycor)=dict[colnames(mycor), 2]
pheatmap(mycor,
         breaks = seq(-1, 1, length.out = 100),
         border = NA, width = 15, height = 15,
         filename = paste0("NMR_descriptive/Figures/Correlations_nmr_dl.pdf")
)

# Computing highest absolute correlation per metab
max_cor <- apply(abs(mycor), 1, max)
max_cor <- sort(max_cor, decreasing = FALSE)
max_cor <- round(max_cor, digits = 3)
write.table(cbind(names(max_cor), max_cor), "NMR_descriptive/Max_correlation_with_biochem_haem.txt",
  col.names = FALSE, row.names = FALSE
)

max_cor <- max_cor[!grepl("DL", names(max_cor))]
write.table(cbind(names(max_cor), max_cor), "NMR_descriptive/Max_correlation_with_biochem_haem_noDL.txt",
  col.names = FALSE, row.names = FALSE
)

tokeep <- names(max_cor)[max_cor < 0.7] # NB: also manually removing redundant ones
