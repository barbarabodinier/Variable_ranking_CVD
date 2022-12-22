rm(list = ls())

library(openxlsx)
library(abind)

dir.create("Data/Split", showWarnings = FALSE)

for (seed in 0:100) {
  print(seed)
  for (gender in c(0, 1)) { # male is 1
    print(c("Female", "Male")[gender + 1])

    # Loading the data
    mydata <- data.frame(readRDS("Data/CVD_imputed_MW_updated.rds")) # TODO: make use of CAD and CVD
    mydict <- read.xlsx("Data/Data_dictionary.xlsx")

    # Stratified by gender
    mydata <- mydata[mydata$sex == gender, ]

    # Selection set
    tau <- 0.4
    eids_no_nmr <- rownames(mydata)[is.na(mydata$X23400)]
    eids_nmr <- rownames(mydata)[!is.na(mydata$X23400)]
    table(rownames(mydata) %in% c(eids_no_nmr, eids_nmr))
    set.seed(seed)
    s0 <- sample(eids_no_nmr[which(mydata[eids_no_nmr, "case"] == 0)],
      size = tau * sum(mydata[eids_no_nmr, "case"] == 0)
    )
    s1 <- sample(eids_no_nmr[which(mydata[eids_no_nmr, "case"] == 1)],
      size = tau * sum(mydata[eids_no_nmr, "case"] == 1)
    )
    s2 <- sample(eids_nmr[which(mydata[eids_nmr, "case"] == 0)],
      size = tau * sum(mydata[eids_nmr, "case"] == 0)
    )
    s3 <- sample(eids_nmr[which(mydata[eids_nmr, "case"] == 1)],
      size = tau * sum(mydata[eids_nmr, "case"] == 1)
    )
    selection_eids <- c(s0, s1, s2, s3)
    print(table(mydata[selection_eids, "case"]))
    saveRDS(selection_eids, paste0("Data/Split/selection_set_eids_", gender, "_", seed, ".rds"))

    # Estimation set
    mydata <- mydata[!as.character(rownames(mydata)) %in% selection_eids, ]
    tau <- 0.5
    eids_no_nmr <- rownames(mydata)[is.na(mydata$X23400)]
    eids_nmr <- rownames(mydata)[!is.na(mydata$X23400)]
    table(rownames(mydata) %in% c(eids_no_nmr, eids_nmr))
    set.seed(seed)
    s0 <- sample(eids_no_nmr[which(mydata[eids_no_nmr, "case"] == 0)],
      size = tau * sum(mydata[eids_no_nmr, "case"] == 0)
    )
    s1 <- sample(eids_no_nmr[which(mydata[eids_no_nmr, "case"] == 1)],
      size = tau * sum(mydata[eids_no_nmr, "case"] == 1)
    )
    s2 <- sample(eids_nmr[which(mydata[eids_nmr, "case"] == 0)],
      size = tau * sum(mydata[eids_nmr, "case"] == 0)
    )
    s3 <- sample(eids_nmr[which(mydata[eids_nmr, "case"] == 1)],
      size = tau * sum(mydata[eids_nmr, "case"] == 1)
    )
    estimation_eids <- c(s0, s1, s2, s3)
    print(table(mydata[estimation_eids, "case"]))
    saveRDS(estimation_eids, paste0("Data/Split/estimation_set_eids_", gender, "_", seed, ".rds"))

    # Performance set
    mydata <- mydata[!as.character(rownames(mydata)) %in% estimation_eids, ]
    perf_eids <- rownames(mydata)
    print(table(mydata[perf_eids, "case"]))
    saveRDS(perf_eids, paste0("Data/Split/performance_set_eids_", gender, "_", seed, ".rds"))

    # Checking
    print(length(intersect(selection_eids, estimation_eids)))
    print(length(intersect(selection_eids, perf_eids)))
    print(length(intersect(estimation_eids, perf_eids)))
    cat("\n")
  }
}
