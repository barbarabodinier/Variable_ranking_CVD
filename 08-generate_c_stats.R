rm(list = ls())

library(data.table)
library(glmnet)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)
library(c060)
library(plotrix)

source("Scripts/functions.R")

# outcomes=c("_cad", "_cvd", "_cad_pce", "_cvd_pce")
outcomes <- c("cvd")
genders <- c("male", "female")
mymethod <- "lasso"
data_input <- "updated"

# Variable annotation
mydict <- read.xlsx("Data/Data_dictionary.xlsx")

# Extracting ever-used variables
mydict <- mydict[which(apply(mydict[, 4:ncol(mydict)], 1, FUN = function(x) {
  any(x == "X")
})), ]
rownames(mydict) <- mydict[, 1]
variable_cat <- mydict[, 3]
names(variable_cat) <- mydict[, 2]


for (i in 1:length(outcomes)) {
  outcome <- outcomes[i]
  print(outcome)

  for (model_id in c(1:4,7)) {
    # Initialising empty matrices
    perf <- matrix(NA, nrow = 2, ncol = 3)
    N <- rep(NA, 3)
    names(N) <- c("merged", "male", "female")
    hr <- matrix(NA, nrow = nrow(mydict), ncol = 2)
    rownames(hr) <- rownames(mydict)
    colnames(hr) <- c("male", "female")

    for (gender in c("male", "female", "merged")) {

      # Loading the data
      mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_", data_input, ".rds")))

      # Extracting training set
      eids0 <- readRDS("Data/Split/performance_set_eids_0_1.rds")
      eids1 <- readRDS("Data/Split/performance_set_eids_1_1.rds")
      mydata_test <- mydata[c(eids0, eids1), ]

      # Extracting selection set
      eids0 <- readRDS("Data/Split/estimation_set_eids_0_1.rds")
      eids1 <- readRDS("Data/Split/estimation_set_eids_1_1.rds")
      mydata <- mydata[c(eids0, eids1), ]


      ### PCE

      cox0 <- readRDS(paste0("Results/Cox_models/cox_pce_female_cvd.rds"))
      cox1 <- readRDS(paste0("Results/Cox_models/cox_pce_male_cvd.rds"))


      # ## Training
      #
      # # Get probabilities
      # if (gender=="merged"){
      # myS=c(cox0$S_training, cox1$S_training)
      # }
      # if (gender=="female"){
      #   myS=c(cox0$S_training)
      # }
      # if (gender=="male"){
      #   myS=c(cox1$S_training)
      # }
      # mydata=mydata[names(myS),]
      # print(all(names(myS)==rownames(mydata)))
      #
      # # Compute survival time
      # survobject=Surv(mydata$time, mydata$case)
      #
      # # Concordance over full sample
      # perf[1,1]=GetConcordance(survobject=survobject, S=myS)


      ## Test

      # Get probabilities
      if (gender == "merged") {
        myS <- c(cox0$S_test, cox1$S_test)
      }
      if (gender == "female") {
        myS <- c(cox0$S_test)
      }
      if (gender == "male") {
        myS <- c(cox1$S_test)
      }
      mydata_test <- mydata_test[names(myS), ]
      N[gender] <- paste(formatC(cumsum(rev(table(mydata_test$case))), big.mark = ","), collapse = "/")
      print(all(names(myS) == rownames(mydata_test)))

      # Compute survival time
      survobject <- Surv(mydata_test$time, mydata_test$case)

      # Concordance over full sample
      if (gender == "merged") {
        perf[1, 1] <- GetConcordance(survobject = survobject, S = myS)
      }
      if (gender == "female") {
        perf[1, 3] <- GetConcordance(survobject = survobject, S = myS)
      }
      if (gender == "male") {
        perf[1, 2] <- GetConcordance(survobject = survobject, S = myS)
      }

      ### LASSO stability selection

      cox0 <- readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_female.rds"))
      cox1 <- readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_male.rds"))

      if (gender == "merged") {
        tmp <- formatC(summary(cox0$cox_model)$coefficients[, 2], format = "f", digits = 2)
        hr[names(tmp), "female"] <- tmp
        tmp <- formatC(summary(cox1$cox_model)$coefficients[, 2], format = "f", digits = 2)
        hr[names(tmp), "male"] <- tmp
        hr <- hr[which(apply(hr, 1, FUN = function(x) {
          sum(is.na(x)) < 2
        })), ]
        rownames(hr) <- mydict[rownames(hr), 2]
        hr <- hr[sort.list(rownames(hr)), ]
        write.xlsx(as.data.frame(hr), paste0("Tables/HR_", outcome, "_m", model_id, "_", data_input, ".xlsx"),
          row.names = TRUE, overwrite = TRUE
        )
      }

      # ## Training
      #
      # # Get probabilities
      # myS=c(cox0$S_training, cox1$S_training)
      # print(all(names(myS)==rownames(mydata)))
      #
      # # Compute survival time
      # survobject=Surv(mydata$time, mydata$case)
      #
      # # Concordance over full sample
      # perf[2,1]=GetConcordance(survobject=survobject, S=myS)

      ## Test

      # Get probabilities
      if (gender == "merged") {
        myS <- c(cox0$S_test, cox1$S_test)
      }
      if (gender == "female") {
        myS <- c(cox0$S_test)
      }
      if (gender == "male") {
        myS <- c(cox1$S_test)
      }
      mydata_test <- mydata_test[names(myS), ]
      print(all(names(myS) == rownames(mydata_test)))

      # Compute survival time
      survobject <- Surv(mydata_test$time, mydata_test$case)

      # Concordance over full sample
      if (gender == "merged") {
        perf[2, 1] <- GetConcordance(survobject = survobject, S = myS)
      }
      if (gender == "female") {
        perf[2, 3] <- GetConcordance(survobject = survobject, S = myS)
      }
      if (gender == "male") {
        perf[2, 2] <- GetConcordance(survobject = survobject, S = myS)
      }


      # ### Random Forest
      #
      # cox0=readRDS(paste0("Results/cox_model_vi_random_forest_",0,outcome,".rds"))
      # cox1=readRDS(paste0("Results/cox_model_vi_random_forest_",1,outcome,".rds"))
      #
      #
      # ## Training
      #
      # # Get probabilities
      # myS=c(cox0$S_training, cox1$S_training)
      # print(all(names(myS)==rownames(mydata)))
      #
      # # Compute survival time
      # survobject=Surv(mydata$time, mydata$case)
      #
      # # Concordance over full sample
      # perf[3,1]=GetConcordance(survobject=survobject, S=myS)
      #
      # ## Test
      #
      # # Get probabilities
      # myS=c(cox0$S_test, cox1$S_test)
      # print(all(names(myS)==rownames(mydata_test)))
      #
      # # Compute survival time
      # survobject=Surv(mydata_test$time, mydata_test$case)
      #
      # # Concordance over full sample
      # perf[3,2]=GetConcordance(survobject=survobject, S=myS)
    }
    ### Saving table
    rownames(perf) <- c("PCE", "LASSO")
    colnames(perf) <- c("Merged", "Men", "Women")
    perf <- rbind(N, perf)
    write.xlsx(as.data.frame(perf), paste0("Tables/C_statistics_", outcome, "_m", model_id, "_", data_input, ".xlsx"),
      row.names = TRUE, col.names = TRUE, overwrite = TRUE
    )

    # Differences in C statistics
    diff <- rbind(c(
      DeltaCStat(perf[2, 1], perf[3, 1]),
      DeltaCStat(perf[2, 2], perf[3, 2]),
      DeltaCStat(perf[2, 3], perf[3, 3])
    ))
    colnames(diff) <- c("Merged", "Men", "Women")
    diff <- rbind(N, diff)
    write.xlsx(as.data.frame(diff), paste0("Tables/Difference_c_statistics_", outcome, "_m", model_id, "_", data_input, ".xlsx"),
      overwrite = TRUE
    )
  }
}
