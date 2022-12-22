rm(list = ls())

library(sharp)
library(survival)
library(openxlsx)
library(plotrix)
library(colorspace)

source("Scripts/functions.R")
source("Scripts/gnd_functions.R")

# Parameters
outcome <- "cvd"
data_input <- "updated"
model_id <- 1

# Variable annotation
dict <- read.xlsx("Data/Data_dictionary.xlsx")
rownames(dict) <- dict[, 1]
variable_cat <- dict[, 3]
names(variable_cat) <- dict[, 2]

{
  pdf(paste0("Figures/Recalibration_PCE_", outcome, ".pdf"), width = 8, height = 16)
  par(mfrow = c(4, 2), mar = c(5, 5, 1, 1))
  for (j in 0:1) {
    gender <- c("male", "female")[j + 1]
    print(gender)

    # Loading the data
    mydata_imputed <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_", data_input, ".rds")))
    mydata_unimputed <- data.frame(readRDS(paste0("Data/Raw_MW/", toupper(outcome), "_unimputed_MW.rds")))
    rownames(mydata_unimputed) <- as.character(mydata_unimputed$eid)
    ids <- intersect(rownames(mydata_imputed), rownames(mydata_unimputed))
    mydata_imputed <- mydata_imputed[ids, ]
    mydata_unimputed <- mydata_unimputed[ids, ]
    mydata <- mydata_unimputed
    mydata$case <- mydata_imputed$case
    mydata$time <- mydata_imputed$time
    mydata <- mydata[, c("logHR", "pce_score", "time", "case")]

    # Extracting training set
    if (gender == "female") {
      eids <- readRDS("Data/Split/performance_set_eids_0_1.rds")
    } else {
      eids <- readRDS("Data/Split/performance_set_eids_1_1.rds")
    }
    mydata_test <- mydata[eids, ]

    # Extracting selection set
    if (gender == "female") {
      eids <- readRDS("Data/Split/estimation_set_eids_0_1.rds")
    } else {
      eids <- readRDS("Data/Split/estimation_set_eids_1_1.rds")
    }
    mydata <- mydata[eids, ]

    # Reading PCE model
    cox_pce <- readRDS(paste0("Results/Refitted_models/cox_pce_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
    print(all(names(cox_pce$S_train) == rownames(mydata)))
    print(all(names(cox_pce$S_test) == rownames(mydata_test)))


    ### Training

    # Not recalibrated
    survobject <- Surv(mydata$time, mydata$case)
    out <- GetObsCumInc(survobject = survobject, linear.predictors = mydata$logHR)
    S_obs_sp <- out$S_obs
    S_pred <- unlist(lapply(split(1 - mydata$pce_score, f = out$deciles), mean))
    beta <- CalibPlot(S_pred = S_pred, S_obs = S_obs_sp, xylim = 0.4, name = "")
    GND <- GetGNDStat(survobject = survobject, predicted_S = 1 - mydata$pce_score, mydata = mydata)
    if (GND[3] == 0) {
      legend("topleft",
        cex = 1.2, bty = "n",
        legend = c(
          eval(parse(text = paste0('expression(beta*"="*', beta, ")")))
          # eval(parse(text = paste0('expression(p[GND]*"<"*10^-20)')))
        )
      )
    } else {
      p <- GND[3]
      p <- formatC(p, format = "e", digits = 2)
      p <- gsub("e", "*'x'*10^", p, fixed = FALSE)
      legend("topleft",
        cex = 1.2, bty = "n",
        legend = c(
          eval(parse(text = paste0('expression(beta*"="*', beta, ")")))
          # eval(parse(text = paste0('expression(p[GND]*"<"*', p, ")")))
        )
      )
    }
    par(xpd = NA)
    if (gender == 1) {
      text("A", x = -0.08, y = 0.4, cex = 4)
    } else {
      text("E", x = -0.08, y = 0.4, cex = 4)
    }
    par(xpd = FALSE)

    # Recalibrated
    survobject <- Surv(mydata$time, mydata$case)
    out <- GetObsCumInc(survobject = survobject, linear.predictors = mydata$logHR)
    S_obs_sp <- out$S_obs
    S_pred <- unlist(lapply(split(cox_pce$S_train, f = out$deciles), mean))
    beta <- CalibPlot(S_pred = S_pred, S_obs = S_obs_sp, xylim = 0.4, name = "")
    GND <- GetGNDStat(survobject = survobject, predicted_S = 1 - cox_pce$S_train, mydata = mydata)
    if (GND[3] == 0) {
      legend("topleft",
        cex = 1.2, bty = "n",
        legend = c(
          eval(parse(text = paste0('expression(beta*"="*', beta, ")")))
          # eval(parse(text = paste0('expression(p[GND]*"<"*10^-20)')))
        )
      )
    } else {
      p <- GND[3]
      p <- formatC(p, format = "e", digits = 2)
      p <- gsub("e", "*'x'*10^", p, fixed = FALSE)
      legend("topleft",
        cex = 1.2, bty = "n",
        legend = c(
          eval(parse(text = paste0('expression(beta*"="*', beta, ")")))
          # eval(parse(text = paste0('expression(p[GND]*"<"*', p, ")")))
        )
      )
    }
    par(xpd = NA)
    if (gender == 1) {
      text("B", x = -0.08, y = 0.4, cex = 4)
    } else {
      text("F", x = -0.08, y = 0.4, cex = 4)
    }
    par(xpd = FALSE)

    ### Test

    # Not recalibrated
    survobject <- Surv(mydata_test$time, mydata_test$case)
    out <- GetObsCumInc(survobject = survobject, linear.predictors = mydata_test$logHR)
    S_obs_sp <- out$S_obs
    S_pred <- unlist(lapply(split(1 - mydata_test$pce_score, f = out$deciles), mean))
    beta <- CalibPlot(S_pred = S_pred, S_obs = S_obs_sp, xylim = 0.4, name = "")
    GND <- GetGNDStat(survobject = survobject, predicted_S = 1 - mydata_test$pce_score, mydata = mydata_test)
    if (GND[3] == 0) {
      legend("topleft",
        cex = 1.2, bty = "n",
        legend = c(
          eval(parse(text = paste0('expression(beta*"="*', beta, ")")))
          # eval(parse(text = paste0('expression(p[GND]*"<"*10^-20)')))
        )
      )
    } else {
      p <- GND[3]
      p <- formatC(p, format = "e", digits = 2)
      p <- gsub("e", "*'x'*10^", p, fixed = FALSE)
      legend("topleft",
        cex = 1.2, bty = "n",
        legend = c(
          eval(parse(text = paste0('expression(beta*"="*', beta, ")")))
          # eval(parse(text = paste0('expression(p[GND]*"<"*', p, ")")))
        )
      )
    }
    par(xpd = NA)
    if (gender == 1) {
      text("C", x = -0.08, y = 0.4, cex = 4)
    } else {
      text("G", x = -0.08, y = 0.4, cex = 4)
    }
    par(xpd = FALSE)

    # Recalibrated
    survobject <- Surv(mydata_test$time, mydata_test$case)
    out <- GetObsCumInc(survobject = survobject, linear.predictors = mydata_test$logHR)
    S_obs_sp <- out$S_obs
    S_pred <- unlist(lapply(split(cox_pce$S_test, f = out$deciles), mean))
    beta <- CalibPlot(S_pred = S_pred, S_obs = S_obs_sp, xylim = 0.4, name = "")
    GND <- GetGNDStat(survobject = survobject, predicted_S = 1 - cox_pce$S_test, mydata = mydata_test)
    if (GND[3] == 0) {
      legend("topleft",
        cex = 1.2, bty = "n",
        legend = c(
          eval(parse(text = paste0('expression(beta*"="*', beta, ")")))
          # eval(parse(text = paste0('expression(p[GND]*"<"*10^-20)')))
        )
      )
    } else {
      p <- GND[3]
      p <- formatC(p, format = "e", digits = 2)
      p <- gsub("e", "*'x'*10^", p, fixed = FALSE)
      legend("topleft",
        cex = 1.2, bty = "n",
        legend = c(
          eval(parse(text = paste0('expression(beta*"="*', beta, ")")))
          # eval(parse(text = paste0('expression(p[GND]*"<"*', p, ")")))
        )
      )
    }
    par(xpd = NA)
    if (gender == 1) {
      text("D", x = -0.08, y = 0.4, cex = 4)
    } else {
      text("H", x = -0.08, y = 0.4, cex = 4)
    }
    par(xpd = FALSE)
  }
  dev.off()
}
