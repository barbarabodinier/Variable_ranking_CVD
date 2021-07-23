rm(list = ls())

library(focus)
library(data.table)
library(glmnet)
library(pheatmap)
library(survival)
library(abind)
library(openxlsx)
library(c060)
library(plotrix)
library(PredictABEL)
library(Hmisc)
library(colorspace)
source("Scripts/reclassification_PredictABEL_function.R")
source("Scripts/functions.R")
source("Scripts/roc_functions.R")

dir.create("Figures", showWarnings = FALSE)
dir.create("Results", showWarnings = FALSE)
dir.create("Figures/Exploration", showWarnings = FALSE)

outcomes <- c("cvd")
model_id <- 1
data_input <- "updated"
nbreaks <- 500
legend_letters=c("A","B","C")
names(legend_letters)=c("merged", "male", "female")

for (outcome in outcomes) {
  print(outcome)
  
  for (model_id in 1) {
    {
      pdf(paste0("Figures/ROC_curve_", outcome, "_m", model_id, "_", data_input, ".pdf"),
          width=15, height=5)
      par(mar = c(5, 5, 1, 1), mfrow=c(1,3))
      for (gender in c("merged", "male", "female")) {
        
        # Loading the data
        mydata <- data.frame(readRDS("Data/CVD_imputed_MW_updated.rds"))
        
        eids_test <- c(
          as.character(readRDS(paste0("Data/Split/performance_set_eids_", 0, "_1.rds"))),
          as.character(readRDS(paste0("Data/Split/performance_set_eids_", 1, "_1.rds")))
        )
        mydata_test <- mydata[eids_test, ]
        
        # Loading the models
        cox_stable0 <- readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_female.rds"))
        cox_stable1 <- readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_male.rds"))
        
        if (gender == "merged") {
          S_lasso <- c(cox_stable0$S_test, cox_stable1$S_test)
        }
        if (gender == "female") {
          S_lasso <- c(cox_stable0$S_test)
        }
        if (gender == "male") {
          S_lasso <- c(cox_stable1$S_test)
        }
        mydata_test <- mydata_test[names(S_lasso), ]
        
        # cox_stable0=readRDS(paste0("Results/cox_model_vi_random_forest_",0,outcome,".rds"))
        # cox_stable1=readRDS(paste0("Results/cox_model_vi_random_forest_",1,outcome,".rds"))
        # S_rf=c(cox_stable0$S_test, cox_stable1$S_test)
        
        cox_pce0 <- readRDS(paste0("Results/Cox_models/cox_pce_female_", outcome, ".rds"))
        cox_pce1 <- readRDS(paste0("Results/Cox_models/cox_pce_male_", outcome, ".rds"))
        
        if (gender == "merged") {
          S_pce <- c(cox_pce0$S_test, cox_pce1$S_test)
        }
        if (gender == "female") {
          S_pce <- c(cox_pce0$S_test)
        }
        if (gender == "male") {
          S_pce <- c(cox_pce1$S_test)
        }
        
        print(all(names(S_pce) == rownames(mydata_test)))
        print(all(names(S_lasso) == rownames(mydata_test)))
        # print(all(names(S_rf)==rownames(mydata_test)))
        
        # Compute sensitivity, specificity, AUC
        roc_pce <- GetROC(y = mydata_test$case, x = 1 - S_pce, nbreaks = nbreaks)
        roc_lasso <- GetROC(y = mydata_test$case, x = 1 - S_lasso, nbreaks = nbreaks)
        # roc_rf=GetROC(y=mydata_test$case, x=1-S_rf, nbreaks=nbreaks)
        
        # ROC curve
        plot(roc_pce$fpr, roc_pce$tpr,
             type = "l", las = 1, cex.lab = 1.5,
             xlab = "False Positive Rate", ylab = "True Positive Rate",
             col = "navy", lwd = 1.5
        )
        mtext(legend_letters[gender], side = 2, line = 2.5, at=1, cex=2, las=1)
        # lines(roc_rf$fpr, roc_rf$tpr, type="l", las=1, cex.lab=1.5,
        #       col="forestgreen", lwd=1.5)
        lines(roc_lasso$fpr, roc_lasso$tpr,
              type = "l", las = 1, cex.lab = 1.5,
              col = "tomato", lwd = 1.5, lty = 2
        )
        abline(0, 1, lty = 3)
        legend("bottomright",
               lty = c(1, 2, 1), col = c("navy", "tomato", "forestgreen"), lwd = 1.5, bty="n", cex=1.3,
               legend = c(
                 paste0("PCE: AUC = ", formatC(roc_pce$auc, format = "f", digits = 2)),
                 paste0("LASSO Stability Selection: AUC = ", formatC(roc_lasso$auc, format = "f", digits = 2))
               )
        )
      }
      dev.off()
    }
  }
}
