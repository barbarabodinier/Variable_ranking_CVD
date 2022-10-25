rm(list = ls())
# setwd("~/Dropbox/UK_Biobank/Data/GRS/cvd_genotyped/")
setwd("~/Dropbox/UK_Biobank/Data/GRS/cvd_1m/")

library(parallel)
library(survival)
library(survminer)
library(splines)
source("~/Dropbox/UK_Biobank/Scripts/QRISK_algorithm.R")

covars <- readRDS("Data/covariates.rds")
covars$smoking <- covars$smoking + 1

covars$ethnicity <- with(covars, ifelse(ethnicity == "White/missing", 1,
  ifelse(ethnicity == "Indian", 2,
    ifelse(ethnicity == "Pakistani", 3,
      ifelse(ethnicity == "Bangladeshi", 4,
        ifelse(ethnicity == "Other Asian", 5,
          ifelse(ethnicity == "Caribbean", 6,
            ifelse(ethnicity == "African", 7,
              ifelse(ethnicity == "Chinese", 8, 9)
            )
          )
        )
      )
    )
  )
))
# covars$T1DM = ifelse(covars$diabetes == 1, 1, 0)
# covars$T2DM = ifelse(covars$diabetes == 2, 1, 0)
covars$T1DM <- ifelse(covars$diabetes_hba1c == 1, 1, 0)
covars$T2DM <- ifelse(covars$diabetes_hba1c == 2, 1, 0)

# subset to remove prevalent CVD and their matched controls
incident_controls <- subset(covars, prevalent_CVD == 0 & prev_CVD_controls == 0)

# subset by sex
females <- subset(incident_controls, sex.0.0 == 0)
males <- subset(incident_controls, sex.0.0 == 1)

QRISK_female <- t(sapply(1:nrow(females), function(i) {
  with(
    females,
    QRiskFemale(
      age = age_attending.0.0[i],
      b_AF = AF[i],
      b_atypicalantipsy = atypicals[i],
      b_corticosteroids = steroid_use[i],
      b_migraine = migraine[i],
      b_ra = RA[i],
      b_renal = CKD[i],
      b_semi = mental[i],
      b_sle = SLE[i],
      b_treatedhyp = HTN_med[i],
      b_type1 = T1DM[i],
      b_type2 = T2DM[i],
      bmi = bmi1.0.0[i],
      ethrisk = ethnicity[i],
      fh_cvd = first_degree_rel_CHD[i],
      rati = adj_chol_HDL_ratio[i],
      sbp = SBP.mean[i],
      sbps5 = SBP.sd[i],
      smoke_cat = smoking[i],
      town = deprivation.0.0[i]
    )
  )
}))
colnames(QRISK_female) <- c("log_HR", "Score")
rownames(QRISK_female) <- females$eid

# Male QRISK score / log(HR)
QRISK_male <- t(sapply(1:nrow(males), function(i) {
  with(
    males,
    QRiskMale(
      age = age_attending.0.0[i],
      b_AF = AF[i],
      b_atypicalantipsy = atypicals[i],
      b_corticosteroids = steroid_use[i],
      b_impotence2 = ED[i],
      b_migraine = migraine[i],
      b_ra = RA[i],
      b_renal = CKD[i],
      b_semi = mental[i],
      b_sle = SLE[i],
      b_treatedhyp = HTN_med[i],
      b_type1 = T1DM[i],
      b_type2 = T2DM[i],
      bmi = bmi1.0.0[i],
      ethrisk = ethnicity[i],
      fh_cvd = first_degree_rel_CHD[i],
      rati = adj_chol_HDL_ratio[i],
      sbp = SBP.mean[i],
      sbps5 = SBP.sd[i],
      smoke_cat = smoking[i],
      town = deprivation.0.0[i]
    )
  )
}))
colnames(QRISK_male) <- c("log_HR", "Score")
rownames(QRISK_male) <- males$eid

# combine male and females QRISK
QRISK_both <- data.frame(rbind(QRISK_male, QRISK_female))
QRISK_both$eid <- rownames(QRISK_both)
QRISK_both <- merge(QRISK_both, incident_controls, by = "eid")

mytable_qrisk <- QRISK_both[, 2:4]
rownames(mytable_qrisk) <- mytable_qrisk[, 3]
mytable_qrisk <- mytable_qrisk[, -3]
apply(mytable_qrisk, 2, FUN = function(x) {
  sum(is.na(x))
})

saveRDS(mytable_qrisk, "Data/predicted_qrisk_adj.rds")
