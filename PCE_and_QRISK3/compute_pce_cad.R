rm(list=ls())
setwd("~/Dropbox/UK_Biobank/Data/GRS/cad_pce_1million/")
# setwd("~/Dropbox/UK_Biobank/Data/GRS/cvd_pce_1million/")

library(parallel)
library(survival)
library(splines)
source('~/Dropbox/UK_Biobank/Scripts/PCE_algorithm.R')

covars=readRDS("Data/covariates.rds")
covars=data.frame(covars)
rownames(covars)=covars$eid_236
# covars$smoking = covars$smoking + 1
# covars=covars[covars$cholesterol_med==0,]


### Removing missing values

# Removing individuals with missing sbp
table(is.na(covars$SBP.mean))
covars=covars[-which(is.na(covars$SBP.mean)),]
covars$untreated_sbp=ifelse(covars$cholesterol_med==0, yes=covars$SBP.mean, no=0)
covars$treated_sbp=ifelse(covars$cholesterol_med==1, yes=covars$SBP.mean, no=0)

# Removing individuals with missing eth 
table(is.na(covars$ethnicity))
# covars=covars[-which(is.na(covars$ethnicity)),]

# Removing individuals with missing cholesterol
table(is.na(covars$Cholesterol))
covars=covars[-which(is.na(covars$Cholesterol)),]

# Removing individuals with missing HDL cholesterol
table(is.na(covars$HDL_cholesterol))
covars=covars[-which(is.na(covars$HDL_cholesterol)),]

# Removing prevalent CAD/CVD and their matched controls
eid_cad=read.table("Data/eid_validation.txt")
eid_cvd=read.table("../cvd_1million/Data/eid_validation.txt")
eid=unique(c(eid_cad[,1], eid_cvd[,1]))
covars=covars[which(!covars$eid_236%in%eid),]

# # Retrieving info on "prefer not to say" smokers
# table(covars$smoking_status.0.0==-3)
# covars$smoking_status.0.0[covars$smoking_status.0.0==-3]=NA
# table(is.na(covars$smoking_status.0.0))
# covars$smoker=ifelse(covars$smoking_status.0.0=="2", yes=1, no=0)
# table(is.na(covars$smoker))
covars$smoker=ifelse(covars$smoking>0, yes=1, no=0)
table(covars$smoker)

# covars_qrisk=data.frame(readRDS("../cad_1million/Data/covariates.rds"))
# rownames(covars_qrisk)=covars_qrisk$eid_236
# ids=intersect(rownames(covars), rownames(covars_qrisk))
# # table(covars[ids,"smoking_status.0.0"], covars_qrisk[ids,"smoking"])
# 
# table(covars[ids,"any_CAD"], covars_qrisk[ids,"any_CAD"])
# table(covars[ids,"incident_CAD"], covars_qrisk[ids,"incident_CAD"])
# table(covars[ids,"prevalent_CAD"], covars_qrisk[ids,"prevalent_CAD"])

# eids_qrisk=

# eid_missing=covars[is.na(covars$smoker), "eid_236"]
# table(eid_missing%in%rownames(covars_qrisk))
# eid_missing=as.character(eid_missing[eid_missing%in%rownames(covars_qrisk)])
# covars[eid_missing, "smoker"]=ifelse(covars_qrisk[eid_missing, "smoking"]>0, yes=1, no=0)

# table(covars$smoker)
table(is.na(covars$smoker))
# covars=covars[-which(is.na(covars$smoker)),] # N=1,778
# covars$smoker=ifelse(covars$smoking_status.0.0=="2", yes=1, no=0)
# table(is.na(covars$smoker))


### Compute PCE

covars_females=covars[which(covars$sex.0.0==0),]
covars_males=covars[which(covars$sex.0.0==1),]

pce_females=ComputeFemalePCE(covars_females)
colnames(pce_females) = c('log_HR','Score')
rownames(pce_females) = covars_females$eid
ids=sample(nrow(pce_females), size=1000)
# par(mfrow=c(3,3))
# for(S0 in seq(0.7,0.75,by=0.01)){
# # S0=0.85
# plot(pce_females[ids,2], S0^exp(pce_females[ids,1]))
# abline(0,1,lty=2)
# }

pce_males=ComputeFemalePCE(covars_males)
colnames(pce_males) = c('log_HR','Score')
rownames(pce_males) = covars_males$eid


#combine male and females QRISK
pce_both = data.frame(rbind(pce_males,pce_females))
pce_both$eid = rownames(pce_both)
pce_both = merge(pce_both, covars, by = 'eid')

mytable_pce=pce_both[,2:4]
rownames(mytable_pce)=mytable_pce[,3]
mytable_pce=mytable_pce[,-3]
apply(mytable_pce, 2, FUN=function(x){sum(is.na(x))})

saveRDS(covars, "Data/covars_prospective.rds")
saveRDS(mytable_pce, "Data/predicted_pce_adj.rds")


# ### Compare with QRISK dataset
# dim(covars)
# covars_test_qrisk=readRDS("../cad_1million/Data/covars_test.rds")
# dim(covars_test_qrisk)
# eids_test_qrisk=rownames(covars_test_qrisk)
# eids_test_pce=covars$eid_236
# 
# table(eids_test_qrisk%in%eids_test_pce)
# table(covars_test_qrisk$case, eids_test_qrisk%in%eids_test_pce)
# 
# eids_lost=eids_test_qrisk[eids_test_qrisk%in%eids_test_pce]
# covars_qrisk_original=data.frame(readRDS("../cad_1million/Data/covariates.rds"))
# rownames(covars_qrisk_original)=covars_qrisk_original$eid_236
# View(covars_qrisk_original[eids_lost,])
# 
# write.xlsx(covars_qrisk_original[eids_lost, c("age_attending.0.0", "sex.0.0", "ethnicity",
#                                               "SBP.mean", "HTN_med", "diabetes", "smoking",
#                                               "Cholesterol", "HDL_cholesterol")], 
#            "Data/characteristics_lost_ind.xlsx", col.names=TRUE, row.names=TRUE)
# write.table(covars_qrisk_original[eids_lost, c("age_attending.0.0", "sex.0.0", "ethnicity",
#                                               "SBP.mean", "HTN_med", "diabetes", "smoking",
#                                               "Cholesterol", "HDL_cholesterol")], 
#            "Data/characteristics_lost_ind.txt", col.names=TRUE, row.names=TRUE)
# 
# 
# View(covars_qrisk_original[eids_lost, c("age_attending.0.0", "sex.0.0", "ethnicity",
#                                    "SBP.mean", "HTN_med", "diabetes", "smoking",
#                                    "Cholesterol", "HDL_cholesterol")])
# 



