rm(list = ls())

outcome_name <- "CVD"

# Loading main dataset
mydata <- readRDS(paste0("Data/", outcome_name, "_imputed_MW_full.rds"))

# Loading updated outcome data
outcome <- readRDS(paste0("Data/outcome_", outcome_name, ".rds"))

# Merging the datasets
eids <- intersect(rownames(mydata), rownames(outcome))
outcome <- outcome[eids, ]
mydata <- mydata[eids, ]
print(all(rownames(mydata) == rownames(outcome)))

table(mydata$case, outcome[rownames(mydata), "prevalent_case"], useNA = "ifany")
table(mydata$case, outcome[rownames(mydata), "incident_case"], useNA = "ifany")
table(mydata$case, outcome[rownames(mydata), "case"], useNA = "ifany")

# Excluding prevalent cases
eids <- eids[outcome$prevalent_case == 0]
outcome <- outcome[eids, ]
mydata <- mydata[eids, ]
print(all(rownames(mydata) == rownames(outcome)))

table(mydata$case, outcome[rownames(mydata), "prevalent_case"], useNA = "ifany")
table(mydata$case, outcome[rownames(mydata), "incident_case"], useNA = "ifany")
table(mydata$case, outcome[rownames(mydata), "case"], useNA = "ifany")

# Comparisons between previous and updated times to diagnosis
eids <- rownames(mydata)[mydata$case == 1]
plot(mydata[eids, "time"], outcome[eids, "time_to_diagnosis"]) # due to date of death

# Updating case/control status
print(all(rownames(mydata) == rownames(outcome)))
mydata$case <- outcome$incident_case

# Updating survival time
censorship_date <- as.Date("2021-04-07")
mydata$time <- outcome$time_to_diagnosis # time for incident cases
mydata$time[is.na(mydata$time)] <- outcome[is.na(mydata$time), "date_death"] - outcome[is.na(mydata$time), "date_recr"] # time for dead controls
mydata$time[is.na(mydata$time)] <- censorship_date - outcome[is.na(mydata$time), "date_recr"] # time for living controls
eids <- rownames(mydata)[mydata$case == 1]
plot(mydata[eids, "time"], outcome[eids, "time_to_diagnosis"]) # checking survival time for cases
median(mydata$time) / 365.25

# Saving updated dataset
saveRDS(mydata, paste0("Data/", outcome_name, "_imputed_MW_updated.rds"))

# Using censorship at 10 years
mydata <- readRDS(paste0("Data/", outcome_name, "_imputed_MW_updated.rds"))
mydata$case[which(outcome$time_to_diagnosis > 10 * 365.25)] <- 0
mydata$time[which(mydata$time > 10 * 365.25)] <- round(10 * 365.25)
saveRDS(mydata, paste0("Data/", outcome_name, "_imputed_MW_censored.rds"))

# Using censorship in 2017
censorship_date <- as.Date("2017-03-31")
mydata <- readRDS(paste0("Data/", outcome_name, "_imputed_MW_updated.rds"))
mydata$case[which(outcome$date_diagnosis > censorship_date)] <- 0
outcome[which(outcome$date_diagnosis > censorship_date), "time"] <- NA
mydata$time <- outcome$time_to_diagnosis
mydata$time[is.na(mydata$time)] <- censorship_date - outcome[is.na(mydata$time), "date_recr"]
saveRDS(mydata, paste0("Data/", outcome_name, "_imputed_MW_2017.rds"))
