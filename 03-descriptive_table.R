rm(list = ls())

library(openxlsx)

# Parameters
outcome <- "cvd"
data_input <- "updated"
genders <- c("male", "female")
model_id <- 1

# Variable annotation
dict <- read.xlsx("Data/Data_dictionary.xlsx")
rownames(dict) <- dict[, 1]
variable_cat <- dict[, 3]
names(variable_cat) <- dict[, 2]
nightingale <- names(which(variable_cat == "Nightingale"))

# Looping over models
for (model_id in c(1, 3)) {
  print(model_id)
  for (gender in genders) {
    print(gender)

    mytable_cont <- mytable_bin <- NULL
    for (status in c("0", "1")) {
      # Loading the data
      mydata <- data.frame(readRDS(paste0("Data/", toupper(outcome), "_imputed_MW_", data_input, ".rds")))

      # Extracting stratum
      ids_rows <- which((mydata$sex == ifelse(gender == "male", yes = 1, no = 0)) & (mydata$case == status))
      tmpdata <- mydata[ids_rows, ]

      # Extracting variables
      if (gender == "male") {
        ids_columns <- rownames(dict)[which(dict[, paste0("Model.", model_id, ".(Male)")] == "X")]
      } else {
        ids_columns <- rownames(dict)[which(dict[, paste0("Model.", model_id, ".(Female)")] == "X")]
      }
      tmpdata <- tmpdata[, ids_columns]

      # Using complete case
      tmpdata <- na.exclude(tmpdata)

      # Identifying binary variables
      binary <- ifelse(apply(tmpdata, 2, FUN = function(x) {
        length(unique(x))
      }) < 3, yes = 1, no = 0)

      # Identifying log-transformed variables
      logged <- ifelse(dict[ids_columns, 3] %in% c("Biochemistry", "Haematology", "Nightingale"), yes = 1, no = 0)

      # Continuous variables
      continuous_mean <- c(
        formatC(apply(tmpdata[, which((binary == 0) & (logged == 0))], 2, mean),
          format = "f", digits = 2
        ),
        formatC(apply(tmpdata[, which((binary == 0) & (logged == 1))], 2, FUN = function(x) {
          mean(exp(x))
        }),
        format = "f", digits = 2
        )
      )
      continuous_sd <- c(
        formatC(apply(tmpdata[, which((binary == 0) & (logged == 0))], 2, sd),
          format = "f", digits = 2
        ),
        formatC(apply(tmpdata[, which((binary == 0) & (logged == 1))], 2, FUN = function(x) {
          sd(exp(x))
        }),
        format = "f", digits = 2
        )
      )
      tmpcont <- c(
        formatC(nrow(tmpdata), format = "f", digits = 0, big.mark = ","),
        paste0(continuous_mean, " (", continuous_sd, ")")
      )

      # Binary variables
      count <- formatC(apply(tmpdata[, which(binary == 1)], 2, sum), format = "f", digits = 0, big.mark = ",")
      percentage <- formatC(apply(tmpdata[, which(binary == 1)], 2, sum) / nrow(tmpdata) * 100, format = "f", digits = 2)
      tmpbin <- c(
        paste0(count, " (", percentage, ")")
      )

      # Setting row names
      if (status == "0") {
        mytable_cont <- rbind(rep("", 2), cbind(dict[names(continuous_mean), c(3, 2)]))
        mytable_bin <- cbind(dict[names(count), c(3, 2)])
      }

      # Concatenating the columns
      mytable_cont <- cbind(mytable_cont, tmpcont)
      mytable_bin <- cbind(mytable_bin, tmpbin)
    }
    colnames(mytable_cont) <- colnames(mytable_bin) <- c("Category", "Variable", "Non-cases", "Cases")
    mytable <- rbind(mytable_cont, rep("", 4), mytable_bin)

    write.xlsx(as.data.frame(mytable), paste0("Tables/Descriptive_table_", outcome, "_m", model_id, "_", data_input, "_", gender, ".xlsx"),
      rowNames = FALSE, colNames = TRUE
    )
  }
}
