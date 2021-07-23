rm(list = ls())

library(openxlsx)


## Proportions of missing values

# Read supplementary tables by Matt
for (gender in c("m", "f")){
  supp=read.csv(paste0("Results/Supplementary_material_MW/missing_data_",gender,".csv"))
  supp=supp[sort.list(supp$Description),]
  supp=supp[,c(2,3,5)]
  supp[,1]=formatC(supp[,1], format="f", digits=2)
  supp[,2]=formatC(supp[,2], format="f", digits=2)
  assign(paste0("supp_",gender), supp)
}

# No missing for alanine aminotransferase in Women
supp_f=rbind(c("0.00","0.00",supp_m$Description[1]),supp_f)

# Merging results in men and women
print(all(supp_m$Description==supp_f$Description))
supp=cbind(Description=supp_m$Description, supp_m[,1:2], supp_f[,1:2])

# Saving the table
write.xlsx(supp, "Tables/Supp_table_1.xlsx", col.names=TRUE, row.names=FALSE, overwrite = TRUE)

