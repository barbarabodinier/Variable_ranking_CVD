setwd("../Results/Sensitivity")

args=commandArgs(trailingOnly=TRUE)
model_id=args[1]
gender=as.character(args[2])

selprop=selected=NULL
for (k in 1:100){
tmp=readRDS(paste0("stability_cvd_m",model_id,"_updated_",gender,"_",k,".rds"))
selprop=cbind(selprop, tmp[,1])
selected=cbind(selected, tmp[,2])
}

saveRDS(selprop, paste0("stability_proportions_cvd_m",model_id,"_updated_",gender,"_merged.rds"))
saveRDS(selected, paste0("stability_selected_cvd_m",model_id,"_updated_",gender,"_merged.rds"))

