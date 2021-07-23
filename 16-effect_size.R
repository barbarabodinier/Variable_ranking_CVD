rm(list=ls())

library(basicPlotteR)

data_input="updated"
outcome="cvd"
model_id=1

gender="male"
mycox=readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
male=sort(ifelse(exp(coef(mycox$cox_model))<1, yes=1/exp(coef(mycox$cox_model)), no=exp(coef(mycox$cox_model))), decreasing=TRUE)

gender="female"
mycox=readRDS(paste0("Results/Cox_models/cox_model_lasso_stable_", outcome, "_m", model_id, "_", data_input, "_", gender, ".rds"))
female=sort(ifelse(exp(coef(mycox$cox_model))<1, yes=1/exp(coef(mycox$cox_model)), no=exp(coef(mycox$cox_model))), decreasing=TRUE)

ids=intersect(names(male), names(female))
male=male[ids]
female=female[ids]

{pdf("Figures/Strength_effect_m1.pdf")
  plot(male, female, pch=19, col="navy", las=1,
       xlab="Strengh of the effect (Men)",
       ylab="Strengh of the effect (Women)")
  abline(0,1,lty=2)
  addTextLabels(male,
                female,
                keepLabelsInside = FALSE,
                col.label = "navy",
                col.line = "navy",
                col.background = lighten("navy", amount = 0.8),
                labels = names(male)
  )
  dev.off()}
