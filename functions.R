RunLogistic <- function(covars, predictors, outcome, conf) {
  summary <- matrix(NA, ncol = 3, nrow = ncol(predictors))
  for (k in 1:ncol(predictors)) {
    # print(k)
    mydatauniv <- cbind(covars[, c(outcome, conf)], x = predictors[, k])
    mydatauniv <- na.exclude(mydatauniv)
    mydatauniv$outcome <- as.numeric(eval(parse(text = paste0("mydatauniv$", outcome))))
    mydatauniv$x <- scale(mydatauniv$x)
    f <- paste0("outcome~x+", paste(conf, collapse = "+"))
    mymodel <- glm(as.formula(f), data = mydatauniv, family = "binomial")
    f0 <- paste0("outcome~", paste(conf, collapse = "+"))
    mymodel0 <- glm(as.formula(f0), data = mydatauniv, family = "binomial")
    myanova <- anova(mymodel0, mymodel, test = "LRT")
    summary[k, ] <- c(nrow(mydatauniv), formatC(exp(coef(mymodel)["x"]), format = "f", digits = 2),
      pval = formatC(myanova$`Pr(>Chi)`[2], format = "e", digits = 2)
    )
  }
  colnames(summary) <- c("nobs", "coef", "pval")
  rownames(summary) <- colnames(predictors)
  summary <- as.data.frame(summary)
  return(summary)
}


GetConcordance <- function(coxph = NULL, survobject, newdata = NULL, S = NULL, times = 3653, digits = 3) {
  if (is.null(S)) {
    S0 <- summary(survfit(coxph), times = times, extend = TRUE)$surv
    mylp <- predict(coxph, newdata = newdata, type = "lp")
    S <- S0^exp(mylp)
  }
  c <- survConcordance(survobject ~ S)
  cs <- 1 - c$concordance
  cs_lower <- cs - c$std.err * 1.96
  cs_upper <- cs + c$std.err * 1.96
  print(paste0(
    formatC(cs, format = "f", digits = digits), " [",
    formatC(cs_lower, format = "f", digits = digits), "-",
    formatC(cs_upper, format = "f", digits = digits), "]"
  ))
}


DeltaCStat <- function(c1, c2) {
  mean1 <- as.numeric(gsub(" \\[.*", "", c1))
  se1 <- (mean1 - as.numeric(gsub("-.*", "", gsub(".*\\[", "", c1)))) / 1.96

  mean2 <- as.numeric(gsub(" \\[.*", "", c2))
  se2 <- (mean2 - as.numeric(gsub("-.*", "", gsub(".*\\[", "", c2)))) / 1.96

  return(paste0(
    formatC(mean2 - mean1, format = "f", digits = 3),
    " [", formatC(mean2 - mean1 - 1.96 * sqrt(se2^2 + se1^2), format = "f", digits = 3), "-",
    formatC(mean2 - mean1 + 1.96 * sqrt(se2^2 + se1^2), format = "f", digits = 3), "]"
  ))
}


GetSurvivalProbabilities <- function(coxph, survobject, newdata, times = 3653) {
  S0 <- summary(survfit(coxph), times = times, extend = TRUE)$surv
  mylp <- predict(coxph, newdata = newdata, type = "lp")
  S <- S0^exp(mylp)
  return(S)
}


CategoricalToFactor <- function(x, categorical) {
  mycategorical <- intersect(colnames(x), categorical)
  for (k in which(colnames(x) %in% mycategorical)) {
    x[, k] <- as.factor(x[, k])
  }
  return(x)
}


GetCapacity <- function(proba_of_event, case, N = 1000) {
  # Extracting all visited probabilities of event
  proba_list <- sort(unique(proba_of_event))

  # Reducing the numbers for efficiency
  proba_list <- proba_list[unique(round(seq(1, length(proba_list), length.out = N)))]

  # Initialisation
  myproportion <- mytpr <- NULL
  pb <- utils::txtProgressBar(style = 3)

  for (k in 1:length(proba_list)) {
    # Extracting at-risk participants for a given threshold
    p <- proba_list[k]
    tmppred <- ifelse(proba_of_event >= p, yes = 1, no = 0)
    tmppred <- factor(tmppred, levels = c(0, 1))

    # Comparing at-risk with actual incident cases
    cont <- table(tmppred, case)
    myproportion <- c(myproportion, (cont["1", "0"] + cont["1", "1"]) / sum(cont))
    mytpr <- c(mytpr, cont["1", "1"] / (cont["1", "1"] + cont["0", "1"]))

    # Increment loading bar
    utils::setTxtProgressBar(pb, k / length(proba_list))
  }

  # Including a probability of 0
  myproportion <- c(myproportion, 0)
  mytpr <- c(mytpr, 0)

  return(list(proportion = myproportion, tpr = mytpr))
}


# RunCoxSelected=function(mydata, selected, verbose=TRUE){
#   mydata_pred=mydata[,selected,drop=FALSE]
#   x=CategoricalToFactor(mydata_pred, categorical)
#   x_test=CategoricalToFactor(mydata_test[,colnames(mydata_pred),drop=FALSE], categorical)
#   survobject=Surv(mydata$time, mydata$case)
#   myformula=paste(colnames(mydata_pred), collapse="+")
#   if (verbose){
#     print(myformula)
#   }
#   mycox=coxph(as.formula(paste0("survobject~", myformula)), data=x)
#
#   concordance(mycox)$concordance
#   c_selected_train=GetConcordance(mycox, survobject, x)
#
#   survobject_test=Surv(mydata_test$time, mydata_test$case)
#   c_selected_test=GetConcordance(mycox, survobject_test, x_test)
#
#   return(c(c_selected_train, c_selected_test))
# }


RunCoxSelected <- function(mydata, mydata_test, selected, verbose = TRUE, penalisation = FALSE, seed = 1, alpha = 0, digits = 3) {
  # mydata=data.frame(scale(mydata))
  # mydata_test=data.frame(scale(mydata_test))
  x <- mydata
  mydata_pred <- mydata[, selected, drop = FALSE]
  x_test <- mydata_test
  # x=CategoricalToFactor(mydata_pred, categorical)
  # x_test=CategoricalToFactor(mydata_test[,colnames(mydata_pred),drop=FALSE], categorical)
  survobject <- Surv(mydata$time, mydata$case)

  if (!penalisation) {
    myformula <- paste(colnames(mydata_pred), collapse = "+")
    if (verbose) {
      print(myformula)
    }
    mycox <- coxph(as.formula(paste0("survobject~", myformula)), data = x)
  } else {
    y <- as.matrix(mydata[, c("time", "case")])
    colnames(y) <- c("time", "status")
    set.seed(seed)
    mypenalisedcox <- cv.glmnet(data.matrix(x), y, family = "cox", alpha = alpha)
    mycalibratedpenalisedcox <- glmnet(data.matrix(x), y, family = "cox", lambda = mypenalisedcox$lambda.1se, alpha = alpha)

    beta <- as.vector(mycalibratedpenalisedcox$beta)
    names(beta) <- rownames(mycalibratedpenalisedcox$beta)
    nonzero <- names(beta[beta != 0])

    mycox <- coxph(as.formula(paste0("survobject~", paste(nonzero, collapse = "+"))),
      data = mydata,
      init = beta[nonzero], iter = 0
    ) # setting the coefficients to the ones estimated in penalised Cox
  }

  c_selected_train <- GetConcordance(mycox, survobject, x, digits = digits)
  S_train <- GetSurvivalProbabilities(mycox, survobject, x)
  survobject_test <- Surv(mydata_test$time, mydata_test$case)
  c_selected_test <- GetConcordance(mycox, survobject_test, x_test, digits = digits)
  S_test <- GetSurvivalProbabilities(mycox, survobject_test, x_test)

  if (!penalisation) {
    return(list(
      c_training = c_selected_train, c_test = c_selected_test, cox_model = mycox,
      S_training = S_train, S_test = S_test
    ))
  } else {
    return(list(
      c_training = c_selected_train, c_test = c_selected_test, cox_model = mycox, penalised_cox = mypenalisedcox,
      S_training = S_train, S_test = S_test
    ))
  }
}


CreateFigComp <- function(tmp_Q_sex_rec, tmp_Q_GRS_sex_rec, name, myrange = c(0, 70), myrange2 = c(-50, 50)) {{
  postscript(paste0("Figures/Survival_diff_", name, ".eps"),
    width = 7, height = 7
  )
  mat <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
  layout(mat, height = c(1, 4), width = c(4, 1))
  par(mar = c(0, 5, 1, 1), lend = 2)
  h <- hist(100 * (1 - tmp_Q_sex_rec), plot = FALSE, breaks = 100)
  plot(NA,
    type = "n", axes = FALSE, yaxt = "n",
    ylab = expression(Counts ~ (x10^5)), xlab = "", main = NA,
    xlim = myrange,
    ylim = c(0, max(h$counts))
  )
  # xlim=range(h$breaks))
  axis(2, las = 1, at = axTicks(2), labels = axTicks(2) / 10000)
  for (k in 1:length(h$counts)) {
    polygon(
      x = c(h$breaks[c(k, k + 1)], h$breaks[c(k + 1, k, k)]) + c(0.3, -0.3, -0.3, 0.3, 0.3),
      y = c(0, 0, h$counts[k], h$counts[k], 0), col = "grey30", border = NA
    )
  }
  plot.new()
  par(mar = c(5, 5, 1, 1))

  tmp_Q_sex_rec1 <- tmp_Q_sex_rec
  tmp_Q_GRS_sex_rec1 <- tmp_Q_GRS_sex_rec
  set.seed(1)
  # ids=sample(length(tmp_Q_sex_rec), size=0.01*length(tmp_Q_sex_rec))
  ids <- 1:length(tmp_Q_sex_rec)
  tmp_Q_sex_rec1 <- tmp_Q_sex_rec1[ids]
  tmp_Q_GRS_sex_rec1 <- tmp_Q_GRS_sex_rec1[ids]
  plot(100 * (1 - tmp_Q_sex_rec1), 100 * (1 - tmp_Q_GRS_sex_rec1) - 100 * (1 - tmp_Q_sex_rec1),
    pch = 19, cex = 0.1,
    # col=adjustcolor("grey30", alpha.f = 0.3),
    col = "grey30",
    xlim = myrange, ylim = myrange2,
    las = 1, xlab = "PCE 10-year risk (%)",
    ylab = "Difference between 10-year risk by PCE and selected (%)",
    panel.first = c(
      abline(h = seq(-1, 1, by = 0.1) * 100, col = "grey", lty = 3),
      abline(v = seq(-1, 1, by = 0.1) * 100, col = "grey", lty = 3)
    )
  )

  par(mar = c(5, 0, 1, 1), lend = 2)
  h <- hist((100 * (1 - tmp_Q_GRS_sex_rec) - 100 * (1 - tmp_Q_sex_rec)), plot = FALSE, breaks = 100)
  plot(NA,
    type = "n", axes = FALSE, yaxt = "n",
    xlab = expression(Counts ~ (x10^5)), ylab = NA, main = NA,
    ylim = myrange2,
    xlim = c(0, max(h$counts))
  )
  # ylim=range(h$breaks))
  axis(1, las = 1, at = axTicks(1), labels = axTicks(1) / 10000)
  for (k in 1:length(h$counts)) {
    polygon(
      y = c(h$breaks[c(k, k + 1)], h$breaks[c(k + 1, k, k)]) + c(0.3, -0.3, -0.3, 0.3, 0.3),
      x = c(0, 0, h$counts[k], h$counts[k], 0), col = "grey30", border = NA
    )
  }
  dev.off()
}

{
  png(paste0("Figures/Survival_diff_", name, ".png"), width = 7, height = 7, unit = "in", res = 500, pointsize = 12)
  mat <- matrix(c(1, 2, 3, 4), ncol = 2, byrow = TRUE)
  layout(mat, height = c(1, 4), width = c(4, 1))
  par(mar = c(0, 5, 1, 1), lend = 2)
  h <- hist(100 * (1 - tmp_Q_sex_rec), plot = FALSE, breaks = 100)
  plot(NA,
    type = "n", axes = FALSE, yaxt = "n",
    ylab = expression(Counts ~ (x10^5)), xlab = "", main = NA,
    xlim = myrange,
    ylim = c(0, max(h$counts))
  )
  # xlim=range(h$breaks))
  axis(2, las = 1, at = axTicks(2), labels = axTicks(2) / 10000)
  for (k in 1:length(h$counts)) {
    polygon(
      x = c(h$breaks[c(k, k + 1)], h$breaks[c(k + 1, k, k)]) + c(0.3, -0.3, -0.3, 0.3, 0.3),
      y = c(0, 0, h$counts[k], h$counts[k], 0), col = "grey30", border = NA
    )
  }
  plot.new()
  par(mar = c(5, 5, 1, 1))

  tmp_Q_sex_rec1 <- tmp_Q_sex_rec
  tmp_Q_GRS_sex_rec1 <- tmp_Q_GRS_sex_rec
  set.seed(1)
  # ids=sample(length(tmp_Q_sex_rec), size=0.01*length(tmp_Q_sex_rec))
  ids <- 1:length(tmp_Q_sex_rec)
  tmp_Q_sex_rec1 <- tmp_Q_sex_rec1[ids]
  tmp_Q_GRS_sex_rec1 <- tmp_Q_GRS_sex_rec1[ids]
  plot(100 * (1 - tmp_Q_sex_rec1), 100 * (1 - tmp_Q_GRS_sex_rec1) - 100 * (1 - tmp_Q_sex_rec1),
    pch = 19, cex = 0.1,
    # col=adjustcolor("grey30", alpha.f = 0.3),
    col = "grey30",
    xlim = myrange, ylim = myrange2,
    las = 1, xlab = "PCE 10-year risk (%)",
    ylab = "Difference between 10-year risk by PCE and selected (%)",
    panel.first = c(
      abline(h = seq(-1, 1, by = 0.1) * 100, col = "grey", lty = 3),
      abline(v = seq(-1, 1, by = 0.1) * 100, col = "grey", lty = 3)
    )
  )

  par(mar = c(5, 0, 1, 1), lend = 2)
  h <- hist((100 * (1 - tmp_Q_GRS_sex_rec) - 100 * (1 - tmp_Q_sex_rec)), plot = FALSE, breaks = 100)
  plot(NA,
    type = "n", axes = FALSE, yaxt = "n",
    xlab = expression(Counts ~ (x10^5)), ylab = NA, main = NA,
    ylim = myrange2,
    xlim = c(0, max(h$counts))
  )
  # ylim=range(h$breaks))
  axis(1, las = 1, at = axTicks(1), labels = axTicks(1) / 10000)
  for (k in 1:length(h$counts)) {
    polygon(
      y = c(h$breaks[c(k, k + 1)], h$breaks[c(k + 1, k, k)]) + c(0.3, -0.3, -0.3, 0.3, 0.3),
      x = c(0, 0, h$counts[k], h$counts[k], 0), col = "grey30", border = NA
    )
  }
  dev.off()
}}


### Recalibration plots

GetObsCumInc <- function(survobject, linear.predictors) {
  ## KM deciles
  decilesb <- quantile(linear.predictors, probs = seq(0, 1, by = 0.1)) # decile boundaries
  deciles <- cut(linear.predictors,
    breaks = decilesb,
    include.lowest = TRUE, labels = seq(1:10)
  )

  survfit <- survfit(survobject ~ deciles)
  S_obs <- summary(survfit, times = 3653, extend = TRUE)[c("surv", "lower", "upper")]

  return(list(S_obs = S_obs, deciles = deciles))
}


CalibPlot <- function(S_pred, S_obs, name = "", xylim = 1) {
  # S_obs is a vector with the observed, lower and upper S(10)
  plotCI(
    x = 1 - S_pred, y = 1 - S_obs$surv,
    li = 1 - S_obs$upper, ui = 1 - S_obs$lower,
    pch = 19, las = 1, col = "navy", main = name,
    ylim = c(0, xylim), xlim = c(0, xylim), cex.lab = 1.5,
    xlab = "Cumulative incidence (predicted)", ylab = "Cumulative incidence (observed)"
  )
  y <- (1 - S_obs$surv)
  x <- (1 - S_pred)
  slopemodel <- lm(y ~ x)
  abline(0, 1, lty = 2, col = "grey30")
  return(formatC(coef(slopemodel)[2], format = "f", digits = 2))
}


GetGNDStat <- function(survobject, predicted_S, mydata) {
  ## KM deciles
  deciles <- as.numeric(cut2(predicted_S, g = 10))
  censt <- max(mydata$time) + 1
  gnd <- GND.calib(
    pred = predicted_S,
    tvar = mydata$time,
    out = mydata$case,
    cens.t = censt,
    groups = as.numeric(as.character(deciles)),
    adm.cens = censt
  )
  return(gnd)
}


ReformatAUC=function(auc, digits=2){
  return(paste0(formatC(auc[2], format="f", digits=digits),
                " [", 
                formatC(auc[1], format="f", digits=digits), 
                "-", 
                formatC(auc[3], format="f", digits=digits), 
                "]"))
}
