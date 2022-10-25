ComputeFemalePCE <- function(input) {
  input$age <- input$age_attending.0.0
  input$black <- ifelse(input$ethnicity == "Black", yes = 1, no = 0)
  input$sysbp <- input$SBP.mean
  input$rxbp <- input$HTN_med
  input$dm <- input$diabetes
  input$cursmoke <- input$smoker
  input$totchol <- input$Cholesterol
  input$hdlc <- input$HDL_cholesterol

  loghr <- (
    0.106501 * as.numeric(input$age) +
      0.119399 * as.numeric(input$black) +
      0.000056 * (as.numeric(input$sysbp)^2) +
      0.017666 * as.numeric(input$sysbp) +
      0.731678 * as.numeric(input$rxbp) +
      0.943970 * as.numeric(input$dm) +
      1.009790 * as.numeric(input$cursmoke) +
      0.151318 * (as.numeric(input$totchol) / as.numeric(input$hdlc)) +
      -0.008580 * as.numeric(input$age) * as.numeric(input$black) +
      -0.003647 * as.numeric(input$sysbp) * as.numeric(input$rxbp) +
      0.006208 * as.numeric(input$sysbp) * as.numeric(input$black) +
      0.152968 * as.numeric(input$black) * as.numeric(input$rxbp) +
      -0.000153 * as.numeric(input$age) * as.numeric(input$sysbp) +
      0.115232 * as.numeric(input$black) * as.numeric(input$dm) +
      -0.092231 * as.numeric(input$black) * as.numeric(input$cursmoke) +
      0.070498 * as.numeric(input$black) * (as.numeric(input$totchol) / as.numeric(input$hdlc)) +
      -0.000173 * as.numeric(input$black) * as.numeric(input$sysbp) * as.numeric(input$rxbp) +
      -0.000094 * as.numeric(input$age) * as.numeric(input$sysbp) * as.numeric(input$black))

  female.risk <- 1.0 / (1.0 + exp(-(
    -11.938468 +
      0.106501 * as.numeric(input$age) +
      0.119399 * as.numeric(input$black) +
      0.000056 * (as.numeric(input$sysbp)^2) +
      0.017666 * as.numeric(input$sysbp) +
      0.731678 * as.numeric(input$rxbp) +
      0.943970 * as.numeric(input$dm) +
      1.009790 * as.numeric(input$cursmoke) +
      0.151318 * (as.numeric(input$totchol) / as.numeric(input$hdlc)) +
      -0.008580 * as.numeric(input$age) * as.numeric(input$black) +
      -0.003647 * as.numeric(input$sysbp) * as.numeric(input$rxbp) +
      0.006208 * as.numeric(input$sysbp) * as.numeric(input$black) +
      0.152968 * as.numeric(input$black) * as.numeric(input$rxbp) +
      -0.000153 * as.numeric(input$age) * as.numeric(input$sysbp) +
      0.115232 * as.numeric(input$black) * as.numeric(input$dm) +
      -0.092231 * as.numeric(input$black) * as.numeric(input$cursmoke) +
      0.070498 * as.numeric(input$black) * (as.numeric(input$totchol) / as.numeric(input$hdlc)) +
      -0.000173 * as.numeric(input$black) * as.numeric(input$sysbp) * as.numeric(input$rxbp) +
      -0.000094 * as.numeric(input$age) * as.numeric(input$sysbp) * as.numeric(input$black))))

  return(cbind(loghr, female.risk))
}


ComputeMalePCE <- function(input) {
  input$age <- input$age_attending.0.0
  input$black <- ifelse(input$ethnicity == "Black", yes = 1, no = 0)
  input$sysbp <- input$SBP.mean
  input$rxbp <- input$HTN_med
  input$dm <- input$diabetes
  input$cursmoke <- input$smoker
  input$totchol <- input$Cholesterol
  input$hdlc <- input$HDL_cholesterol

  loghr <- (
    0.064200 * as.numeric(input$age) +
      0.315986 * as.numeric(input$black) +
      -0.000061 * (as.numeric(input$sysbp)^2) +
      0.038950 * as.numeric(input$sysbp) +
      2.055533 * as.numeric(input$rxbp) +
      0.842209 * as.numeric(input$dm) +
      0.895589 * as.numeric(input$cursmoke) +
      0.193307 * (as.numeric(input$totchol) / as.numeric(input$hdlc)) +
      -0.014207 * as.numeric(input$sysbp) * as.numeric(input$rxbp) +
      0.011609 * as.numeric(input$sysbp) * as.numeric(input$black) +
      -0.119460 * as.numeric(input$rxbp) * as.numeric(input$black) +
      0.000025 * as.numeric(input$age) * as.numeric(input$sysbp) +
      -0.077214 * as.numeric(input$black) * as.numeric(input$dm) +
      -0.226771 * as.numeric(input$black) * as.numeric(input$cursmoke) +
      -0.117749 * (as.numeric(input$totchol) / as.numeric(input$hdlc)) * as.numeric(input$black) +
      0.004190 * as.numeric(input$black) * as.numeric(input$rxbp) * as.numeric(input$sysbp) +
      -0.000199 * as.numeric(input$black) * as.numeric(input$age) * as.numeric(input$sysbp))

  male.risk <- 1.0 / (1.0 + exp(-(
    -11.219734 +
      0.064200 * as.numeric(input$age) +
      0.315986 * as.numeric(input$black) +
      -0.000061 * (as.numeric(input$sysbp)^2) +
      0.038950 * as.numeric(input$sysbp) +
      2.055533 * as.numeric(input$rxbp) +
      0.842209 * as.numeric(input$dm) +
      0.895589 * as.numeric(input$cursmoke) +
      0.193307 * (as.numeric(input$totchol) / as.numeric(input$hdlc)) +
      -0.014207 * as.numeric(input$sysbp) * as.numeric(input$rxbp) +
      0.011609 * as.numeric(input$sysbp) * as.numeric(input$black) +
      -0.119460 * as.numeric(input$rxbp) * as.numeric(input$black) +
      0.000025 * as.numeric(input$age) * as.numeric(input$sysbp) +
      -0.077214 * as.numeric(input$black) * as.numeric(input$dm) +
      -0.226771 * as.numeric(input$black) * as.numeric(input$cursmoke) +
      -0.117749 * (as.numeric(input$totchol) / as.numeric(input$hdlc)) * as.numeric(input$black) +
      0.004190 * as.numeric(input$black) * as.numeric(input$rxbp) * as.numeric(input$sysbp) +
      -0.000199 * as.numeric(input$black) * as.numeric(input$age) * as.numeric(input$sysbp))))

  return(cbind(loghr, male.risk))
}
