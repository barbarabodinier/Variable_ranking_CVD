reclassification <- function(data, cOutcome, predrisk1, predrisk2, cutoff) {
  c1 <- cut(predrisk1,
    breaks = cutoff, include.lowest = TRUE,
    right = FALSE
  )
  c2 <- cut(predrisk2,
    breaks = cutoff, include.lowest = TRUE,
    right = FALSE
  )
  tabReclas <- table(`Initial Model` = c1, `Updated Model` = c2)
  cat(" _________________________________________\n")
  cat(" \n     Reclassification table    \n")
  cat(" _________________________________________\n")
  ta <- table(c1, c2, data[, cOutcome])
  cat("\n Outcome: absent \n  \n")
  TabAbs <- ta[, , 1]
  tab1 <- cbind(TabAbs, ` % reclassified` = round((rowSums(TabAbs) -
    diag(TabAbs)) / rowSums(TabAbs), 2) * 100)
  names(dimnames(tab1)) <- c("Initial Model", "Updated Model")
  print(tab1)
  cat("\n \n Outcome: present \n  \n")
  TabPre <- ta[, , 2]
  tab2 <- cbind(TabPre, ` % reclassified` = round((rowSums(TabPre) -
    diag(TabPre)) / rowSums(TabPre), 2) * 100)
  names(dimnames(tab2)) <- c("Initial Model", "Updated Model")
  print(tab2)
  cat("\n \n Combined Data \n  \n")
  Tab <- tabReclas
  tab <- cbind(Tab, ` % reclassified` = round((rowSums(Tab) -
    diag(Tab)) / rowSums(Tab), 2) * 100)
  names(dimnames(tab)) <- c("Initial Model", "Updated Model")
  print(tab)
  cat(" _________________________________________\n")
  c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
  c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
  x <- improveProb(
    x1 = as.numeric(c11) * (1 / (length(levels(c11)))),
    x2 = as.numeric(c22) * (1 / (length(levels(c22)))), y = data[
      ,
      cOutcome
    ]
  )
  y <- improveProb(x1 = predrisk1, x2 = predrisk2, y = data[
    ,
    cOutcome
  ])
  cat(
    "\n NRI(Categorical) [95% CI]:", round(x$nri, 4), "[",
    round(x$nri - 1.96 * x$se.nri, 4), "-", round(x$nri +
      1.96 * x$se.nri, 4), "]", "; p-value:", round(2 *
      pnorm(-abs(x$z.nri)), 5), "\n"
  )
  cat(" NRI(Continuous) [95% CI]:", round(y$nri, 4), "[", round(y$nri -
    1.96 * y$se.nri, 4), "-", round(
    y$nri + 1.96 * y$se.nri,
    4
  ), "]", "; p-value:", round(
    2 * pnorm(-abs(y$z.nri)),
    5
  ), "\n")
  cat(" IDI [95% CI]:", round(y$idi, 4), "[", round(y$idi -
    1.96 * y$se.idi, 4), "-", round(
    y$idi + 1.96 * y$se.idi,
    4
  ), "]", "; p-value:", round(
    2 * pnorm(-abs(y$z.idi)),
    5
  ), "\n")
  return(list(cat_nri = x$nri, cat_nri_se = x$se.nri, cont_nri = y$nri, cont_nri_se = y$se.nri, idi = y$idi, idi_se = y$se.idi))
}
