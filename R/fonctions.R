# Fonctions partagées par les scripts d'analyse.

#' DFG estimé par l'équation CKD-EPI 2021 (sans coefficient ethnique).
#' @param scr créatinine sérique en mg/dL
#' @param age âge en années
#' @param sexe "F" ou "M"
#' @return DFG estimé en mL/min/1,73 m²
ckdepi2021 <- function(scr, age, sexe) {
  k <- ifelse(sexe == "F", 0.7, 0.9)
  a <- ifelse(sexe == "F", -0.241, -0.302)
  142 * pmin(scr / k, 1)^a * pmax(scr / k, 1)^(-1.200) *
    0.9938^age * ifelse(sexe == "F", 1.012, 1)
}

#' Nombre d'événements requis (Schoenfeld) pour une covariable continue.
#' @param hr rapport de risques pour une unité de la covariable
#' @param sigma écart-type de la covariable
#' @param r2 part de variance de la covariable expliquée par les autres
schoenfeld_events <- function(hr, sigma, r2 = 0, alpha = 0.05, puissance = 0.90) {
  (qnorm(1 - alpha / 2) + qnorm(puissance))^2 /
    (sigma^2 * log(hr)^2 * (1 - r2))
}

#' Plus petit HR détectable pour un nombre d'événements donné.
schoenfeld_hr_min <- function(evenements, sigma, r2 = 0, alpha = 0.05,
                              puissance = 0.90) {
  exp((qnorm(1 - alpha / 2) + qnorm(puissance)) /
        (sigma * sqrt(evenements * (1 - r2))))
}

#' Combinaison de Rubin pour un paramètre scalaire.
#' @param est estimations par imputation
#' @param var variances correspondantes
#' @return data.frame avec estimation, erreur type, ddl, p
rubin <- function(est, var) {
  m <- length(est)
  q <- mean(est); w <- mean(var); b <- stats::var(est)
  t <- w + (1 + 1 / m) * b
  r <- (1 + 1 / m) * b / w
  ddl <- (m - 1) * (1 + 1 / r)^2
  data.frame(est = q, se = sqrt(t), ddl = ddl,
             p = 2 * pt(-abs(q / sqrt(t)), ddl),
             fmi = (r + 2 / (ddl + 3)) / (r + 1))
}

#' Ajoute un IC à 95 % exponentié à une sortie de rubin().
hr_ic <- function(p) {
  q <- qt(.975, p$ddl)
  data.frame(hr = exp(p$est), bas = exp(p$est - q * p$se),
             haut = exp(p$est + q * p$se), p = p$p, fmi = p$fmi)
}

#' E-value pour un HR, avec un événement fréquent (VanderWeele & Ding, 2017).
#' Le HR est d'abord converti en risque relatif approché.
#' @param hr HR observé ; un HR < 1 est traité par son inverse
evalue_hr <- function(hr) {
  hr <- ifelse(hr < 1, 1 / hr, hr)
  rr <- (1 - 0.5^sqrt(hr)) / (1 - 0.5^sqrt(1 / hr))
  rr + sqrt(rr * (rr - 1))
}

#' Formate un rapport de risques et son IC.
fmt_hr <- function(hr, bas, haut, d = 2) {
  sprintf(paste0("%.", d, "f (%.", d, "f à %.", d, "f)"), hr, bas, haut) |>
    gsub(pattern = ".", replacement = ",", fixed = TRUE)
}

#' Formate une valeur de p.
fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0,001", gsub(".", ",", sprintf("%.3f", p), fixed = TRUE))
}
