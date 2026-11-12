# 02 - Taille d'échantillon
# But : les effectifs sont imposés par des données existantes. On calcule donc
# ce qu'ils permettent de détecter (partie A) ou de modéliser (partie B), sans
# jamais utiliser l'effet du biomarqueur sur le critère de jugement
# (pas de « puissance a posteriori »).
suppressPackageStartupMessages({library(pmsampsize); library(dplyr)})
source("R/fonctions.R")
a <- readRDS("results/derive/flchain.rds")
b <- readRDS("results/derive/mgus2.rds")

# Partie A : formule de Schoenfeld (Hsieh et Lavori) pour une covariable continue.
# Il faut l'écart-type du marqueur (sigma) et la part de sa variance expliquée
# par les covariables d'ajustement (R²) : plus le marqueur est corrélé à l'âge
# ou à la créatinine, moins il apporte d'information propre.
cc <- filter(a, !is.na(lcreat))
sigma_a <- sd(cc$log2flc)
r2_a <- summary(lm(log2flc ~ splines::ns(age, 3) + sex + splines::ns(lcreat, 2),
                   data = cc))$r.squared
d_a <- sum(a$death)
res_a <- tibble(
  sigma_log2flc = sigma_a, r2_covariables = r2_a, deces = d_a,
  deces_requis_hr_1_20 = ceiling(schoenfeld_events(1.2, sigma_a, r2_a)),
  deces_requis_hr_1_10 = ceiling(schoenfeld_events(1.1, sigma_a, r2_a)),
  hr_min_detectable = schoenfeld_hr_min(d_a, sigma_a, r2_a)
)

# Partie B : événements par paramètre et critères de Riley (pmsampsize).
# Le R² de Cox-Snell attendu n'est pas connu : trois scénarios sont testés.
cc_b <- filter(b, !is.na(mspike), !is.na(hgb), !is.na(creat))
e_prog <- sum(cc_b$evt == "progression")
pers_ans <- sum(cc_b$ans)
parametres <- 5
scenarios <- lapply(c(0.03, 0.05, 0.10), function(r2) {
  s <- pmsampsize(type = "s", csrsquared = r2, parameters = parametres,
                  rate = e_prog / pers_ans, timepoint = 10,
                  meanfup = pers_ans / nrow(cc_b))
  tibble(r2cs_hypothese = r2, n_requis = s$sample_size,
         evenements_requis = s$events)
}) |> bind_rows()
res_b <- tibble(n = nrow(cc_b), progressions = e_prog, parametres = parametres,
                epv = e_prog / parametres,
                taux_pour_100_pa = 100 * e_prog / pers_ans,
                suivi_moyen_ans = pers_ans / nrow(cc_b))
write.csv(res_a, "results/tables/A_taille_echantillon.csv", row.names = FALSE)
write.csv(res_b, "results/tables/B_taille_echantillon.csv", row.names = FALSE)
write.csv(scenarios, "results/tables/B_riley.csv", row.names = FALSE)
print(res_a); print(res_b); print(scenarios)
