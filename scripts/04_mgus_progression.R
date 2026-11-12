# 04 - Partie B : pic monoclonal et progression de la MGUS
# But : estimer le risque de progression en traitant le décès sans progression
# comme un risque concurrent (et non comme une censure), mesurer l'effet du pic
# monoclonal avec deux modèles complémentaires (Cox par cause pour l'étiologie,
# Fine-Gray pour le risque absolu), puis sa valeur ajoutée pour la prédiction.
suppressPackageStartupMessages({
  library(survival); library(dplyr); library(tidyr); library(ggplot2)
  library(cmprsk); library(riskRegression); library(prodlim); library(purrr); library(splines)
})
source("R/fonctions.R"); source("R/theme.R")
set.seed(20261104)
b <- readRDS("results/derive/mgus2.rds")
tab <- function(x, f) write.csv(x, file.path("results/tables", f), row.names = FALSE)

# 1. Incidence cumulée --------------------------------------------------------
# survfit() avec un facteur à plusieurs états et id = Aalen-Johansen. La courbe
# 1 - Kaplan-Meier, qui censure les décès, est gardée pour montrer le biais.
aj <- survfit(Surv(ans, evt) ~ 1, data = b, id = id)
km <- survfit(Surv(ans, evt == "progression") ~ 1, data = b)
etats <- aj$states
idx_p <- match("progression", etats); idx_d <- match("deces", etats)
extrait_cif <- function(t) {
  s <- summary(aj, times = t, extend = TRUE)
  k <- summary(km, times = t, extend = TRUE)
  tibble(horizon = t,
         prog_aj = s$pstate[, idx_p], prog_bas = s$lower[, idx_p], prog_haut = s$upper[, idx_p],
         deces_aj = s$pstate[, idx_d], deces_bas = s$lower[, idx_d], deces_haut = s$upper[, idx_d],
         prog_1_km = 1 - k$surv, a_risque = k$n.risk)
}
cif <- extrait_cif(c(5, 10, 20, 30))
tab(cif, "B_incidence.csv")

courbe_aj <- tibble(t = aj$time, `Progression (Aalen-Johansen)` = aj$pstate[, idx_p],
                    `Décès sans progression` = aj$pstate[, idx_d]) |>
  pivot_longer(-t, names_to = "courbe", values_to = "p")
courbe_km <- tibble(t = km$time, courbe = "Progression (1 − Kaplan-Meier, décès censurés)", p = 1 - km$surv)
d_fig <- bind_rows(courbe_aj, courbe_km) |> filter(t <= 30)
p1 <- ggplot(d_fig, aes(t, p, colour = courbe, linetype = courbe)) +
  geom_step(linewidth = .8) +
  scale_colour_manual(values = c("Progression (Aalen-Johansen)" = coul$bleu,
                                 "Progression (1 − Kaplan-Meier, décès censurés)" = coul$orange,
                                 "Décès sans progression" = coul$gris)) +
  scale_linetype_manual(values = c("Progression (Aalen-Johansen)" = 1,
                                   "Progression (1 − Kaplan-Meier, décès censurés)" = 2,
                                   "Décès sans progression" = 1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Années depuis le diagnostic de MGUS", y = "Incidence cumulée")
sauve(p1, "B_incidence.png", 6.8, 4.3)

# 2. Selon le pic monoclonal ---------------------------------------------------
# Test de Gray (cmprsk::cuminc) : compare les incidences cumulées entre groupes
# en présence d'un risque concurrent (le log-rank ne convient pas ici).
bp <- filter(b, !is.na(pic15))
cu <- cuminc(bp$ans, bp$evt_num, group = bp$pic15, cencode = 0)
gray <- tibble(evenement = c("progression", "deces"), stat = cu$Tests[, "stat"], p = cu$Tests[, "pv"])
tab(gray, "B_gray.csv")
ajg <- survfit(Surv(ans, evt) ~ pic15, data = bp, id = id)
sg <- summary(ajg, times = c(10, 20), extend = TRUE)
ip <- match("progression", ajg$states)
cif_g <- tibble(groupe = gsub("pic15=", "", as.character(sg$strata)), horizon = sg$time,
                p = sg$pstate[, ip], bas = sg$lower[, ip], haut = sg$upper[, ip])
tab(cif_g, "B_incidence_groupes.csv")
fg_d <- map_dfr(seq_along(ajg$strata), function(k) {
  x <- ajg[k, ]
  tibble(t = x$time, p = x$pstate[, ip], groupe = gsub("pic15=", "", names(ajg$strata)[k]))
}) |> filter(t <= 30)
p2 <- ggplot(fg_d, aes(t, p, colour = groupe)) + geom_step(linewidth = .85) +
  scale_colour_manual(values = c("Moins de 1,5 g/dL" = coul$gris, "1,5 g/dL ou plus" = coul$bleu)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Années depuis le diagnostic de MGUS", y = "Incidence cumulée de progression")
sauve(p2, "B_incidence_pic.png", 6.2, 3.9)

# 3. Modèles par cause et Fine-Gray --------------------------------------------
# Cox par cause : les décès sont censurés, on estime l'effet sur la vitesse de
# progression chez les patients encore vivants (question étiologique).
# Fine-Gray : finegray() garde les décédés « à risque » avec un poids
# décroissant (fgwt) ; le HR porte sur l'incidence cumulée (risque absolu).
# cluster = id donne une variance robuste, les sujets ayant plusieurs lignes.
cc <- filter(b, !is.na(mspike), !is.na(hgb), !is.na(creat)) |>
  mutate(sexeM = as.integer(sex == "M"))
rhs <- "mspike + age10 + sexeM + hgb + lcreat2"
fp <- coxph(as.formula(paste("Surv(ans, evt == 'progression') ~", rhs)), data = cc)
fd <- coxph(as.formula(paste("Surv(ans, evt == 'deces') ~", rhs)), data = cc)
fgd <- finegray(Surv(ans, evt) ~ ., data = cc |> select(id, ans, evt, mspike, age10, sexeM, hgb, lcreat2),
                etype = "progression", id = id)
fg <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", rhs)), data = fgd,
            weights = fgwt, cluster = id)
lib <- c(mspike = "Pic monoclonal, par g/dL", age10 = "Âge, par 10 ans", sexeM = "Homme",
         hgb = "Hémoglobine, par g/dL", lcreat2 = "Créatinine, par doublement")
tidy <- function(f, m) {
  s <- summary(f); ci <- s$conf.int
  tibble(modele = m, terme = lib[rownames(ci)], hr = ci[, 1], bas = ci[, 3], haut = ci[, 4],
         p = s$coefficients[, ncol(s$coefficients)])
}
modeles <- bind_rows(tidy(fp, "Cox par cause : progression"),
                     tidy(fd, "Cox par cause : décès sans progression"),
                     tidy(fg, "Fine-Gray : progression"))
tab(modeles, "B_modeles.csv")

# Hypothèses et forme
z <- cox.zph(fp)
hyp <- tibble(terme = rownames(z$table), chi2 = z$table[, 1], p = z$table[, 3])
tab(hyp, "B_schoenfeld.csv")
f_ns <- coxph(as.formula("Surv(ans, evt == 'progression') ~ ns(mspike, 3) + age10 + sexeM + hgb + lcreat2"), data = cc)
forme <- tibble(p_non_lineaire = anova(fp, f_ns)[2, "Pr(>|Chi|)"])
tab(forme, "B_forme.csv")

# Sous-groupes
inter <- map_dfr(c("sex", "age70"), function(v) {
  # Pour le sexe, sexeM est retiré : il serait identique à la variable `sex`.
  ajust <- if (v == "sex") "age10 + hgb + lcreat2" else "age10 + sexeM + hgb + lcreat2"
  f <- coxph(as.formula(sprintf("Surv(ans, evt == 'progression') ~ mspike * %s + %s", v, ajust)), data = cc)
  t <- grep("mspike:", names(coef(f)), value = TRUE)
  niv <- levels(cc[[v]])
  fits <- map(niv, function(n) {
    x <- cc[cc[[v]] == n, ]
    rr <- if (v == "sex") "mspike + age10 + hgb + lcreat2" else "mspike + age10 + sexeM + hgb + lcreat2"
    list(f = coxph(as.formula(paste("Surv(ans, evt == 'progression') ~", rr)), data = x), n = nrow(x),
         e = sum(x$evt == "progression"))
  })
  dif <- coef(fits[[2]]$f)["mspike"] - coef(fits[[1]]$f)["mspike"]
  sed <- sqrt(vcov(fits[[2]]$f)["mspike", "mspike"] + vcov(fits[[1]]$f)["mspike", "mspike"])
  map2_dfr(niv, fits, function(n, r) {
    ci <- exp(confint(r$f))["mspike", ]
    tibble(variable = v, niveau = as.character(n), n = r$n, evenements = r$e,
           hr = exp(coef(r$f)["mspike"]), bas = ci[1], haut = ci[2],
           p_commun = summary(f)$coefficients[t, ncol(summary(f)$coefficients)],
           p_stratifie = 2 * pnorm(-abs(dif / sed)))
  })
})
p_inter <- inter |> distinct(variable, p_commun, p_stratifie) |>
  mutate(holm_commun = p.adjust(p_commun, "holm"), holm_stratifie = p.adjust(p_stratifie, "holm"))
inter <- left_join(inter, p_inter, by = c("variable", "p_commun", "p_stratifie"))
tab(inter, "B_sous_groupes.csv")

# Sensibilité : année de diagnostic
fs <- coxph(as.formula(paste("Surv(ans, evt == 'progression') ~", rhs, "+ dxyr")), data = cc)
ci <- exp(confint(fs))["mspike", ]
sens <- tibble(analyse = c("Principale", "Ajustée sur l'année de diagnostic"),
               hr = c(exp(coef(fp)["mspike"]), exp(coef(fs)["mspike"])),
               bas = c(exp(confint(fp))["mspike", 1], ci[1]),
               haut = c(exp(confint(fp))["mspike", 2], ci[2]))
tab(sens, "B_sensibilite.csv")

# 4. Valeur ajoutée pour la prédiction à 10 et 20 ans --------------------------
# CSC() combine les deux Cox par cause pour prédire un risque absolu de
# progression. Score() compare trois modèles emboîtés par AUC et Brier, en
# validation croisée bootstrap (100 réplications) pour limiter l'optimisme
# avec seulement 112 progressions.
cc_df <- as.data.frame(cc)
cc_df$evt_h <- as.integer(cc_df$evt == "progression") + 2L * as.integer(cc_df$evt == "deces")
csc0 <- CSC(Hist(ans, evt_h) ~ age10 + sexeM, data = cc_df, cause = 1)
csc1 <- CSC(Hist(ans, evt_h) ~ mspike + age10 + sexeM, data = cc_df, cause = 1)
csc2 <- CSC(Hist(ans, evt_h) ~ mspike + age10 + sexeM + hgb + lcreat2, data = cc_df, cause = 1)
sc <- Score(list(`Âge et sexe` = csc0, `+ pic monoclonal` = csc1, `+ pic, Hb, créatinine` = csc2),
            formula = Hist(ans, evt_h) ~ 1, data = cc_df, times = c(10, 20), cause = 1,
            metrics = c("auc", "brier"), summary = NULL, plots = NULL,
            cens.model = "km", split.method = "bootcv", B = 100,
            null.model = TRUE, progress.bar = NULL, seed = 20261104)
tab(as.data.frame(sc$AUC$score), "B_auc.csv")
tab(as.data.frame(sc$AUC$contrasts), "B_auc_contraste.csv")
tab(as.data.frame(sc$Brier$score), "B_brier.csv")
tab(as.data.frame(sc$Brier$contrasts), "B_brier_contraste.csv")

# 5. Risques absolus pour des profils types (Fine-Gray) ------------------------
# Sur un modèle Fine-Gray, 1 - survfit() donne directement l'incidence cumulée.
prof <- expand.grid(mspike = c(0.5, 1.5, 2.5), age10 = c(6.0, 7.5), sexeM = 0,
                    hgb = 13, lcreat2 = 0)
sb <- survfit(fg, newdata = prof, se.fit = FALSE)
s_t <- summary(sb, times = c(10, 20), extend = TRUE)$surv
prof$risque10 <- 1 - s_t[1, ]; prof$risque20 <- 1 - s_t[2, ]
tab(prof, "B_profils.csv")

# 6. Décision ------------------------------------------------------------------
ac <- as.data.frame(sc$AUC$contrasts)
ligne <- ac |> filter(times == 10, model == "+ pic monoclonal", reference == "Âge et sexe")
pk <- modeles |> filter(modele == "Cox par cause : progression", terme == lib["mspike"])
dec <- tibble(
  critere = c("Association", "Gain d'AUC à 10 ans"),
  regle = c("IC à 95 % du HR excluant 1", "Gain ≥ 0,05 et IC excluant 0"),
  rempli = c(pk$bas > 1, ligne$delta.AUC >= 0.05 & ligne$lower > 0))
tab(dec, "B_decision.csv")

print(cif); print(gray); print(cif_g); print(modeles, n = 20); print(hyp); print(forme)
print(inter |> select(variable, niveau, n, evenements, hr, bas, haut, p_commun, holm_commun, p_stratifie, holm_stratifie))
print(sens); print(as.data.frame(sc$AUC$score)); print(ac); print(as.data.frame(sc$Brier$score)); print(prof); print(dec)
