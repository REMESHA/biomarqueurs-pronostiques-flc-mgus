# 03 - Partie A : valeur pronostique des FLC pour la mortalité toutes causes
# But : répondre à trois questions successives.
#   1. Association : le FLC est-il lié au décès, à âge, sexe et rein comparables ?
#   2. Valeur ajoutée : améliore-t-il la prédiction par rapport au modèle de base ?
#   3. Robustesse : le résultat tient-il face aux biais anticipés dans le plan ?
# Puis appliquer les critères de décision fixés à l'avance (C1 à C5).
suppressPackageStartupMessages({
  library(survival); library(dplyr); library(tidyr); library(ggplot2)
  library(mice); library(splines); library(riskRegression); library(purrr)
})
source("R/fonctions.R")
set.seed(20261103)
a <- readRDS("results/derive/flchain.rds")
tab <- function(x, f) write.csv(x, file.path("results/tables", f), row.names = FALSE)

# 1. Imputation multiple de la créatinine --------------------------------------
# 17 % de créatinine manquante, surtout en 2002-2003 : les cas complets seraient
# sélectionnés. Le modèle d'imputation contient le critère de jugement sous la
# forme recommandée pour un Cox (indicateur de décès + risque cumulé de
# Nelson-Aalen, White et Royston 2009), sinon l'association serait atténuée.
a$H0 <- nelsonaalen(a, annees, death)
imp_vars <- a |>
  transmute(age, age2 = age^2, sex, periode, log2flc, log2flc2 = log2flc^2,
            ratio_anormal, mgus, death, H0, lcreat)
# Seule lcreat est imputée (pmm : valeurs tirées parmi des donneurs observés,
# donc plausibles). Une seule variable incomplète : une itération suffit.
pred <- quickpred(imp_vars, mincor = 0); pred[, ] <- 0
pred["lcreat", setdiff(names(imp_vars), "lcreat")] <- 1
meth <- make.method(imp_vars); meth[] <- ""; meth["lcreat"] <- "pmm"
M <- 20
imp <- mice(imp_vars, m = M, method = meth, predictorMatrix = pred, maxit = 1,
            printFlag = FALSE, seed = 20261103)
impute <- lapply(1:M, function(k) {
  d <- a
  d$lcreat <- mice::complete(imp, k)$lcreat
  d$creat_i <- exp(d$lcreat)
  d$dfg_i <- ckdepi2021(d$creat_i, d$age, as.character(d$sex))
  d$dfg60 <- factor(ifelse(d$dfg_i < 60, "DFG < 60", "DFG ≥ 60"),
                    levels = c("DFG ≥ 60", "DFG < 60"))
  d
})

base_rhs <- "ns(age, 3) + sex + ns(lcreat, 2) + strata(periode)"
f_base <- as.formula(paste("Surv(annees, death) ~", base_rhs))
f_plein <- update(f_base, . ~ . + log2flc)

# Ajuste le même Cox sur les 20 jeux imputés et combine par les règles de Rubin.
pool_mi <- function(formule, terme, donnees = impute) {
  fits <- lapply(donnees, function(d) coxph(formule, data = d))
  e <- sapply(fits, function(f) coef(f)[terme])
  v <- sapply(fits, function(f) diag(vcov(f))[terme])
  hr_ic(rubin(e, v))
}

# 2. Analyse principale --------------------------------------------------------
# Trois niveaux d'ajustement pour montrer le rôle de la confusion (brut -> âge
# et sexe -> âge, sexe, créatinine), puis le rapport κ/λ ajusté sur le FLC.
principal <- pool_mi(f_plein, "log2flc")
brut <- {
  f <- coxph(Surv(annees, death) ~ log2flc + strata(periode), data = a)
  ci <- exp(confint(f))
  data.frame(hr = exp(coef(f)), bas = ci[1], haut = ci[2],
             p = summary(f)$coefficients[, 5], fmi = NA)
}
ajust_age <- {
  f <- coxph(Surv(annees, death) ~ log2flc + ns(age, 3) + sex + strata(periode),
             data = a)
  ci <- exp(confint(f))["log2flc", ]
  data.frame(hr = exp(coef(f)["log2flc"]), bas = ci[1], haut = ci[2],
             p = summary(f)$coefficients["log2flc", 5], fmi = NA)
}
ratio_sec <- pool_mi(update(f_plein, . ~ . + ratio_anormal), "ratio_anormal")
effets <- bind_rows(
  mutate(brut, modele = "Stratifié sur la période seule"),
  mutate(ajust_age, modele = "Ajusté sur l'âge et le sexe"),
  mutate(principal, modele = "Analyse principale : âge, sexe, créatinine (IM)"),
  mutate(ratio_sec, modele = "Rapport κ/λ anormal, ajusté sur le FLC (IM)")
)
rownames(effets) <- NULL
tab(effets, "A_effets.csv")

# 3. Forme de la relation (cas complets) ---------------------------------------
# Test du rapport de vraisemblance : terme linéaire contre spline à 3 ddl
# (les deux modèles sont emboîtés). La courbe donne le HR par rapport au FLC
# médian ; l'IC vient de la matrice de covariance appliquée aux différences de
# lignes de la matrice de plan (même nœuds de spline que le modèle).
cc <- filter(a, !is.na(lcreat))
m_base <- coxph(f_base, data = cc, x = TRUE, y = TRUE)
m_lin <- coxph(f_plein, data = cc, x = TRUE, y = TRUE)
m_nl <- coxph(update(f_base, . ~ . + ns(log2flc, 3)), data = cc)
lrt <- anova(m_lin, m_nl)
lrt_base <- anova(m_base, m_lin)
forme <- tibble(p_non_lineaire = lrt[2, "Pr(>|Chi|)"],
                chi2_non_lineaire = lrt[2, "Chisq"],
                chi2_flc = lrt_base[2, "Chisq"])
ref <- median(cc$log2flc)
g <- tibble(flc = exp(seq(log(quantile(a$flc, .005)), log(quantile(a$flc, .995)),
                          length.out = 80)))
nd <- tibble(age = 65, sex = factor("F", levels = levels(cc$sex)), lcreat = median(cc$lcreat),
             periode = factor("1996", levels = levels(cc$periode)),
             log2flc = log2(g$flc))
nd0 <- nd[1, ]; nd0$log2flc <- ref
# HR relatif à la médiane, avec IC par la matrice de plan
X <- model.matrix(delete.response(terms(m_nl)), rbind(nd, nd0) |> mutate(periode = "1996"),
                  xlev = m_nl$xlevels)
X <- X[, names(coef(m_nl)), drop = FALSE]
dX <- sweep(X[1:80, , drop = FALSE], 2, X[81, ])
lp <- as.numeric(dX %*% coef(m_nl))
se <- sqrt(rowSums((dX %*% vcov(m_nl)) * dX))
courbe <- g |> mutate(hr = exp(lp), bas = exp(lp - 1.96 * se), haut = exp(lp + 1.96 * se))
tab(courbe, "A_courbe_flc.csv")

# 4. Proportionnalité et effet selon le suivi ----------------------------------
# cox.zph teste si les HR sont constants dans le temps (résidus de Schoenfeld).
# survSplit découpe le suivi en intervalles pour estimer un HR par période ;
# pas de terme principal pour l'intervalle, absorbé par le risque de base.
zph <- cox.zph(m_lin)
split <- survSplit(Surv(annees, death) ~ ., data = cc, cut = c(2, 5, 10),
                   episode = "inter") |>
  mutate(inter = factor(inter, labels = c("0-2 ans", "2-5 ans", "5-10 ans", "> 10 ans")))
m_t <- coxph(Surv(tstart, annees, death) ~ ns(age, 3) + sex + ns(lcreat, 2) +
               log2flc:inter + strata(periode), data = split)
nm <- grep("log2flc", names(coef(m_t)), value = TRUE)
ci <- exp(confint(m_t))[nm, ]
par_suivi <- tibble(intervalle = levels(split$inter), hr = exp(coef(m_t)[nm]),
                    bas = ci[, 1], haut = ci[, 2],
                    deces = as.integer(tapply(split$death, split$inter, sum)))
tab(par_suivi, "A_hr_par_suivi.csv")
hyp <- tibble(terme = rownames(zph$table), chi2 = zph$table[, 1], ddl = zph$table[, 2],
              p = zph$table[, 3])
tab(hyp, "A_schoenfeld.csv")

# 5. Sous-groupes --------------------------------------------------------------
# Deux tests d'interaction sont rapportés :
# - commun : terme produit dans un seul modèle, ajustement commun aux sous-groupes ;
# - stratifié : différence des HR de modèles séparés, où tous les ajustements
#   varient selon le sous-groupe (équivalent à une interaction complète).
sg <- function(var) {
  # Pour le sexe, la variable de sous-groupe remplace `sex` dans l'ajustement :
  # garder les deux crée deux colonnes identiques, dont coxph élimine une.
  ajust <- if (var == "sexe") "ns(age, 3) + ns(lcreat, 2)" else "ns(age, 3) + sex + ns(lcreat, 2)"
  f_sg <- if (var == "sexe") update(f_plein, . ~ . - sex) else f_plein
  fm <- as.formula(sprintf("Surv(annees, death) ~ %s + log2flc * %s + strata(periode)", ajust, var))
  fits <- lapply(impute, function(d) coxph(fm, data = d))
  terme <- grep("log2flc:", names(coef(fits[[1]])), value = TRUE)
  stopifnot(length(terme) == 1)
  p_commun <- rubin(sapply(fits, \(f) coef(f)[terme]), sapply(fits, \(f) diag(vcov(f))[terme]))$p
  niv <- levels(impute[[1]][[var]]); stopifnot(length(niv) == 2)
  par_imp <- lapply(impute, function(d) lapply(niv, function(n) {
    x <- d[d[[var]] == n, ]
    f <- coxph(f_sg, data = x)
    list(e = coef(f)["log2flc"], v = vcov(f)["log2flc", "log2flc"],
         n = nrow(x), ev = sum(x$death))
  }))
  p_strat <- rubin(sapply(par_imp, \(r) r[[2]]$e - r[[1]]$e),
                   sapply(par_imp, \(r) r[[2]]$v + r[[1]]$v))$p
  map_dfr(seq_along(niv), function(k) {
    r <- lapply(par_imp, `[[`, k)
    hr_ic(rubin(sapply(r, `[[`, "e"), sapply(r, `[[`, "v"))) |>
      mutate(variable = var, niveau = niv[k], n = mean(sapply(r, `[[`, "n")),
             evenements = mean(sapply(r, `[[`, "ev")),
             p_commun = p_commun, p_stratifie = p_strat)
  })
}
impute <- lapply(impute, function(d) { d$sexe <- factor(ifelse(d$sex == "F", "Femmes", "Hommes"), levels = c("Femmes", "Hommes")); d })
sous_g <- bind_rows(sg("sexe"), sg("age70"), sg("dfg60"))
p_inter <- sous_g |> distinct(variable, p_commun, p_stratifie) |>
  mutate(holm_commun = p.adjust(p_commun, "holm"),
         holm_stratifie = p.adjust(p_stratifie, "holm"))
sous_g <- left_join(sous_g, p_inter, by = c("variable", "p_commun", "p_stratifie"))
rownames(sous_g) <- NULL
tab(sous_g, "A_sous_groupes.csv")

# 6. Sensibilité ---------------------------------------------------------------
# Chaque analyse répond à un biais du plan : MGUS (clone plasmocytaire), repère
# à 2 ans (causalité inverse), âge comme échelle de temps (forme de
# l'ajustement sur l'âge, qui ne respecte pas la proportionnalité), cas
# complets (imputation), DFG à la place de la créatinine (mesure du rein).
s1 <- pool_mi(f_plein, "log2flc", lapply(impute, \(d) filter(d, mgus == 0)))
rep2 <- lapply(impute, \(d) d |> filter(annees > 2) |> mutate(annees = annees - 2))
s2 <- pool_mi(f_plein, "log2flc", rep2)
s3 <- pool_mi(Surv(age, age + annees, death) ~ sex + ns(lcreat, 2) + log2flc + strata(periode),
              "log2flc")
s4 <- {
  ci <- exp(confint(m_lin))["log2flc", ]
  data.frame(hr = exp(coef(m_lin)["log2flc"]), bas = ci[1], haut = ci[2],
             p = summary(m_lin)$coefficients["log2flc", 5], fmi = NA)
}
s5 <- pool_mi(Surv(annees, death) ~ ns(age, 3) + sex + ns(log(dfg_i), 2) + log2flc + strata(periode),
              "log2flc")
sens <- bind_rows(
  mutate(principal, analyse = "Analyse principale"),
  mutate(s1, analyse = "Sans MGUS"),
  mutate(s2, analyse = "Repère à 2 ans"),
  mutate(s3, analyse = "Âge comme échelle de temps"),
  mutate(s4, analyse = "Cas complets"),
  mutate(s5, analyse = "DFG estimé à la place de la créatinine")
)
rownames(sens) <- NULL
tab(sens, "A_sensibilite.csv")

# E-value : force minimale (échelle RR) qu'un facteur de confusion non mesuré
# devrait avoir avec le FLC et le décès pour expliquer tout le HR par doublement.
ev <- tibble(hr = principal$hr, bas = principal$bas,
             evalue_point = evalue_hr(principal$hr),
             evalue_bas = evalue_hr(principal$bas))
tab(ev, "A_evalue.csv")

# 7. Causes de décès (exploratoire) --------------------------------------------
top <- names(sort(table(a$chapter), decreasing = TRUE))[1:3]
cause <- map_dfr(top, function(ch) {
  fits <- lapply(impute, function(d) {
    d$ev_cause <- as.integer(d$death == 1 & !is.na(d$chapter) & d$chapter == ch)
    coxph(update(f_plein, Surv(annees, ev_cause) ~ .), data = d)
  })
  hr_ic(rubin(sapply(fits, \(x) coef(x)["log2flc"]), sapply(fits, \(x) vcov(x)["log2flc", "log2flc"]))) |>
    mutate(cause = ch, deces = sum(a$chapter == ch, na.rm = TRUE))
}) |> mutate(p_holm = p.adjust(p, "holm"))
rownames(cause) <- NULL
tab(cause, "A_causes.csv")

# 8. Discrimination : C de Harrell ---------------------------------------------
# concordance(m_base, m_lin) estime les deux C sur les mêmes paires avec leur
# covariance : l'IC de la différence en tient compte. Le bootstrap (200
# réplications) mesure l'optimisme : performance sur l'échantillon bootstrap
# moins performance du même modèle sur les données d'origine.
cmp <- concordance(m_base, m_lin)
cc_est <- coef(cmp); V <- vcov(cmp)
dc <- unname(cc_est[2] - cc_est[1]); se_dc <- sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2])
B <- 200
# Boucle au niveau global : les formules y résolvent `bd` lors de la prédiction.
opt <- tibble(o0 = numeric(B), o1 = numeric(B))
for (b in 1:B) {
  bd <- cc[sample(nrow(cc), replace = TRUE), ]
  f0 <- coxph(f_base, data = bd); f1 <- coxph(f_plein, data = bd)
  opt$o0[b] <- concordance(f0, newdata = bd)$concordance -
    concordance(f0, newdata = cc)$concordance
  opt$o1[b] <- concordance(f1, newdata = bd)$concordance -
    concordance(f1, newdata = cc)$concordance
}
disc <- tibble(
  modele = c("Base", "Base + FLC"), c_apparent = cc_est,
  optimisme = c(mean(opt$o0), mean(opt$o1)),
  c_corrige = cc_est - c(mean(opt$o0), mean(opt$o1)),
  se = sqrt(diag(V)))
dc_tab <- tibble(delta_c = dc, se = se_dc, bas = dc - 1.96 * se, haut = dc + 1.96 * se,
                 delta_c_corrige = disc$c_corrige[2] - disc$c_corrige[1])
tab(disc, "A_concordance.csv"); tab(dc_tab, "A_delta_c.csv")

# Exploratoire, décidé après avoir vu la non-linéarité : gain de C avec une spline.
cmp_nl <- concordance(m_base, m_nl); e_nl <- coef(cmp_nl); V_nl <- vcov(cmp_nl)
d_nl <- unname(e_nl[2] - e_nl[1]); se_nl <- sqrt(V_nl[1, 1] + V_nl[2, 2] - 2 * V_nl[1, 2])
tab(tibble(delta_c = d_nl, bas = d_nl - 1.96 * se_nl, haut = d_nl + 1.96 * se_nl),
    "A_delta_c_spline.csv")

# 9. AUC dépendante du temps et Brier ------------------------------------------
# Score() (riskRegression) pondère par l'inverse de la probabilité de censure
# (IPCW). Mesures apparentes, sans validation croisée : l'optimisme mesuré sur
# le C est inférieur à 0,001 (7 paramètres, près de 2 000 décès).
# La strate 2002-2003 (25 cas complets) n'a pas 10 ans de suivi : pas de risque
# prédictible à cet horizon. Elle est retirée de cette section et de la suivante.
cc_perf <- filter(cc, periode != "2002-2003") |> droplevels()
m_base_p <- coxph(f_base, data = cc_perf, x = TRUE, y = TRUE)
m_lin_p <- coxph(f_plein, data = cc_perf, x = TRUE, y = TRUE)
sc <- Score(list(Base = m_base_p, `Base + FLC` = m_lin_p),
            formula = Surv(annees, death) ~ 1, data = cc_perf, times = c(5, 10),
            metrics = c("auc", "brier"), summary = NULL, plots = NULL,
            cens.model = "km", conf.int = TRUE, null.model = TRUE, se.fit = TRUE)
auc <- as.data.frame(sc$AUC$score)
auc_c <- as.data.frame(sc$AUC$contrasts)
bri <- as.data.frame(sc$Brier$score)
bri_c <- as.data.frame(sc$Brier$contrasts)
tab(auc, "A_auc_t.csv"); tab(auc_c, "A_auc_contraste.csv")
tab(bri, "A_brier.csv"); tab(bri_c, "A_brier_contraste.csv")

# 10. Calibration en validation croisée à 10 plis ------------------------------
# Chaque sujet reçoit un risque prédit par un modèle qui ne l'a pas vu. On
# compare ensuite, par décile, le risque moyen prédit au risque observé (KM).
pli <- sample(rep(1:10, length.out = nrow(cc_perf)))
pred_cv <- matrix(NA_real_, nrow(cc_perf), 2, dimnames = list(NULL, c("t5", "t10")))
for (k in 1:10) {
  ent <- cc_perf[pli != k, ]; test <- cc_perf[pli == k, ]
  f <- coxph(f_plein, data = ent, x = TRUE, y = TRUE)
  pred_cv[pli == k, ] <- predictRisk(f, newdata = test, times = c(5, 10))
}
calib <- map_dfr(c(5, 10), function(t) {
  p <- pred_cv[, paste0("t", t)]
  dec <- cut(p, unique(quantile(p, 0:10 / 10, na.rm = TRUE)), include.lowest = TRUE, labels = FALSE)
  map_dfr(sort(unique(dec)), function(q) {
    i <- which(dec == q)
    s <- summary(survfit(Surv(annees, death) ~ 1, data = cc_perf[i, ]), times = t, extend = TRUE)
    tibble(horizon = t, decile = q, predit = mean(p[i]), observe = 1 - s$surv,
           bas = 1 - s$upper, haut = 1 - s$lower, n = length(i))
  })
})
tab(calib, "A_calibration.csv")
# Pente : régression du log(-log(1 - risque CV)) sur la survie jusqu'à l'horizon.
pente <- map_dfr(c(5, 10), function(t) {
  p <- pmin(pmax(pred_cv[, paste0("t", t)], 1e-6), 1 - 1e-6)
  cll <- log(-log(1 - p))
  d <- cc_perf |> mutate(tt = pmin(annees, t), ev = as.integer(death == 1 & annees <= t), cll = cll)
  f <- coxph(Surv(tt, ev) ~ cll + strata(periode), data = d)
  ci <- confint(f)
  tibble(horizon = t, pente = unname(coef(f)), bas = ci[1], haut = ci[2],
         oe = (1 - summary(survfit(Surv(annees, death) ~ 1, data = cc_perf), times = t)$surv) / mean(p))
})
tab(pente, "A_pente_calibration.csv")

# 11. Risques absolus pour un profil type --------------------------------------
pq <- quantile(cc$flc, c(.25, .5, .75, .95))
prof <- tibble(flc = pq, quantile = names(pq), age = 65,
               sex = factor("F", levels = levels(cc$sex)),
               lcreat = log(1), periode = factor("1996", levels = levels(cc$periode)),
               log2flc = log2(pq))
prof$risque10 <- as.numeric(predictRisk(m_lin_p, newdata = as.data.frame(prof), times = 10))
tab(prof |> select(quantile, flc, risque10), "A_profils.csv")

# 12. Décision -----------------------------------------------------------------
# Application mécanique des seuils du plan : aucune décision n'est prise à la main.
crit <- tibble(
  critere = c("C1. Association", "C2. Ampleur", "C3. Valeur ajoutée", "C4. Calibration", "C5. Robustesse"),
  regle = c("IC à 95 % du HR excluant 1", "HR par doublement ≥ 1,20",
            "Gain de C ≥ 0,01 et IC excluant 0", "Pente entre 0,90 et 1,10 à 5 et 10 ans",
            "HR ≥ 1,10 et IC excluant 1 dans toutes les sensibilités"),
  rempli = c(principal$bas > 1, principal$hr >= 1.2,
             dc >= 0.01 & (dc - 1.96 * se_dc) > 0,
             all(pente$pente >= .9 & pente$pente <= 1.1),
             all(sens$hr >= 1.1 & sens$bas > 1))
)
tab(crit, "A_decision.csv")
tab(forme, "A_forme.csv")

# 13. Comparaison à la publication source (Dispenzieri, 2012) --------------------
# Sous-échantillon public : on vérifie l'ordre de grandeur, pas l'égalité exacte.
pub <- filter(cc, mgus == 0) |> mutate(haut_flc = as.integer(flc >= 4.72))
m_p1 <- coxph(Surv(annees, death) ~ haut_flc, data = pub)
m_p2 <- coxph(Surv(annees, death) ~ haut_flc + age + sex + creatinine, data = pub)
ci1 <- exp(confint(m_p1)); ci2 <- exp(confint(m_p2))
tab(tibble(analyse = c("Univariée", "Ajustée sur âge, sexe, créatinine"),
           hr = c(exp(coef(m_p1)), exp(coef(m_p2)["haut_flc"])),
           bas = c(ci1[1, 1], ci2["haut_flc", 1]), haut = c(ci1[1, 2], ci2["haut_flc", 2]),
           hr_publie = c(4.4, 2.1), bas_publie = c(4.1, 1.9), haut_publie = c(4.7, 2.2),
           n = nrow(pub), part_haut = mean(pub$haut_flc)), "A_comparaison_publication.csv")

# Figures ----------------------------------------------------------------------
source("R/theme.R")
p1 <- ggplot(courbe, aes(flc, hr)) +
  geom_hline(yintercept = 1, colour = "grey60", linewidth = .4) +
  geom_ribbon(aes(ymin = bas, ymax = haut), fill = coul$bleu, alpha = .18) +
  geom_line(colour = coul$bleu, linewidth = .9) +
  scale_x_log10(breaks = c(1.5, 2, 3, 5, 8, 12, 20)) + scale_y_log10(breaks = c(.8, 1, 1.25, 1.5, 2, 3, 4)) +
  labs(x = "FLC total, κ + λ (mg/dL, échelle log)", y = "HR par rapport au FLC médian (échelle log)")
sauve(p1, "A_courbe_flc.png", 6.5, 4)

cal <- calib |> mutate(h = factor(paste(horizon, "ans"), levels = c("5 ans", "10 ans")))
p2 <- ggplot(cal, aes(predit, observe)) +
  geom_abline(colour = "grey60", linewidth = .4) +
  geom_errorbar(aes(ymin = bas, ymax = haut), width = 0, colour = coul$bleu, alpha = .6) +
  geom_point(colour = coul$bleu, size = 1.8) +
  facet_wrap(~h) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_x_continuous(labels = scales::percent) + scale_y_continuous(labels = scales::percent) +
  labs(x = "Risque de décès prédit (validation croisée)", y = "Risque observé (Kaplan-Meier)")
sauve(p2, "A_calibration.png", 7, 3.8)

fp <- bind_rows(
  sens |> transmute(groupe = "Sensibilité", libelle = analyse, hr, bas, haut),
  sous_g |> transmute(groupe = "Sous-groupes", libelle = niveau, hr, bas, haut),
  par_suivi |> transmute(groupe = "Temps de suivi", libelle = intervalle, hr, bas, haut)
) |> mutate(libelle = factor(libelle, levels = rev(unique(libelle))),
            groupe = factor(groupe, levels = c("Sensibilité", "Sous-groupes", "Temps de suivi")))
p3 <- ggplot(fp, aes(hr, libelle)) +
  geom_vline(xintercept = 1, colour = "grey60", linewidth = .4) +
  geom_errorbar(aes(xmin = bas, xmax = haut), width = 0, colour = coul$bleu, orientation = "y") +
  geom_point(colour = coul$bleu, size = 2) +
  facet_grid(groupe ~ ., scales = "free_y", space = "free_y") +
  scale_x_log10(breaks = c(1, 1.25, 1.5, 2, 2.5, 3)) +
  labs(x = "HR de décès par doublement du FLC (IC à 95 %, échelle log)", y = NULL)
sauve(p3, "A_foret.png", 7.2, 5)

saveRDS(list(zph = zph$table), "results/derive/A_divers.rds")
cat("\n---\n"); print(effets); print(forme); print(par_suivi); print(sens); print(sous_g |> select(variable, niveau, hr, bas, haut, p_commun, holm_commun, p_stratifie, holm_stratifie)); print(cause |> select(cause, hr, bas, haut, deces, p_holm)); print(disc); print(dc_tab); print(d_nl); print(auc); print(auc_c); print(bri); print(pente); print(crit); print(prof |> select(quantile, flc, risque10)); print(ev); print(hyp)
