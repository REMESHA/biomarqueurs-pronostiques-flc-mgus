# 01 - Préparation des données
# But : construire les deux jeux d'analyse (flchain pour la partie A, mgus2 pour
# la partie B), créer les variables dérivées, documenter le flux de sujets et
# produire les tableaux descriptifs. Les jeux préparés sont sauvegardés en .rds
# pour que les scripts suivants partent exactement des mêmes données.
suppressPackageStartupMessages({library(survival); library(dplyr)})
source("R/fonctions.R")
dir.create("results/derive", showWarnings = FALSE, recursive = TRUE)

# Partie A : flchain ----------------------------------------------------------
flchain <- survival::flchain
flux_a <- c(`Sujets dans flchain` = nrow(flchain),
            `Exclus : suivi nul` = sum(flchain$futime == 0))
a <- flchain |>
  # Un suivi nul (prélèvement le jour du décès) n'apporte aucune information de
  # survie et pose problème dans un modèle de Cox : ces sujets sont exclus.
  filter(futime > 0) |>
  mutate(
    id = row_number(),
    annees = futime / 365.25,                  # jours -> années
    flc = kappa + lambda,                      # biomarqueur : FLC total
    log2flc = log2(flc),                       # effet lu par doublement
    ratio = kappa / lambda,
    ratio_anormal = as.integer(ratio < 0.26 | ratio > 1.65),  # bornes usuelles
    lcreat = log(creatinine),                  # créatinine très asymétrique
    # Période de prélèvement : le mode de recrutement a changé au fil du temps.
    periode = cut(sample.yr, c(1994, 1995, 1996, 1997, 2001, 2003),
                  labels = c("1995", "1996", "1997", "1998-2001", "2002-2003")),
    dfg = ckdepi2021(creatinine, age, as.character(sex)),
    age70 = factor(ifelse(age >= 70, "70 ans ou plus", "Moins de 70 ans"),
                   levels = c("Moins de 70 ans", "70 ans ou plus")),
    groupe_flc = cut(flc, quantile(flc, 0:4 / 4), include.lowest = TRUE,
                     labels = paste0("Q", 1:4))
  )
flux_a <- c(flux_a, `Inclus` = nrow(a),
            `Décès` = sum(a$death),
            `Créatinine manquante` = sum(is.na(a$creatinine)),
            `Cas complets` = sum(!is.na(a$creatinine)),
            `Décès, cas complets` = sum(a$death[!is.na(a$creatinine)]))

# Contrôles de cohérence : si le package survival change ses données, le
# pipeline s'arrête ici au lieu de produire un rapport aux chiffres faux.
stopifnot(nrow(a) == 7871, sum(a$death) == 2166, all(a$flc > 0),
          all(is.na(a$chapter) == (a$death == 0)))

decrit <- function(d) tibble(
  n = nrow(d),
  `Âge médian (ans)` = median(d$age),
  `Femmes (%)` = 100 * mean(d$sex == "F"),
  `FLC médian (mg/dL)` = median(d$flc),
  `Créatinine médiane (mg/dL)` = median(d$creatinine, na.rm = TRUE),
  `Rapport κ/λ anormal (%)` = 100 * mean(d$ratio_anormal),
  `MGUS (%)` = 100 * mean(d$mgus),
  `Créatinine manquante (%)` = 100 * mean(is.na(d$creatinine)),
  `Décès (%)` = 100 * mean(d$death),
  `Suivi médian (ans)` = median(d$annees)
)
tab_a <- bind_rows(mutate(decrit(a), groupe = "Ensemble"),
                   a |> group_by(groupe = as.character(groupe_flc)) |>
                     group_modify(~decrit(.x)) |> ungroup())
write.csv(tab_a, "results/tables/A_descriptif.csv", row.names = FALSE)

# Données manquantes selon la période : justifie l'imputation et la strate.
mq <- a |> group_by(periode) |>
  summarise(n = n(), manquante_pct = 100 * mean(is.na(creatinine)),
            age_median = median(age), deces_pct = 100 * mean(death))
write.csv(mq, "results/tables/A_manquants_par_periode.csv", row.names = FALSE)

# Rapport κ/λ anormal et MGUS : dans ce fichier, les deux se recouvrent presque
# entièrement. Le tableau sert à interpréter l'analyse secondaire du rapport.
rm_tab <- a |> count(ratio_anormal, mgus) |>
  mutate(flc_median = sapply(seq_len(n()), function(i)
    median(a$flc[a$ratio_anormal == ratio_anormal[i] & a$mgus == mgus[i]])))
write.csv(rm_tab, "results/tables/A_ratio_mgus.csv", row.names = FALSE)

# Partie B : mgus2 ------------------------------------------------------------
mgus2 <- survival::mgus2
stopifnot(all(mgus2$ptime <= mgus2$futime))
b <- mgus2 |>
  mutate(
    # Temps et type du premier événement : la progression si elle a lieu,
    # sinon le décès ou la dernière nouvelle.
    etime = ifelse(pstat == 1, ptime, futime),
    evt = case_when(pstat == 1 ~ "progression",
                    death == 1 ~ "deces",
                    TRUE ~ "censure"),
    # Le premier niveau du facteur est la censure : survfit() le reconnaît
    # comme tel pour l'estimateur d'Aalen-Johansen.
    evt = factor(evt, levels = c("censure", "progression", "deces")),
    evt_num = as.integer(evt) - 1L,            # 0, 1, 2 pour cmprsk::cuminc
    ans = etime / 12,                          # mois -> années
    age10 = age / 10,
    lcreat2 = log2(creat),
    pic15 = factor(ifelse(mspike >= 1.5, "1,5 g/dL ou plus", "Moins de 1,5 g/dL"),
                   levels = c("Moins de 1,5 g/dL", "1,5 g/dL ou plus")),
    age70 = factor(ifelse(age >= 70, "70 ans ou plus", "Moins de 70 ans"),
                   levels = c("Moins de 70 ans", "70 ans ou plus"))
  )
flux_b <- c(`Patients dans mgus2` = nrow(b),
            `Exclus : temps nul` = sum(b$etime <= 0))
b <- filter(b, etime > 0)
cc_b <- !is.na(b$mspike) & !is.na(b$hgb) & !is.na(b$creat)
flux_b <- c(flux_b, `Inclus` = nrow(b),
            `Progressions` = sum(b$evt == "progression"),
            `Décès sans progression` = sum(b$evt == "deces"),
            `Censurés` = sum(b$evt == "censure"),
            `Cas complets` = sum(cc_b),
            `Progressions, cas complets` = sum(b$evt[cc_b] == "progression"),
            `Décès sans progression, cas complets` = sum(b$evt[cc_b] == "deces"))
stopifnot(nrow(b) == 1384, sum(b$evt == "progression") == 115, sum(cc_b) == 1338)

decrit_b <- function(d) tibble(
  n = nrow(d),
  `Âge médian (ans)` = median(d$age),
  `Femmes (%)` = 100 * mean(d$sex == "F"),
  `M-spike médian (g/dL)` = median(d$mspike, na.rm = TRUE),
  `Hémoglobine médiane (g/dL)` = median(d$hgb, na.rm = TRUE),
  `Créatinine médiane (mg/dL)` = median(d$creat, na.rm = TRUE),
  `Progression (%)` = 100 * mean(d$evt == "progression"),
  `Décès sans progression (%)` = 100 * mean(d$evt == "deces")
)
tab_b <- bind_rows(
  mutate(decrit_b(b), groupe = "Ensemble"),
  b |> filter(!is.na(pic15)) |> group_by(groupe = as.character(pic15)) |>
    group_modify(~decrit_b(.x)) |> ungroup()
)
write.csv(tab_b, "results/tables/B_descriptif.csv", row.names = FALSE)

flux <- rbind(data.frame(partie = "A", etape = names(flux_a), n = unname(flux_a)),
              data.frame(partie = "B", etape = names(flux_b), n = unname(flux_b)))
write.csv(flux, "results/tables/flux.csv", row.names = FALSE)
saveRDS(a, "results/derive/flchain.rds")
saveRDS(b, "results/derive/mgus2.rds")
print(flux); print(mq); print(rm_tab)
