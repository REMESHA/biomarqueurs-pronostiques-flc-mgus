# Reproduit l'ensemble du projet depuis les données du package survival.
# À lancer depuis la racine du dépôt : Rscript run_all.R

# Les scripts contiennent des caractères accentués et des symboles (κ, λ, ≥).
# Dans une session R non UTF-8, source() tronque silencieusement les fichiers :
# les scripts s'exécutent en partie et les résultats restent ceux d'une
# exécution précédente. On arrête donc tout de suite dans ce cas.
if (!isTRUE(l10n_info()$`UTF-8`)) {
  stop("La session R n'est pas en UTF-8. Lancer par exemple : ",
       "LANG=C.UTF-8 LC_ALL=C.UTF-8 Rscript run_all.R")
}
options(warn = 1)  # avertissements affichés au moment où ils surviennent

# Les fonctions partagées sont testées avant de s'en servir.
testthat::test_dir("tests/testthat", stop_on_failure = TRUE)

# On repart de zéro : aucun fichier intermédiaire d'une exécution précédente.
unlink("results/derive", recursive = TRUE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

for (s in c("01_donnees", "02_taille_echantillon", "03_flc_mortalite", "04_mgus_progression")) {
  message("== ", s)
  source(file.path("scripts", paste0(s, ".R")), local = new.env(), encoding = "UTF-8")
}
writeLines(capture.output(sessionInfo()), "session_info.txt")
