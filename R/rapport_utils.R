# Formatage des nombres pour le rapport (notation française).
virg <- function(x, d = 2) gsub(".", ",", sprintf(paste0("%.", d, "f"), x), fixed = TRUE)
pct <- function(x, d = 1) paste0(virg(100 * x, d), " %")
ic <- function(b, h, d = 2) paste0(virg(b, d), " à ", virg(h, d))
lire <- function(f) read.csv(file.path("..", "results", "tables", f), check.names = FALSE,
                             encoding = "UTF-8")
pv <- function(p) ifelse(p < 0.001, "< 0,001", virg(p, 3))
oui_non <- function(x) ifelse(as.logical(x), "Rempli", "Non rempli")
