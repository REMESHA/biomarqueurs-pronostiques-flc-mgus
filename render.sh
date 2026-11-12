#!/usr/bin/env bash
# Produit le site publié dans docs/ à partir des résultats.
set -euo pipefail
cd "$(dirname "$0")/rapport"
quarto render rapport.qmd --to html --output-dir ../docs --output index.html
cp ../results/figures/B_incidence.png fig_incidence.png  # typst ne lit pas hors du dossier
quarto render synthese.qmd --to typst --output-dir ../docs
rm -f fig_incidence.png
cd ../protocole
quarto render plan_analyse_statistique.md --to html --output-dir ../docs \
  -M embed-resources:true -M lang:fr -M title:"Plan d'analyse statistique"
