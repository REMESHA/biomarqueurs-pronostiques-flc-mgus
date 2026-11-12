# Valeur pronostique de biomarqueurs sanguins

Évaluation de deux biomarqueurs sur des données publiques de la Mayo Clinic : les chaînes légères libres (FLC) pour la mortalité, et le pic monoclonal pour la progression de la MGUS. Le projet suit une démarche d'étude complète : plan d'analyse, biais anticipés, valeur ajoutée, critères de décision et communication des résultats.

**À consulter :** [rapport](https://remesha.github.io/biomarqueurs-pronostiques-flc-mgus/) · [plan d'analyse](https://remesha.github.io/biomarqueurs-pronostiques-flc-mgus/plan_analyse_statistique.html) · [synthèse d'une page](https://remesha.github.io/biomarqueurs-pronostiques-flc-mgus/synthese.pdf)

## Résultats principaux

- **FLC :** HR de décès de 1,71 par doublement après ajustement (IC à 95 % 1,58 à 1,85), mais un gain de C de 0,009, sous le seuil prévu de 0,01. Décision selon la règle fixée à l'avance : réorienter.
- **Pic monoclonal :** HR de progression de 2,50 par g/dL (IC à 95 % 1,81 à 3,46). Le gain d'AUC à 10 ans reste trop imprécis pour conclure avec 112 progressions.
- **Risques concurrents :** à 20 ans, censurer les décès double le risque estimé de progression (21 % au lieu de 10 %).

## Méthodes

Modèle de Cox, imputation multiple, splines, correction de l'optimisme par bootstrap, AUC dépendante du temps, score de Brier, calibration en validation croisée, E-value, Aalen-Johansen, Cox par cause et Fine-Gray. Le rapport suit REMARK.

## Reproduire

Avec R (≥ 4.4), Quarto et les packages listés dans `install.R` :

```bash
Rscript install.R
Rscript run_all.R
./render.sh
```

La session R doit être en UTF-8 (sinon `run_all.R` s'arrête avec un message). Sous Linux : `LANG=C.UTF-8 LC_ALL=C.UTF-8 Rscript run_all.R`. Les tests unitaires sont lancés avant les analyses.

Les données proviennent du package R `survival`. Les graines sont fixées dans les scripts.

## Structure

```text
protocole/   plan d'analyse statistique
R/           fonctions et thème des figures
scripts/     préparation, taille d'échantillon, analyses
tests/       tests unitaires
results/     tableaux et figures
rapport/     sources Quarto
docs/        site publié
checklists/  liste de contrôle REMARK
```

## Limites

Analyse rétrospective de données publiques, à visée de démonstration méthodologique. Les biomarqueurs ne sont étudiés que sous l'angle pronostique, et aucune validation externe n'a été réalisée.

Code sous licence MIT.
