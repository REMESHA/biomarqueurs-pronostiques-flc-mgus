
**Projet :** valeur pronostique de biomarqueurs sanguins dans la population générale âgée (chaînes légères libres) et dans la MGUS (pic monoclonal)  
**Auteur :** Guy Darcy Remesha  
**Version :** 1.0, 3 novembre 2026

> **Statut de ce plan.** Il a été rédigé avant l'ajustement de tout modèle sur les critères de jugement. Seules la documentation et la structure des données avaient été consultées (variables, données manquantes, répartition des années de prélèvement). Les données étant publiques et les résultats déjà publiés, il ne s'agit pas d'un plan préenregistré au sens strict. Toute déviation est consignée à la fin du rapport.

## 1. Contexte et cas d'usage

Deux cas d'usage pronostiques sont étudiés :

| | Partie A | Partie B |
|---|---|---|
| **Données** | `survival::flchain` | `survival::mgus2` |
| **Population** | Résidents du comté d'Olmsted de 50 ans ou plus, échantillon de 7 874 sujets | 1 384 patients consécutifs atteints de MGUS |
| **Biomarqueur principal** | Chaînes légères libres sériques totales (FLC = κ + λ) | Taille du pic monoclonal (M-spike) |
| **Question** | Le FLC apporte-t-il une information sur la mortalité au-delà de l'âge, du sexe et de la fonction rénale ? | Le pic monoclonal prédit-il la progression, en tenant compte du décès survenu avant celle-ci ? |

Il s'agit d'évaluations **pronostiques**. L'absence de traitement randomisé interdit toute conclusion sur une valeur **prédictive** de la réponse à un traitement.

## 2. Partie A : FLC et mortalité

### 2.1 Population

Sont inclus tous les sujets ayant un suivi strictement positif. Les 3 sujets prélevés le jour de leur décès sont exclus, car leur temps de suivi est nul.

### 2.2 Critères et variables

- **Critère principal :** décès toutes causes, en années depuis le prélèvement.
- **Biomarqueur principal :** log2(κ + λ). L'effet est exprimé pour un doublement du FLC.
- **Biomarqueur secondaire :** rapport κ/λ anormal (< 0,26 ou > 1,65), ajusté sur le FLC.
- **Ajustement :** âge (spline naturelle, 3 ddl), sexe, log(créatinine) (spline naturelle, 2 ddl), avec stratification sur la période de prélèvement (1995, 1996, 1997, 1998-2001, 2002-2003).

La stratification sur la période est justifiée par le mode de recrutement. Les sujets ont d'abord été recrutés lors de consultations, puis par échantillonnage aléatoire de ceux qui n'étaient pas encore venus.

### 2.3 Biais et facteurs de confusion anticipés

| Risque | Mécanisme | Réponse prévue |
|---|---|---|
| Confusion par la fonction rénale | L'insuffisance rénale réduit l'élimination des FLC et augmente la mortalité | Ajustement flexible sur la créatinine, sous-groupes selon le DFG estimé |
| Causalité inverse | Une maladie préclinique peut élever le FLC | Analyse de repère à 2 ans |
| Sélection par données manquantes | La créatinine manque pour environ 90 % des sujets prélevés en 2002-2003, contre 11 à 21 % avant | Imputation multiple incluant la période et le critère de jugement |
| Confusion résiduelle | Comorbidités non mesurées | E-value, interprétation non causale |
| Hétérogénéité dans le temps | L'effet peut décroître avec le suivi | Test de Schoenfeld, HR par intervalle de suivi |
| Population spécifique | Population surtout blanche d'un seul comté américain | Limite déclarée, validation externe requise |

### 2.4 Analyses

**Analyse principale.** Modèle de Cox multivarié après imputation multiple de la créatinine (20 imputations, *predictive mean matching*). Le modèle d'imputation inclut l'âge, l'âge au carré, le sexe, la période, le FLC, son carré, le rapport κ/λ, la MGUS, l'indicateur de décès et le risque cumulé de Nelson-Aalen. Les résultats sont combinés selon les règles de Rubin.

**Forme de la relation.** Comparaison du terme linéaire à une spline à 3 ddl (test du rapport de vraisemblance sur cas complets).

**Valeur ajoutée, sur cas complets pour la créatinine :**

- C de Harrell du modèle de base et du modèle avec FLC, avec différence et correction de l'optimisme par bootstrap (200 réplications) ;
- AUC dépendante du temps et score de Brier à 5 et 10 ans (IPCW) ;
- calibration à 5 et 10 ans sur prédictions en validation croisée à 10 plis (déciles, pente de calibration).

Le NRI n'est pas utilisé, car il peut suggérer un gain pour un marqueur sans valeur. La courbe de décision n'est pas retenue, faute de décision clinique définie dans ce cas d'usage.

**Sous-groupes.** Interactions prédéfinies avec le sexe, l'âge (< 70 / ≥ 70 ans) et le DFG estimé CKD-EPI 2021 (< 60 / ≥ 60 mL/min/1,73 m²). Les p d'interaction sont corrigés par la méthode de Holm.

**Sensibilité.** (1) exclusion de la MGUS ; (2) repère à 2 ans ; (3) âge comme échelle de temps ; (4) cas complets ; (5) DFG estimé à la place de la créatinine.

**Exploratoire.** HR par cause de décès pour les trois chapitres les plus fréquents (Holm).

### 2.5 Critères de décision définis à l'avance

| Critère | Seuil |
|---|---|
| C1. Association | IC à 95 % du HR par doublement excluant 1 |
| C2. Ampleur | HR par doublement ≥ 1,20 |
| C3. Valeur ajoutée | Gain de C ≥ 0,01 avec IC à 95 % excluant 0 |
| C4. Calibration | Pente de calibration entre 0,90 et 1,10 à 5 et 10 ans |
| C5. Robustesse | Dans toutes les analyses de sensibilité, HR ≥ 1,10 et IC excluant 1 |

- **Poursuivre** vers une validation externe si les cinq critères sont remplis.
- **Réorienter** si C1 et C2 sont remplis mais pas C3 : le marqueur est associé au risque, sans utilité démontrée pour la prédiction individuelle dans cette population.
- **Arrêter** si C1 ou C2 n'est pas rempli.
- Si C1, C2 et C3 sont remplis mais pas C4 ou C5, le signal est conservé et le point en échec doit être investigué avant de poursuivre.

Ces seuils sont des choix de travail. Le seuil de 0,01 sur le C est une convention, à adapter au cas d'usage.

## 3. Partie B : pic monoclonal et progression dans la MGUS

### 3.1 Population et critères

Tous les patients sont inclus. L'analyse multivariée porte sur les cas complets, les données manquantes attendues restant inférieures à 5 %.

- **Événement d'intérêt :** progression vers une hémopathie plasmocytaire (`pstat = 1`).
- **Événement concurrent :** décès sans progression préalable.
- **Censure :** vivant sans progression à la dernière nouvelle.
- **Horizons :** 10 et 20 ans.

### 3.2 Pourquoi un modèle à risques concurrents

Censurer les décès revient à estimer un risque dans un monde où personne ne meurt d'autre chose. Dans une population âgée, cela surestime le risque réel de progression.

### 3.3 Analyses

1. **Incidence cumulée :** estimateur d'Aalen-Johansen, comparé à 1 − Kaplan-Meier. Comparaison selon M-spike ≥ 1,5 g/dL par le test de Gray.
2. **Question étiologique :** modèles de Cox par cause pour la progression et le décès sans progression, avec M-spike, âge, sexe, hémoglobine et log2(créatinine). Chaque covariable consomme 1 ddl.
3. **Prédiction du risque absolu :** modèle de Fine-Gray avec les mêmes covariables.
4. **Valeur ajoutée :** AUC dépendante du temps et score de Brier à 10 et 20 ans pour l'âge et le sexe seuls, puis avec les biomarqueurs (Cox par cause, validation croisée bootstrap).
5. **Forme et hypothèses :** spline à 3 ddl pour le M-spike, résidus de Schoenfeld.
6. **Sous-groupes :** interactions du M-spike avec le sexe et l'âge (< 70 / ≥ 70 ans), avec correction de Holm.
7. **Sensibilité :** ajustement sur l'année de diagnostic.

### 3.4 Critère de décision

Le pic monoclonal est retenu comme marqueur de stratification si le HR par cause par g/dL a un IC à 95 % excluant 1 **et** si l'AUC à 10 ans gagne au moins 0,05 par rapport au modèle âge et sexe, avec un IC de la différence excluant 0.

## 4. Taille d'échantillon

Les effectifs sont imposés par des données existantes. La justification porte donc sur ce qu'ils permettent de détecter ou de modéliser :

- **Partie A :** formule de Schoenfeld pour le nombre de décès nécessaire et le plus petit HR détectable, avec une puissance de 90 % et α bilatéral de 5 %, en tenant compte de la corrélation du biomarqueur avec les covariables ;
- **Partie B :** nombre d'événements par paramètre et critères de Riley (`pmsampsize`) pour un modèle de survie, avec un R² de Cox-Snell hypothétique.

Aucune puissance a posteriori n'est calculée à partir de l'effet observé.

## 5. Rapport et logiciel

- Recommandations REMARK, avec une liste de contrôle dans `checklists/`.
- R, avec packages et graine fixés. Le test α est bilatéral à 5 %.
- Les intervalles de confiance sont privilégiés. Les analyses secondaires, de sous-groupes et exploratoires sont présentées comme telles.
