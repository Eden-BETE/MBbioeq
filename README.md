# MBbioeq — Model-Based Bioequivalence in Pharmacokinetics

An R package providing a unified, open-source framework for model-based bioequivalence (BE) assessment using nonlinear mixed-effects models (NLMEM).

## Overview

**MBbioeq** implements model-based bioequivalence analysis for parallel-group pharmacokinetic studies. It combines:

- **NLMEM estimation** via the SAEM algorithm (through `saemix`)
- **Standard error computation** using the Fisher Information Matrix and the delta method
- **Bioequivalence testing** via TOST (Two One-Sided Tests) and BOT (Bioequivalence Optimal Test)
- **Gallant-corrected confidence intervals** for small-sample uncertainty

## Installation

```r
# From GitHub
devtools::install_github("Eden-BETE/MBbioeq")
```

## Usage

```r
library(MBbioeq)

# Two One-Sided Tests (TOST)
result <- run_MBTost(
  data_path    = "your_data.txt",
  model        = "1cpt",          # or "1cpt_tlag"
  psi0         = c(CL = 5, V = 50, ka = 1),
  error_model  = "combined"
)

# Bioequivalence Optimal Test (BOT)
result <- run_MBbot(
  data_path   = "your_data.txt",
  model       = "1cpt",
  psi0        = c(CL = 5, V = 50, ka = 1),
  error_model = "combined"
)
```

## Supported Models

| Model | Description |
|-------|-------------|
| `1cpt` | One-compartment, first-order absorption |
| `1cpt_tlag` | One-compartment with absorption lag time |

## Key Features

- Works with **sparse sampling designs** (few observations per subject)
- **No commercial software required** (unlike mbbe which requires NONMEM)
- Reports both **Normal** and **Gallant-corrected** 90% confidence intervals
- Separate BE decisions for **AUC** and **Cmax**
- Fully reproducible workflow in R

## References

The methodological framework implemented in **MBbioeq** is based on the following works:

### Core methodology

- **Dubois A, Lavielle M, Gsteiger S, Pigeolet E, Mentré F** (2011).
  *Model-based analyses of bioequivalence crossover trials using the stochastic approximation expectation maximisation algorithm.*
  Statistics in Medicine, **30**(21), 2582–2600.
  https://doi.org/10.1002/sim.4286

- **Loingeville F, Rakez M, Nguyen TT, Donnelly M, Fang L, Feng K, Zhao L, Grosser S, Sun G, Sun W, Mentré F, Bertrand J** (2025).
  *Model-based approach for two-stage group sequential or adaptive designs in bioequivalence studies using parallel and crossover designs.*
  Statistical Methods in Medical Research, pp. 1–14.
  https://doi.org/10.1177/09622802251534925

### Bioequivalence testing

- **Schuirmann DJ** (1987).
  *A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability.*
  Journal of Pharmacokinetics and Biopharmaceutics, **15**(6), 657–680.

- **US Food and Drug Administration** (2021).
  *Bioequivalence Studies with Pharmacokinetic Endpoints for Drugs Submitted under an ANDA.*
  FDA Guidance for Industry.
  https://www.fda.gov/media/87219/download

- **US Food and Drug Administration** (2022).
  *Statistical Approaches to Establishing Bioequivalence.*
  FDA Guidance for Industry.
  https://www.fda.gov/regulatory-information/search-fda-guidance-documents/statistical-approaches-establishing-bioequivalence-0

- **European Medicines Agency** (2001).
  *Note for Guidance on the Investigation of Bioavailability and Bioequivalence.*
  CPMP/EWP/QWP/1401/98.
  https://www.ema.europa.eu/en/documents/scientific-guideline/note-guidance-investigation-bioavailability-bioequivalence_en.pdf

### SAEM algorithm

- **Delyon B, Lavielle M, Moulines E** (1999).
  *Convergence of a stochastic approximation version of the EM algorithm.*
  The Annals of Statistics, **27**(1), 94–128.

- **Comets E, Brendel K, Mentré F** (2017).
  *saemix: SAEM Algorithm for Nonlinear Mixed-Effects Models.*
  R package version 4.2.0.
  https://CRAN.R-project.org/package=saemix

### Standard errors & delta method

- **Oehlert GW** (1992).
  *A note on the delta method.*
  The American Statistician, **46**(1), 27–29.

- **Bertrand J, Comets E, Chenel M, Mentré F** (2012).
  *Some Alternatives to Asymptotic Tests for the Analysis of Pharmacogenetic Data Using Nonlinear Mixed Effects Models.*
  Biometrics, **68**, 146–155.
  https://doi.org/10.1111/j.1541-0420.2011.01665.x

### Sparse designs & related work

- **Loingeville F, Bertrand J, Nguyen TT, Sharan S, Feng K, Sun W, Han J, Grosser S, Zhao L, Fang L, Dette H, Mentré F** (2020).
  *New model-based bioequivalence statistical approaches for pharmacokinetic studies with sparse sampling.*
  The AAPS Journal, **22**, 1–8.

- **Nguyen TT, Bazzoli C, Mentré F** (2011).
  *Design evaluation and optimisation in crossover pharmacokinetic studies analysed by nonlinear mixed effects models.*
  Statistics in Medicine, **30**(22), 2810–2826.
  https://doi.org/10.1002/sim.4390

---

## Authors

**Eden BETE** · University of Lille · eden.bete.etu@univ-lille.fr
**Célia BENGUENNA** · University of Lille · celia.benguenna.etu@univ-lille.fr
---

*Submitted to the Journal of Statistical Software*
