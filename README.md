# glaucoma-surgery-trends

Code and reproducibility materials for the manuscript:

**Ten-Year Trends in Glaucoma Surgical Choices and Patient Characteristics at a Tertiary Hospital in China: a retrospective hospital information system-based cohort study (2015–2024)**

## Authors and affiliations

- **Shiyu Liao** — Department of Ophthalmology, Huaxia Hospital, Chengdu, Sichuan, China
- **Chunying Liu** — Department of Ophthalmology, West China Hospital, Sichuan University, Chengdu, Sichuan, China
- **Xin Wei** — Department of Ophthalmology, West China Hospital, Sichuan University, Chengdu, Sichuan, China

## What this repository contains

This repository provides the public-facing analysis code, codebook, and reproducibility documentation for a retrospective hospital information system (HIS)-based study of glaucoma surgery trends from 2015 to 2024.

The repository is designed for **method transparency and partial reproducibility**. It does **not** include the underlying patient-level clinical data.

## Repository structure

```text
glaucoma-surgery-trends/
├── README.md
├── LICENSE
├── CITATION.cff
├── zenodo.json
├── codebook/
│   └── glaucoma_surgery_codebook_tableS9.xlsx
├── docs/
│   ├── availability_statements.md
│   └── suggested_release_checklist.md
└── R/
    ├── 00_helpers.R
    ├── 01_main_logistic_analysis.R
    ├── 02_sensitivity_nonlinearity_and_interaction.R
    ├── 03_sensitivity_excluding_lens_only.R
    ├── 04_projection_internal_trend.R
    └── 05_clustering_analysis.R
```

## Analysis workflow

The scripts are organized by analysis type rather than by figure number.

1. **00_helpers.R**
   - shared functions
   - package loading
   - date parsing
   - procedure-code standardization
   - cohort construction
   - codebook-based grouping of surgery types

2. **01_main_logistic_analysis.R**
   - builds the analytic cohort
   - fits the primary multivariable logistic regression
   - exports adjusted odds ratios
   - computes discrimination and calibration metrics
   - generates the main forest plot and calibration plot

3. **02_sensitivity_nonlinearity_and_interaction.R**
   - evaluates linearity assumptions for age and calendar year
   - compares linear vs spline and piecewise-year models
   - estimates subtype-by-year interaction effects

4. **03_sensitivity_excluding_lens_only.R**
   - repeats trend analyses after excluding lens-only procedures (`P` and `P+I`)
   - provides an important mechanism-oriented sensitivity analysis

5. **04_projection_internal_trend.R**
   - projects adjusted probabilities of internal outflow surgery through 2030
   - uses patient-level bootstrap confidence intervals

6. **05_clustering_analysis.R**
   - derives patient-level clinical profiles using Gower distance + PAM clustering
   - links profiles back to procedure-level internal/external approach selection
   - includes a GEE sensitivity analysis

## Data availability

The underlying data are not publicly available because they contain sensitive patient information and are subject to institutional and ethical restrictions. De-identified data may be made available to qualified researchers upon reasonable request to the corresponding author, subject to institutional review and any required data-use agreement.

## Code availability

All custom R scripts used for data processing, statistical analysis, and figure generation are provided in this repository. The archived version of this repository is available on Zenodo: https://doi.org/10.5281/zenodo.19043147

## Expected input data

The scripts assume an Excel file containing the original analytic extract with Chinese variable names as used internally in the study database. The public repository does **not** include that file.

Place an authorized copy of the dataset locally and update the `input_file` path at the top of each script before running. The code was intentionally written with **relative paths** and without local machine identifiers.

## Software environment

Recommended environment:

- R >= 4.2
- Packages used across scripts:
  - `readxl`
  - `dplyr`
  - `tidyr`
  - `stringr`
  - `ggplot2`
  - `broom`
  - `writexl`
  - `pROC`
  - `ResourceSelection`
  - `binom`
  - `scales`
  - `splines`
  - `sandwich`
  - `lmtest`
  - `cluster`
  - `geepack`

## How to cite this code

Please cite both:

1. the associated manuscript, once published; and
2. the archived software release used in your work.

Suggested software citation:

Liu C. glaucoma-surgery-trends (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.19043147

## Public release notes

Before making the repository public, confirm that:

- no patient-level data are included;
- no local desktop paths remain in the code;
- no internal server names, credentials, or database locations remain in the code;
- no intermediate files containing potentially identifiable data are included;
- all repository links and DOI placeholders have been updated.
