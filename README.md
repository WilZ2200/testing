# Breast Cancer Screening Preferences: A Discrete Choice Experiment

## Study Overview

This project implements a Discrete Choice Experiment (DCE) to investigate women's preferences for breast cancer screening programs. The study examines trade-offs between seven key screening attributes: method, frequency, cost, sensitivity, false positive rate, waiting time, and pain/discomfort.

**Protocol:** DCE-BCS-2026-001
**Status:** In development
**Target population:** Women aged 40-74 years, United States

## Research Objectives

1. Estimate the relative importance of breast cancer screening attributes in women's preferences
2. Calculate willingness-to-pay (WTP) for improvements in each attribute
3. Identify preference heterogeneity across demographic and clinical subgroups
4. Estimate predicted uptake rates for different screening program configurations

## DCE Attributes

| Attribute | Levels |
|-----------|--------|
| Screening method | Mammography, MRI, Ultrasound |
| Screening frequency | Every year, Every 2 years, Every 3 years |
| Out-of-pocket cost | $0, $50, $150, $300 |
| Sensitivity | 70%, 85%, 95% |
| False positive rate | 5%, 10%, 15% |
| Waiting time for results | 1 day, 1 week, 3 weeks |
| Pain/discomfort | None, Mild, Moderate |

## Directory Structure

```
project/
├── README.md                  # This file
├── .gitignore                 # Git ignore rules
├── survey/                    # Survey instruments and design
│   ├── dce_questionnaire.tex  # Main DCE survey (LaTeX)
│   ├── experimental_design.md # D-optimal design documentation
│   └── distribution/          # Distribution and ethics
│       ├── distribution_plan.md    # Recruitment and distribution strategy
│       ├── informed_consent.tex    # IRB informed consent form
│       └── irb_outline.md         # IRB application outline
├── data/                      # Study data
│   ├── raw/                   # Raw survey exports (do not edit)
│   └── processed/             # Cleaned and analysis-ready data
├── analysis/                  # R scripts for statistical analysis
├── manuscript/                # Manuscript drafts
│   └── sections/              # Individual manuscript sections
├── figures/                   # Generated figures and plots
├── tables/                    # Generated tables
└── literature/                # Reference papers and literature notes
```

## Methodology

- **Design:** D-optimal fractional factorial, generated with Ngene
- **Choice sets:** 48 total, blocked into 4 versions of 12 tasks each
- **Alternatives:** 2 labeled alternatives + opt-out per choice set
- **Sample size:** Target n = 600 (after quality exclusions)
- **Analysis:** Mixed logit (primary), conditional logit, latent class models

## Software Requirements

- **R** (>= 4.2.0) with packages: `apollo`, `mlogit`, `survival`, `ggplot2`, `dplyr`
- **LaTeX** (TeX Live or MikTeX) for survey compilation
- **Ngene** (>= 1.3) for experimental design generation
- **Qualtrics** for survey hosting and distribution

## Timeline

| Phase | Timeline |
|-------|----------|
| Design and development | March-April 2026 |
| IRB approval | April-May 2026 |
| Pilot testing | May-June 2026 |
| Data collection | July-August 2026 |
| Analysis | September-November 2026 |
| Manuscript | November 2026-January 2027 |

## Team

- **PI:** [Name]
- **Co-Investigators:** [Names]
- **Research Coordinator:** [Name]

## License

This project is for academic research purposes. Data and analysis code will be made available upon publication per journal data-sharing policies.
