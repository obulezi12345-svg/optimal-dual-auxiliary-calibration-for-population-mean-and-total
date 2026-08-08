# optimal-dual-auxiliary-calibration-for-population-mean-and-total
Optimal Dual Auxiliary Calibration Estimation of Population Mean and Total: R Code for Simulations and Real Data Applications

# Optimal Dual Auxiliary Calibration Estimation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R Language](https://img.shields.io/badge/Language-R%20%3E%3D%204.0-blue.svg)](https://www.r-project.org/)

This repository contains the official **R** code implementation for the simulation studies and real data applications presented in the manuscript:

> **Obulezi, O. J. (2026).** *Optimal Dual Auxiliary Calibration Estimation of Population Mean and Total.* Journal of Official Statistics, 42(1), 89–112.

---

## 📌 Overview

Traditional sampling estimators often underperform when multi-dimensional auxiliary data or complex heteroscedasticity are present. This repository implements an optimal dual-auxiliary calibration approach that minimizes a generalized chi-square distance measure under exact dual benchmark constraints ($x_1$ and $x_2$). 

### Key Features Implemented
1. **Closed-Form Calibration Estimators:** Computation of optimal calibration weights ($w_i^*$) and point estimates for population mean and total.
2. **Analytical & Empirical Variances:** Implementation of sample-based $g$-weight adjusted variance estimators.
3. **Comparative Baseline Estimators:**
   * **Horvitz–Thompson Estimator** (Horvitz & Thompson, 1952)
   * **Classical Univariate Calibration** (Deville & Särndal, 1992)
   * **Bivariate Generalized Calibration** (Gupta et al., 2022)
   * **Robust Dual Calibration** (Singh et al., 2023)
   * **Stratified Dual Calibration** (Abubakar et al., 2024)
   * **Proposed Optimal Dual Calibration** (Obulezi, 2026)
4. **Hypothesis Testing & Power Evaluation:** Simulation scripts evaluating empirical Type I error control ($\alpha = 0.05$) and empirical power curves under Pitman local alternative hypotheses.

---

## 📁 Repository Structure

```text
.
├── R/
│   ├── simulation_study.R      # Monte Carlo simulation framework (Table 2, Fig 1, Fig 2)
│   └── real_data_application.R # Real-world data application on `trees` and `swiss` (Tables 3 & 4)
├── output/                     # Generated CSV results and high-resolution plots
│   ├── simulation_results_performance.csv
│   ├── simulation_results_power.csv
│   ├── plot_relative_efficiency.png
│   └── plot_power_curves.png
├── LICENSE                     # MIT License
└── README.md                   # Repository documentation
