# Gastric Cancer Staging Visualization

Supporting analytics package for the proposed PhD dissertation **“Prospective Validation of Station-Specific Risk Guidance for KLASS-Standardized Gastrectomy.”** The repository distils retrospective TCGA evidence into reproducible figures that anchor Aim 1 calibration targets and communicate methodological readiness to admissions reviewers and clinical collaborators.

## Executive Summary
- Reconstructs AJCC TNM distributions and overall survival for the TCGA STAD 
    PanCanAtlas 2018 cohort.
- Demonstrates a validation-first workflow that will later be extended to KLASS video, 
    fluorescence, and telemetry data streams.
- Ships with audited figures (`*.png`) and open-source code so faculty can 
    independently reproduce every analytic claim in the proposal dossier.

## Scientific Context and Objectives
This codebase substantiates the translational premise that station-specific lymphadenectomy guidance must be calibrated against well-characterised baselines before prospective deployment. The TNM visualizations quantify disease burden heterogeneity, while survival panels document outcome gradients that KLASS-referenced metrics must ultimately improve upon. The material is curated for review committees that expect clearly articulated hypotheses, transparent data lineage, and rigorous, open methods.

## Data Provenance and Governance
- **Source:** TCGA STAD PanCanAtlas 2018 clinical release
    (`data\tcga_2018_clinical_data.tsv`).
- **Content:** Patient identifiers (hashed), AJCC TNM annotations (6th/7th editions) 
    and overall survival metadata.
- **Compliance:** Data remain de-identified, publicly licensed, and cited according to 
    TCGA policy. Survival statements in this repository are purely observational and should be validated locally before informing interventional guidance.

## Analytical Workflow
`staging_visualization.py` performs three phases:
1. **Schema validation** – verifies that PanCanAtlas headers match the expected 
    AJCC-concordant fields listed below.
2. **Data harmonisation** – standardises staging labels and survival endpoints without 
    imputing or simulating records.
3. **Figure generation** – exports publication-quality PNGs for committee packets and 
    slide decks.

`risk_calculator.py` extends the repo with a KLASS-inspired risk modelling workflow:
1. **Model configuration** – loads `models/heuristic_klass.json` (logistic recurrence 
    heuristics) and, when available, `models/han2012_jco.json` (Han 2012 Cox survival nomogram).
2. **Cohort harmonisation** – reuses the TCGA TSV, imputing tumor size and LN ratio 
    when clinical fields are missing and flagging every substitution.
3. **Dual-model reporting** – exports recurrence, sensitivity, survival, and    
    calibration figures so faculty can evaluate actuarial readiness alongside staging visuals.

## Generated Figures
**Staging visualization output**
- `tnm_staging_distribution.png` – paired bar panels for AJCC T, N, M, and aggregate 
    stages.
- `tn_heatmap.png` – density heatmap across the T×N grid, highlighting areas where 
    station-level guidance is clinically most consequential.
- `os_by_stage.png` – median overall survival with event rates, underscoring the 
    outcome gradient that Aim 1 seeks to tighten.

**Risk calculator output**
- `risk_predictions.png` – patient-level recurrence estimates with categorical mix.
- `sensitivity_analysis.png` – tornado-style chart showing how LN yield shifts risk.
- `tcga_cohort_summary.png` – histogram + TN heatmap for cohort-wide predictions.
- `calibration_curve.png` – optional if scikit-learn is installed and event labels    
    exist.
- `survival_predictions_han2012.png` and `survival_vs_recurrence_comparison.png` – 
    only when the Han 2012 Cox config is present and `--skip-survival` is not set.

Any additional PNGs present in the repository are archived artifacts from internal discussions and are not part of the automated workflow.

## Usage

### Staging Visualization

```bash
# (Optional) prepare a virtual environment
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Produce the figures using the bundled TSV and repository root as output
python staging_visualization.py

# Override defaults if needed
python staging_visualization.py --data /path/to/tcga.tsv --output-dir figures/
```

### Risk Calculator

```bash
# Run the dual-model risk workflow against the bundled TCGA cohort
python risk_calculator.py --data data/tcga_2018_clinical_data.tsv --output-dir reports/

# Point to custom configs or disable survival outputs if Cox resources are unavailable
python risk_calculator.py \
  --model-config models/heuristic_klass.json \
  --survival-model models/han2012_jco.json \
  --skip-survival
```

Default inputs live in `data/tcga_2018_clinical_data.tsv`, but both the dataset and configuration paths can be swapped to reflect local registries or institution-specific coefficients.

## Software Requirements
- Python 3.9 or newer.
- Packages enumerated in `requirements.txt` (pandas, numpy, matplotlib, seaborn).
- For Debian/Ubuntu hosts lacking tooling, install once:

```bash
sudo apt update && sudo apt install -y python3-pip python3-venv
python3 -m pip install --upgrade pip
```

## Input Validation Schema
The script halts with a descriptive message if any required PanCanAtlas header is absent:

- `Patient ID` → `patient_id`
- `Neoplasm Disease Stage American Joint Committee on Cancer Code` → `ajcc_stage`
- `American Joint Committee on Cancer Tumor Stage Code` → `t_stage`
- `Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code` → 
    `n_stage`
- `American Joint Committee on Cancer Metastasis Stage Code` → `m_stage`
- `Overall Survival (Months)` → `overall_survival_months`
- `Overall Survival Status` → `overall_survival_status`

This safeguard enables screening committees to verify dataset integrity without digging into the source code.

## Clinical Interpretation Notes
- AJCC groupings follow the labels embedded within TCGA; stages are not re-derived.
- Although the cohort was prospectively consented, it remains retrospective 
    observational data and should not be used for causal inference without institutional validation.
- Reported survival summaries reflect overall survival only; disease-free 
    intervals were not uniformly available in the 2018 release.

## Repository Stewardship
Author: **Maximilian Herbert Dressler**

## Acknowledgement
“The results presented here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga.”

## Citations
- Cerami E, Gao J, Dogrusoz U, *et al.* The cBio Cancer Genomics Portal: An Open  
    Platform for Exploring Multidimensional Cancer Genomics Data. *Cancer Discovery.* 2012;2(5):401–404.
- Gao J, Aksoy BA, Dogrusoz U, *et al.* Integrative Analysis of Complex Cancer 
    Genomics and Clinical Profiles Using the cBioPortal. *Science Signaling.* 2013;6(269):pl1.
- Liu J, Lichtenberg T, Hoadley KA, *et al.* An Integrated TCGA Pan-Cancer Clinical 
    Data Resource to Drive High-Quality Survival Outcome Analytics. *Cell.* 2018;173(2):400–416.e11.
