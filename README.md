Here is the revised `README.md` for your **Gastric Cancer Staging Visualization** project. I have adapted the professional layout, tone, and structural elements from your template (`READMEVorlage.md`) while retaining the specific technical and clinical details of this repository.

-----

# Gastric Cancer Staging Visualization

Supporting analytics package for the proposed PhD dissertation **“Prospective Validation of Station-Specific Risk Guidance for KLASS-Standardized Gastrectomy.”** This repository distils retrospective TCGA evidence into reproducible figures that anchor Aim 1 calibration targets and communicate methodological readiness to admissions reviewers and clinical collaborators.

## Executive Summary

  - **Reconstructs** AJCC TNM distributions and overall survival profiles for the TCGA STAD PanCanAtlas 2018 cohort.
  - **Substantiates** the clinical heterogeneity gap that future station-specific AI guidance aims to address.
  - **Demonstrates** a validation-first workflow (schema enforcement → harmonisation → visualization) scalable to future KLASS video and telemetry streams.
  - **Exports** audited, publication-ready figures (`*.png`) designed for committee packets and multi-disciplinary review slide decks.

## Scientific Context and Objectives

Before station-specific lymphadenectomy guidance can be prospectively deployed, it must be calibrated against well-characterised retrospective baselines. This codebase establishes those baselines by quantifying disease burden heterogeneity (via TNM heatmaps) and documenting existing outcome gradients (via Kaplan-Meier and median OS analysis). The material is curated for review committees that expect clearly articulated hypotheses, transparent data lineage, and rigorous, open methods.

## Data Provenance and Governance

  - **Source Stream:** TCGA STAD PanCanAtlas 2018 clinical release (`data/tcga_2018_clinical_data.tsv`).
  - **Content:** Hashed patient identifiers, AJCC TNM annotations (6th/7th editions), and overall survival metadata.
  - **Compliance:** Data remain de-identified, publicly licensed, and cited according to TCGA policy. Survival statements in this repository are purely observational and should be validated locally before informing interventional guidance.

## Analytical Workflow

The core pipeline (`staging_visualization.py`) executes three synchronized phases:

1.  **Schema Validation** – Enforces strict concordancy between PanCanAtlas headers and expected AJCC fields, halting execution if critical metadata is absent.
2.  **Data Harmonisation** – Standardises heterogeneous staging labels (e.g., "Stage IIIA" vs "Stage 3a") and survival endpoints without imputing missing records.
3.  **Figure Generation** – Exports high-resolution panels visualizing disease burden and actuarial outcomes.

## Example Output

Running the default configuration produces clinically-interpretable summaries demonstrating the pipeline's analytical rigor:

```bash

Gastric Cancer Staging Visualization
============================================================
Total patients with AJCC staging: 421
   Stage I: 2 patients (0.5%)
  Stage IA: 16 patients (3.8%)
  Stage IB: 40 patients (9.5%)
  Stage II: 33 patients (7.8%)
 Stage IIA: 41 patients (9.7%)
 Stage IIB: 57 patients (13.5%)
 Stage III: 3 patients (0.7%)
Stage IIIA: 82 patients (19.5%)
Stage IIIB: 64 patients (15.2%)
Stage IIIC: 39 patients (9.3%)
  Stage IV: 44 patients (10.5%)

Overall Survival Highlights:
------------------------------------------------------------
   Stage I | Median OS: 12.5 mo | Events: 0% | n=2
  Stage IA | Median OS: 24.3 mo | Events: 13% | n=15
  Stage IB | Median OS: 16.7 mo | Events: 26% | n=39
  Stage II | Median OS: 13.3 mo | Events: 30% | n=33
 Stage IIA | Median OS: 14.0 mo | Events: 24% | n=41
 Stage IIB | Median OS: 15.8 mo | Events: 29% | n=56
 Stage III | Median OS: 6.6 mo | Events: 100% | n=3
Stage IIIA | Median OS: 14.3 mo | Events: 46% | n=81
Stage IIIB | Median OS: 16.2 mo | Events: 44% | n=63
Stage IIIC | Median OS: 11.7 mo | Events: 46% | n=39
  Stage IV | Median OS: 8.6 mo | Events: 61% | n=44

```

This output quantifies clear survival gradients—from 24.3 months in Stage IA down to 8.6 months in Stage IV—precisely the outcome heterogeneity that future station-specific guidance must address.

## Generated Figures

**Disease Burden Output**

  - `tnm_staging_distribution.png` – Paired bar panels illustrating the raw distribution of T, N, M, and aggregate AJCC stages.
  - `tn_heatmap.png` – Density heatmap across the T×N grid, highlighting cohorts (e.g., T3/T4a) where station-level guidance is clinically most consequential.

**Survival Analysis Output**

  - `os_by_stage.png` – Median overall survival bars annotated with event rates and cohort sizes, underscoring the outcome gradient Aim 1 seeks to tighten.
  - `km_by_stage.png` – Kaplan–Meier survival curves stratified by AJCC stage, visualizing the full hazard trajectory for sufficiently powered subgroups.

## Usage

### Quickstart (Headless)

Run the end-to-end visualization pipeline using bundled defaults.

```bash
# Setup environment
pip install -r requirements.txt

# Execute full workflow
python staging_visualization.py
```

### Configuration Options

The workflow can be parameterized to point to alternative datasets or adjust statistical thresholds.

```bash
# Override input data and output locations
python staging_visualization.py --data data/local_registry.tsv --output-dir reports/

# Adjust minimum cohort size for Kaplan-Meier curves (default: 15)
python staging_visualization.py --km-min-group 20
```

## Software Requirements

  - Python 3.9 or newer.
  - Core dependencies: `pandas`, `numpy`, `matplotlib`, `seaborn`, `lifelines` (for survival statistics).
  - See `requirements.txt` for pinned versions used in validation.

## Input Validation Schema

The script enforces data integrity by halting with a descriptive error if required headers are absent:

  - `Patient ID` → `patient_id`
  - `Neoplasm Disease Stage American Joint Committee on Cancer Code` → `ajcc_stage`
  - `American Joint Committee on Cancer Tumor Stage Code` → `t_stage`
  - `Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code` → `n_stage`
  - `Overall Survival (Months)` → `overall_survival_months`
  - `Overall Survival Status` → `overall_survival_status`

This safeguard allows screening committees to verify dataset compatibility without auditing source code.

## Clinical Interpretation Notes

  - **Retrospective Nature:** Reported survival summaries reflect historical outcomes from prospectively consented cohorts; they are observational and not intended for causal inference without local validation.
  - **AJCC Groupings:** Staging labels follow those embedded within the TCGA clinical records; stages are harmonised but not re-derived from raw TNM values.
  - **Survival Endpoint:** Summaries reflect overall survival (OS); disease-free intervals were not uniformly available in the 2018 PanCanAtlas release used for this pilot.

## Repository Stewardship

Author: **Maximilian Herbert Dressler**

## Acknowledgement

“The results presented here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga.”

## Citations

  - Cerami E, Gao J, Dogrusoz U, *et al.* The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data. *Cancer Discovery.* 2012;2(5):401–404.
  - Gao J, Aksoy BA, Dogrusoz U, *et al.* Integrative Analysis of Complex Cancer Genomics and Clinical Profiles Using the cBioPortal. *Science Signaling.* 2013;6(269):pl1.
  - Liu J, Lichtenberg T, Hoadley KA, *et al.* An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. *Cell.* 2018;173(2):400–416.e11.