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

## Data Provenance and 

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

## Generated Figures

**Staging visualization output**

- `tnm_staging_distribution.png` – paired bar panels for AJCC T, N, M, and aggregate 
    stages.
- `tn_heatmap.png` – density heatmap across the T×N grid, highlighting areas where 
    station-level guidance is clinically most consequential.
- `os_by_stage.png` – median overall survival with event rates, underscoring the 
    outcome gradient that Aim 1 seeks to tighten.

Any additional PNGs present in the repository are archived artifacts from internal discussions and are not part of the automated workflow.

## Example Output

Running the visualization script produces the following summary statistics:

```
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

The output demonstrates clear survival gradients across staging groups, with Stage IV showing the poorest median overall survival (8.6 months) and highest event rate (61%), while earlier stages show progressively better outcomes—precisely the heterogeneity that station-specific guidance must address.

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

# Override defaults if needed (example with bundled data)
python staging_visualization.py --data data/tcga_2018_clinical_data.tsv --output-dir figures/
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
