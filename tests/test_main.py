import copy
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml

import staging_visualization as sv


def test_print_summary_handles_empty_stage_counts(capsys):
    stage_counts = pd.Series(dtype=int)

    text = sv.print_summary(stage_counts, pd.DataFrame())

    assert "No AJCC staging data available" in text
    assert "Gastric Cancer Staging Visualization" in text


def test_print_summary_with_logrank_and_rmst(capsys):
    stage_counts = pd.Series([2, 3], index=["Stage I", "Stage II"])
    survival_summary = pd.DataFrame(
        {
            "median_os": [10.0, 12.0],
            "median_os_ci_lower": [9.0, 11.0],
            "median_os_ci_upper": [11.0, 13.0],
            "patient_count": [2, 3],
            "event_rate": [0.5, 0.66],
            "event_rate_pct": [50.0, 66.0],
        },
        index=["Stage I", "Stage II"],
    )
    logrank_results = {
        "test_statistic": 1.23,
        "p_value": 0.0456,
        "degrees_freedom": 1,
        "correction_method": "Benjamini-Hochberg FDR",
        "pairwise_results": [
            {
                "stages": ("Stage I", "Stage II"),
                "test_statistic": 1.1,
                "p_value": 0.05,
                "p_value_adjusted": 0.06,
            }
        ],
    }
    rmst_summary = pd.DataFrame(
        {
            "ajcc_stage": ["Stage I"],
            "rmst_months": [9.5],
            "rmst_ci_lower": [8.5],
            "rmst_ci_upper": [10.5],
            "events": [1],
            "patient_count": [2],
        }
    )

    text = sv.print_summary(
        stage_counts,
        survival_summary,
        logrank_results=logrank_results,
        rmst_summary=rmst_summary,
        rmst_horizon=12.0,
    )

    assert "Log-Rank Test (Omnibus):" in text
    assert "Restricted Mean Survival Time" in text
    assert "Stage I vs Stage II" in text
    assert "adj. p=" in text
    assert "95% CI" in text


def test_main_end_to_end(monkeypatch, tmp_path):
    data = pd.DataFrame(
        {
            "Patient ID": ["P1", "P2", "P3", "P4"],
            "Neoplasm Disease Stage American Joint Committee on Cancer Code": [
                "Stage I",
                "Stage I",
                "Stage II",
                "Stage II",
            ],
            "American Joint Committee on Cancer Tumor Stage Code": ["T1", "T1", "T2", "T2"],
            "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": [
                "N0",
                "N1",
                "N0",
                "N1",
            ],
            "American Joint Committee on Cancer Metastasis Stage Code": ["M0", "M0", "M0", "M0"],
            "Overall Survival (Months)": [12, 10, 6, 8],
            "Overall Survival Status": ["0:LIVING", "1:DECEASED", "1:DECEASED", "0:LIVING"],
        }
    )
    data_path = tmp_path / "cohort.tsv"
    data.to_csv(data_path, sep="\t", index=False)

    config = copy.deepcopy(sv.DEFAULT_CONFIG)
    output_dir = tmp_path / "outputs"
    config["file_settings"]["default_data_path"] = str(data_path)
    config["file_settings"]["default_output_dir"] = str(output_dir)
    config["file_settings"]["figure_dpi"] = 100
    config["visualization"]["km_min_group_size"] = 1
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    monkeypatch.setattr(sys, "argv", ["prog", "--config", str(config_path), "--rmst-months", "12"])

    sv.main()

    figure_names = config["visualization"]["filenames"]
    expected_files = [
        output_dir / figure_names["stage_distribution"],
        output_dir / figure_names["tn_heatmap"],
        output_dir / figure_names["survival_by_stage"],
        output_dir / figure_names["km_by_stage"],
        output_dir / figure_names["rmst_by_stage"],
        output_dir / figure_names["statistical_summary"],
    ]
    for path in expected_files:
        assert path.exists()
        assert path.stat().st_size > 0
