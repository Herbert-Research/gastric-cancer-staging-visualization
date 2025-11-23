import matplotlib
import numpy as np
import pandas as pd
import pytest

import staging_visualization as sv

# Ensure plotting uses a non-interactive backend in test environments.
matplotlib.use("Agg")


@pytest.fixture
def synthetic_survival_data():
    np.random.seed(42)
    n_per_stage = 30
    stages = ["Stage IA", "Stage IB", "Stage IV"]
    data: list[dict[str, object]] = []
    for stage in stages:
        for i in range(n_per_stage):
            mean_os = 10.0 if stage == "Stage IV" else 30.0
            data.append(
                {
                    "patient_id": f"{stage}_{i}",
                    "ajcc_stage": stage,
                    "overall_survival_months": max(1, np.random.exponential(mean_os)),
                    "survival_event": np.random.choice([True, False], p=[0.6, 0.4]),
                }
            )
    return pd.DataFrame(data)


class TestPlotKaplanMeier:
    def test_returns_false_for_empty_dataframe(self, tmp_path):
        empty_df = pd.DataFrame(
            columns=["overall_survival_months", "ajcc_stage", "survival_event"]
        )
        output_path = tmp_path / "km_test.png"

        result = sv.plot_kaplan_meier(
            empty_df,
            output_path,
            stage_order=["Stage I", "Stage II"],
            min_group_size=15,
            figure_dpi=100,
        )

        assert result is False
        assert not output_path.exists()

    def test_creates_file_for_valid_data(self, synthetic_survival_data, tmp_path):
        output_path = tmp_path / "km_test.png"

        result = sv.plot_kaplan_meier(
            synthetic_survival_data,
            output_path,
            stage_order=["Stage IA", "Stage IB", "Stage IV"],
            min_group_size=15,
            figure_dpi=100,
        )

        assert result is True
        assert output_path.exists()
        assert output_path.stat().st_size > 0

    def test_excludes_small_groups(self, tmp_path):
        data = pd.DataFrame(
            {
                "ajcc_stage": ["Stage I"] * 5 + ["Stage II"] * 30,
                "overall_survival_months": [10.0] * 35,
                "survival_event": [True] * 35,
            }
        )
        output_path = tmp_path / "km_test.png"

        result = sv.plot_kaplan_meier(
            data,
            output_path,
            stage_order=["Stage I", "Stage II"],
            min_group_size=15,
            figure_dpi=100,
        )

        assert result is True


class TestPlotSurvivalByStage:
    def test_returns_false_for_empty_summary(self, tmp_path):
        empty_summary = pd.DataFrame()
        output_path = tmp_path / "survival_test.png"

        result = sv.plot_survival_by_stage(
            empty_summary, output_path, min_group_size=15, figure_dpi=100
        )

        assert result is False

    def test_creates_file_for_valid_summary(self, tmp_path):
        summary = pd.DataFrame(
            {
                "median_os": [10.0, 12.0],
                "patient_count": [5, 6],
                "event_rate": [0.4, 0.5],
                "event_rate_pct": [40.0, 50.0],
            },
            index=["Stage I", "Stage II"],
        )
        output_path = tmp_path / "survival_stage.png"

        result = sv.plot_survival_by_stage(
            summary, output_path, min_group_size=2, figure_dpi=100
        )

        assert result is True
        assert output_path.exists()

    def test_excludes_small_groups_and_logs(self, tmp_path, capsys):
        summary = pd.DataFrame(
            {
                "median_os": [10.0],
                "patient_count": [1],
                "event_rate": [0.2],
                "event_rate_pct": [20.0],
            },
            index=["Stage I"],
        )
        output_path = tmp_path / "survival_excluded.png"

        result = sv.plot_survival_by_stage(
            summary, output_path, min_group_size=5, figure_dpi=100
        )

        captured = capsys.readouterr()
        assert result is False
        assert "no stages met min_group_size" in captured.out


class TestPlotTnHeatmap:
    def test_returns_false_for_missing_tn_data(self, tmp_path):
        df = pd.DataFrame({"t_stage": [None, None], "n_stage": [None, None]})
        output_path = tmp_path / "heatmap_test.png"

        result = sv.plot_tn_heatmap(
            df,
            output_path,
            t_stage_order=["T1", "T2"],
            n_stage_order=["N0", "N1"],
            figure_dpi=100,
        )

        assert result is False

    def test_creates_file_for_valid_tn_data(self, tmp_path):
        df = pd.DataFrame(
            {"t_stage": ["T1", "T2", "T3"] * 20, "n_stage": ["N0", "N1", "N2"] * 20}
        )
        output_path = tmp_path / "heatmap_test.png"

        result = sv.plot_tn_heatmap(
            df,
            output_path,
            t_stage_order=["T1", "T2", "T3"],
            n_stage_order=["N0", "N1", "N2"],
            figure_dpi=100,
        )

        assert result is True
        assert output_path.exists()


class TestPlotStagePanels:
    def test_creates_stage_panel_figure(self, tmp_path):
        t_counts = pd.Series([3, 2], index=["T1", "T2"])
        n_counts = pd.Series([4, 1], index=["N0", "N1"])
        m_counts = pd.Series([5, 0], index=["M0", "M1"])
        stage_counts = pd.Series([3, 2], index=["Stage I", "Stage II"])
        output_path = tmp_path / "panels.png"

        result = sv.plot_stage_panels(
            t_counts, n_counts, m_counts, stage_counts, output_path, figure_dpi=100
        )

        assert result is True
        assert output_path.exists()
