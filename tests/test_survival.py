import numpy as np
import pandas as pd
import pytest

import staging_visualization as sv


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


class TestComputeLogrankStatistics:
    def test_returns_dict_with_expected_keys(self, synthetic_survival_data):
        result = sv.compute_logrank_statistics(
            synthetic_survival_data, min_group_size=15
        )

        assert result is not None
        assert "test_statistic" in result
        assert "p_value" in result
        assert "pairwise_results" in result
        assert "correction_method" in result

    def test_returns_none_for_insufficient_data(self):
        small_df = pd.DataFrame(
            {
                "ajcc_stage": ["Stage I"] * 5,
                "overall_survival_months": [10.0] * 5,
                "survival_event": [True] * 5,
            }
        )

        result = sv.compute_logrank_statistics(small_df, min_group_size=15)

        assert result is None

    def test_returns_none_for_single_stage(self, synthetic_survival_data):
        single_stage = synthetic_survival_data[
            synthetic_survival_data["ajcc_stage"] == "Stage IA"
        ]

        result = sv.compute_logrank_statistics(single_stage, min_group_size=15)

        assert result is None

    def test_returns_none_when_no_valid_rows(self):
        df = pd.DataFrame(
            {
                "ajcc_stage": ["Stage I"],
                "overall_survival_months": [np.nan],
                "survival_event": [True],
            }
        )

        result = sv.compute_logrank_statistics(df, min_group_size=1)

        assert result is None

    def test_pairwise_results_sorted_by_pvalue(self, synthetic_survival_data):
        result = sv.compute_logrank_statistics(
            synthetic_survival_data, min_group_size=15
        )
        assert result is not None
        pairwise = result["pairwise_results"]
        adjusted = [p["p_value_adjusted"] for p in pairwise]

        assert adjusted == sorted(adjusted)
        assert all("p_value" in entry for entry in pairwise)

    def test_test_statistic_is_positive(self, synthetic_survival_data):
        result = sv.compute_logrank_statistics(
            synthetic_survival_data, min_group_size=15
        )

        assert result is not None
        assert result["test_statistic"] > 0

    def test_p_value_in_valid_range(self, synthetic_survival_data):
        result = sv.compute_logrank_statistics(
            synthetic_survival_data, min_group_size=15
        )

        assert result is not None
        assert 0 <= result["p_value"] <= 1


class TestBuildSurvivalSummary:
    def test_returns_empty_for_empty_input(self):
        empty_df = pd.DataFrame(
            columns=[
                "overall_survival_months",
                "ajcc_stage",
                "patient_id",
                "survival_event",
            ]
        )

        result = sv.build_survival_summary(empty_df, ["Stage I"])

        assert result.empty

    def test_returns_expected_columns(self, synthetic_survival_data):
        result = sv.build_survival_summary(
            synthetic_survival_data, ["Stage IA", "Stage IB", "Stage IV"]
        )

        assert "median_os" in result.columns
        assert "median_os_ci_lower" in result.columns
        assert "median_os_ci_upper" in result.columns
        assert "patient_count" in result.columns
        assert "event_rate" in result.columns

    def test_patient_counts_match_input(self, synthetic_survival_data):
        result = sv.build_survival_summary(
            synthetic_survival_data, ["Stage IA", "Stage IB", "Stage IV"]
        )

        for stage in ["Stage IA", "Stage IB", "Stage IV"]:
            expected = len(
                synthetic_survival_data[
                    synthetic_survival_data["ajcc_stage"] == stage
                ]
            )
            assert result.loc[stage, "patient_count"] == expected

    def test_median_ci_contains_true_value_for_exponential_distribution(self):
        np.random.seed(0)
        scale = 24.0
        durations = np.random.exponential(scale, size=500)
        df = pd.DataFrame(
            {
                "patient_id": [f"P{i}" for i in range(len(durations))],
                "ajcc_stage": ["Stage I"] * len(durations),
                "overall_survival_months": durations,
                "survival_event": [True] * len(durations),
            }
        )

        result = sv.build_survival_summary(df, ["Stage I"])
        row = result.loc["Stage I"]
        true_median = np.log(2) * scale

        assert abs(row["median_os"] - true_median) < 2.0
        assert row["median_os_ci_lower"] < true_median < row["median_os_ci_upper"]


class TestComputeRmstByStage:
    def test_rmst_positive_values(self, synthetic_survival_data):
        result = sv.compute_rmst_by_stage(
            synthetic_survival_data,
            horizon_months=60.0,
            stage_order=["Stage IA", "Stage IB", "Stage IV"],
            min_group_size=15,
            n_bootstrap=100,
        )

        assert all(result["rmst_months"] > 0)

    def test_rmst_less_than_horizon(self, synthetic_survival_data):
        horizon = 60.0
        result = sv.compute_rmst_by_stage(
            synthetic_survival_data,
            horizon_months=horizon,
            stage_order=["Stage IA", "Stage IB", "Stage IV"],
            min_group_size=15,
            n_bootstrap=100,
        )

        assert all(result["rmst_months"] <= horizon)

    def test_invalid_horizon_raises_error(self, synthetic_survival_data):
        with pytest.raises(ValueError):
            sv.compute_rmst_by_stage(
                synthetic_survival_data,
                horizon_months=-10.0,
                stage_order=["Stage IA"],
                min_group_size=15,
            )

    def test_empty_dataframe_returns_empty(self):
        empty_df = pd.DataFrame(
            {
                "ajcc_stage": [],
                "overall_survival_months": [],
                "survival_event": [],
            }
        )

        result = sv.compute_rmst_by_stage(
            empty_df, horizon_months=10.0, stage_order=["Stage I"], min_group_size=1
        )

        assert result.empty

    def test_rmst_confidence_intervals_monotonic(self, synthetic_survival_data):
        result = sv.compute_rmst_by_stage(
            synthetic_survival_data,
            horizon_months=60.0,
            stage_order=["Stage IA", "Stage IB", "Stage IV"],
            min_group_size=15,
            n_bootstrap=100,
        )

        assert "rmst_ci_lower" in result.columns
        assert "rmst_ci_upper" in result.columns
        assert all(result["rmst_ci_lower"] <= result["rmst_months"])
        assert all(result["rmst_months"] <= result["rmst_ci_upper"])


def test_bootstrap_rmst_ci_contains_true_value():
    np.random.seed(0)
    scale = 24.0
    durations = pd.Series(np.random.exponential(scale, size=400))
    events = pd.Series([True] * len(durations))
    horizon = 30.0

    point_est, ci_lower, ci_upper = sv.bootstrap_rmst_ci(
        durations, events, horizon=horizon, n_bootstrap=200
    )
    true_rmst = (1 - np.exp(-horizon / scale)) * scale

    assert abs(point_est - true_rmst) < 0.75
    assert ci_lower < true_rmst < ci_upper
