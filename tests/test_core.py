import copy

import pandas as pd
import pytest

import staging_visualization as sv


@pytest.fixture
def minimal_config():
    config = copy.deepcopy(sv.DEFAULT_CONFIG)
    config["reverse_column_mapping"] = {
        value: key for key, value in config["column_mapping"].items()
    }
    return config


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("Stage IIIA", "Stage IIIA"),
        ("stage 3a", "Stage 3A"),
        ("NA", None),
    ],
)
def test_clean_stage_label_handles_common_inputs(raw, expected):
    assert sv.clean_stage_label(raw) == expected


def test_preprocess_missing_required_columns_raises_keyerror(minimal_config):
    df = pd.DataFrame({"Patient ID": ["TCGA-01"]})

    with pytest.raises(KeyError):
        sv.preprocess(df, minimal_config)


def test_rmst_calculation():
    df = pd.DataFrame(
        {
            "ajcc_stage": ["Stage I", "Stage I"],
            "overall_survival_months": [10.0, 10.0],
            "survival_event": [True, True],
        }
    )

    rmst = sv.compute_rmst_by_stage(
        df,
        horizon_months=20.0,
        stage_order=["Stage I"],
        min_group_size=1,
        n_bootstrap=50,
    )

    assert not rmst.empty
    assert rmst["rmst_months"].iloc[0] == 10.0
    assert rmst["rmst_ci_lower"].iloc[0] == 10.0
    assert rmst["rmst_ci_upper"].iloc[0] == 10.0
