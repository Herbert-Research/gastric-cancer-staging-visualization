import copy

import numpy as np
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


class TestCleanStageLabel:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("Stage IIIA", "Stage IIIA"),
            ("stage 3a", "Stage 3A"),
            ("STAGE IB", "Stage IB"),
            ("  Stage II  ", "Stage II"),
            ("NA", None),
            ("N/A", None),
            ("", None),
            (None, None),
            (np.nan, None),
        ],
    )
    def test_handles_common_inputs(self, raw, expected):
        assert sv.clean_stage_label(raw) == expected


class TestNormalizeTnmCode:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("T1", "T1"),
            ("t2a", "T2A"),
            ("  N0  ", "N0"),
            ("NA", None),
            ("", None),
            (None, None),
        ],
    )
    def test_handles_common_inputs(self, raw, expected):
        assert sv.normalize_tnm_code(raw) == expected


class TestIsEvent:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("1:DECEASED", True),
            ("0:LIVING", False),
            ("DECEASED", True),
            ("DEAD", True),
            ("Alive", False),
            ("", False),
            (None, False),
        ],
    )
    def test_handles_common_inputs(self, raw, expected):
        assert sv.is_event(raw) == expected


class TestPreprocess:
    def test_missing_required_columns_raises_keyerror(self, minimal_config):
        df = pd.DataFrame({"Patient ID": ["TCGA-01"]})

        with pytest.raises(KeyError):
            sv.preprocess(df, minimal_config)

    def test_valid_input_returns_dataframe(self, minimal_config):
        df = pd.DataFrame(
            {
                "Patient ID": ["TCGA-01"],
                "Neoplasm Disease Stage American Joint Committee on Cancer Code": [
                    "Stage IIA"
                ],
                "American Joint Committee on Cancer Tumor Stage Code": ["T2"],
                "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": [
                    "N1"
                ],
                "American Joint Committee on Cancer Metastasis Stage Code": ["M0"],
                "Overall Survival (Months)": ["24.5"],
                "Overall Survival Status": ["0:LIVING"],
            }
        )

        result = sv.preprocess(df, minimal_config)

        assert isinstance(result, pd.DataFrame)
        assert "ajcc_stage" in result.columns
        assert "survival_event" in result.columns

    def test_survival_months_converted_to_numeric(self, minimal_config):
        df = pd.DataFrame(
            {
                "Patient ID": ["TCGA-01"],
                "Neoplasm Disease Stage American Joint Committee on Cancer Code": [
                    "Stage I"
                ],
                "American Joint Committee on Cancer Tumor Stage Code": ["T1"],
                "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": [
                    "N0"
                ],
                "American Joint Committee on Cancer Metastasis Stage Code": ["M0"],
                "Overall Survival (Months)": ["24.5"],
                "Overall Survival Status": ["0:LIVING"],
            }
        )

        result = sv.preprocess(df, minimal_config)

        assert result["overall_survival_months"].dtype == np.float64


class TestUtilityFunctions:
    def test_extract_status_text(self):
        assert sv.extract_status_text("0:LIVING") == "Living"
        assert sv.extract_status_text("Deceased") == "Deceased"
        assert sv.extract_status_text(None) is None

    def test_ordered_counts_and_category_order(self):
        series = pd.Series(["Stage II", "Stage I", "Stage II", None])
        order = ["Stage I", "Stage II", "Stage III"]
        counts = sv.ordered_counts(series, order)
        assert list(counts.index)[:2] == ["Stage I", "Stage II"]

        categories = sv.determine_category_order(pd.Series(["N1", "N2", "NX"]), ["N0", "N1"])
        assert categories == ["N1", "N2", "NX"]


class TestSchemaValidation:
    def test_negative_survival_months_raises_error(self, minimal_config):
        df = pd.DataFrame({
            "Patient ID": ["TCGA-01"],
            "Neoplasm Disease Stage American Joint Committee on Cancer Code": ["Stage I"],
            "American Joint Committee on Cancer Tumor Stage Code": ["T1"],
            "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": ["N0"],
            "American Joint Committee on Cancer Metastasis Stage Code": ["M0"],
            "Overall Survival (Months)": ["-5.0"],  # Invalid
            "Overall Survival Status": ["0:LIVING"]
        })
        with pytest.raises(ValueError, match="schema validation"):
            sv.preprocess(df, minimal_config)

    def test_duplicate_patient_ids_raises_error(self, minimal_config):
        df = pd.DataFrame({
            "Patient ID": ["TCGA-01", "TCGA-01"],  # Duplicate
            "Neoplasm Disease Stage American Joint Committee on Cancer Code": ["Stage I", "Stage II"],
            "American Joint Committee on Cancer Tumor Stage Code": ["T1", "T2"],
            "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": ["N0", "N1"],
            "American Joint Committee on Cancer Metastasis Stage Code": ["M0", "M0"],
            "Overall Survival (Months)": ["24.0", "36.0"],
            "Overall Survival Status": ["0:LIVING", "1:DECEASED"]
        })
        with pytest.raises(ValueError, match="schema validation"):
            sv.preprocess(df, minimal_config)

    def test_invalid_stage_code_raises_error(self, minimal_config):
        df = pd.DataFrame({
            "Patient ID": ["TCGA-01"],
            "Neoplasm Disease Stage American Joint Committee on Cancer Code": ["Stage V"],  # Invalid
            "American Joint Committee on Cancer Tumor Stage Code": ["T1"],
            "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": ["N0"],
            "American Joint Committee on Cancer Metastasis Stage Code": ["M0"],
            "Overall Survival (Months)": ["24.0"],
            "Overall Survival Status": ["0:LIVING"]
        })
        with pytest.raises(ValueError, match="schema validation"):
            sv.preprocess(df, minimal_config)
