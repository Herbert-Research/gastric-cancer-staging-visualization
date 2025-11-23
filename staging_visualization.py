"""
Gastric cancer staging visualizer backed by TCGA STAD clinical data.

The script ingests the publicly available TCGA PanCanAtlas 2018 gastric cancer
clinical TSV, harmonises AJCC TNM annotations, and outputs figures that can
support KLASS-standardised workflow discussions.
"""

from __future__ import annotations

import argparse
import copy
import os
import sys
import tempfile

# Configure matplotlib cache directory to suppress warnings
os.environ['MPLCONFIGDIR'] = os.path.join(tempfile.gettempdir(), '.matplotlib_cache')
from pathlib import Path
from typing import Any, Iterable, Optional, TypedDict, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import pandera.pandas as pa
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test, pairwise_logrank_test
from lifelines.utils import median_survival_times, restricted_mean_survival_time
from statsmodels.stats.multitest import multipletests
import yaml
from schemas import PreprocessedCohortSchema

sns.set_style("whitegrid")

BASE_DIR = Path(__file__).resolve().parent
CONFIG_PATH_DEFAULT = BASE_DIR / "config.yaml"

DEFAULT_CONFIG = {
    "file_settings": {
        "default_data_path": "data/tcga_2018_clinical_data.tsv",
        "default_output_dir": ".",
        "figure_dpi": 300,
    },
    "column_mapping": {
        "Patient ID": "patient_id",
        "Sample ID": "sample_id",
        "Neoplasm Disease Stage American Joint Committee on Cancer Code": "ajcc_stage",
        "American Joint Committee on Cancer Tumor Stage Code": "t_stage",
        "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": "n_stage",
        "American Joint Committee on Cancer Metastasis Stage Code": "m_stage",
        "Overall Survival (Months)": "overall_survival_months",
        "Overall Survival Status": "overall_survival_status",
    },
    "required_columns": (
        "patient_id",
        "ajcc_stage",
        "t_stage",
        "n_stage",
        "m_stage",
        "overall_survival_months",
        "overall_survival_status",
    ),
    "stage_ordering": {
        "ajcc_stages": [
            "Stage 0",
            "Stage I",
            "Stage IA",
            "Stage IB",
            "Stage II",
            "Stage IIA",
            "Stage IIB",
            "Stage III",
            "Stage IIIA",
            "Stage IIIB",
            "Stage IIIC",
            "Stage IVA",
            "Stage IVB",
            "Stage IV",
        ],
        "t_stages": [
            "Tis",
            "T1",
            "T1A",
            "T1B",
            "T2",
            "T2A",
            "T2B",
            "T3",
            "T4",
            "T4A",
            "T4B",
            "TX",
        ],
        "n_stages": ["N0", "N1", "N2", "N3", "N3A", "N3B", "NX"],
        "m_stages": ["M0", "M1", "MX"],
    },
    "visualization": {
        "km_min_group_size": 15,
        "filenames": {
            "stage_distribution": "tnm_staging_distribution.png",
            "tn_heatmap": "tn_heatmap.png",
            "survival_by_stage": "os_by_stage.png",
            "km_by_stage": "km_by_stage.png",
            "rmst_by_stage": "rmst_by_stage.csv",
            "statistical_summary": "statistical_summary.txt",
        },
    },
}

class FileSettings(TypedDict):
    default_data_path: Path
    default_output_dir: Path
    figure_dpi: int


class StageOrdering(TypedDict):
    ajcc_stages: list[str]
    t_stages: list[str]
    n_stages: list[str]
    m_stages: list[str]


class VisualizationFilenames(TypedDict):
    stage_distribution: str
    tn_heatmap: str
    survival_by_stage: str
    km_by_stage: str
    rmst_by_stage: str
    statistical_summary: str


class Visualization(TypedDict):
    km_min_group_size: int
    filenames: VisualizationFilenames


class Config(TypedDict):
    file_settings: FileSettings
    column_mapping: dict[str, str]
    reverse_column_mapping: dict[str, str]
    required_columns: tuple[str, ...]
    stage_ordering: StageOrdering
    visualization: Visualization


def deep_merge_dict(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge_dict(merged[key], value)
        else:
            merged[key] = value
    return merged


def resolve_path(path_value: str | Path, relative_to_base: bool = True) -> Path:
    path = Path(path_value)
    if path.is_absolute():
        return path
    return (BASE_DIR / path) if relative_to_base else path.resolve()


def load_config(config_path: Path) -> Config:
    config_path = resolve_path(config_path, relative_to_base=False)
    base_config = copy.deepcopy(DEFAULT_CONFIG)
    if not config_path.exists():
        print(
            f"Could not locate configuration file at {config_path}. Please supply --config.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    try:
        with config_path.open("r", encoding="utf-8") as file:
            try:
                user_config = yaml.safe_load(file) or {}
            except yaml.YAMLError as exc:
                print(
                    f"Failed to parse configuration file at {config_path}: {exc}",
                    file=sys.stderr,
                )
                raise SystemExit(1)
    except OSError as exc:
        print(
            f"Could not read configuration file at {config_path}: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(1)

    if not isinstance(user_config, dict):
        print(
            f"Configuration file at {config_path} must contain a mapping of settings.",
            file=sys.stderr,
        )
        raise SystemExit(1)

    merged = deep_merge_dict(base_config, user_config)

    file_settings = merged.get("file_settings") or {}
    if not isinstance(file_settings, dict):
        print(
            "Configuration error: file_settings must be a mapping of options.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    base_file_settings = cast(dict[str, Any], base_config["file_settings"])
    default_data_path = file_settings.get(
        "default_data_path", base_file_settings["default_data_path"]
    )
    merged["file_settings"]["default_data_path"] = resolve_path(
        str(default_data_path)
    )
    default_output_dir = file_settings.get(
        "default_output_dir", base_file_settings["default_output_dir"]
    )
    merged["file_settings"]["default_output_dir"] = resolve_path(
        str(default_output_dir)
    )
    figure_dpi = file_settings.get("figure_dpi", base_file_settings["figure_dpi"])
    if figure_dpi is not None:
        merged["file_settings"]["figure_dpi"] = int(figure_dpi)

    column_mapping = merged.get("column_mapping") or {}
    if not isinstance(column_mapping, dict):
        print(
            "Configuration error: column_mapping must be a mapping of source to target columns.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    merged["column_mapping"] = dict(column_mapping)
    merged["reverse_column_mapping"] = {
        value: key for key, value in merged["column_mapping"].items()
    }

    required_columns = merged.get("required_columns") or base_config["required_columns"]
    if not isinstance(required_columns, (list, tuple)):
        print(
            "Configuration error: required_columns must be a list of normalized field names.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    merged["required_columns"] = tuple(required_columns)

    stage_ordering = merged.get("stage_ordering") or {}
    if not isinstance(stage_ordering, dict):
        print(
            "Configuration error: stage_ordering must be a mapping with ajcc_stages, t_stages, n_stages, and m_stages.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    base_stage_ordering = cast(dict[str, Any], base_config["stage_ordering"])
    ajcc_stages_value = stage_ordering.get(
        "ajcc_stages", base_stage_ordering["ajcc_stages"]
    )
    t_stages_value = stage_ordering.get("t_stages", base_stage_ordering["t_stages"])
    n_stages_value = stage_ordering.get("n_stages", base_stage_ordering["n_stages"])
    m_stages_value = stage_ordering.get("m_stages", base_stage_ordering["m_stages"])
    merged["stage_ordering"] = {
        "ajcc_stages": list(ajcc_stages_value) if ajcc_stages_value else [],
        "t_stages": list(t_stages_value) if t_stages_value else [],
        "n_stages": list(n_stages_value) if n_stages_value else [],
        "m_stages": list(m_stages_value) if m_stages_value else [],
    }

    visualization = merged.get("visualization") or {}
    if not isinstance(visualization, dict):
        print(
            "Configuration error: visualization must be a mapping of plotting options.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    base_visualization = cast(dict[str, Any], base_config["visualization"])
    km_min_group_size_value = visualization.get(
        "km_min_group_size", base_visualization["km_min_group_size"]
    )
    filenames_value = visualization.get("filenames", base_visualization["filenames"])
    merged["visualization"] = {
        "km_min_group_size": int(km_min_group_size_value) if km_min_group_size_value is not None else 0,
        "filenames": dict(filenames_value) if isinstance(filenames_value, dict) else {},
    }
    return cast(Config, merged)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate AJCC staging visualisations from TCGA STAD cohort."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=CONFIG_PATH_DEFAULT,
        help="Path to configuration YAML (default: %(default)s).",
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=None,
        help="Path to the TCGA clinical TSV (default: configured in config.yaml).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to store generated figures (default: configured in config.yaml).",
    )
    parser.add_argument(
        "--km-min-group",
        type=int,
        default=None,
        help="Minimum patients per stage required to draw a Kaplan-Meier curve (default: configured in config.yaml).",
    )
    parser.add_argument(
        "--rmst-months",
        type=float,
        default=None,
        help="Optional time horizon (months) to compute restricted mean survival time by stage (default: disabled).",
    )
    return parser.parse_args()


def load_cohort(data_path: Path) -> pd.DataFrame:
    data_path = resolve_path(data_path, relative_to_base=False)
    if not data_path.exists():
        print(
            f"Could not locate cohort file at {data_path}. Please supply --data.",
            file=sys.stderr,
        )
        raise SystemExit(1)
    try:
        return pd.read_csv(
            data_path,
            sep="\t",
            na_values=["NA", "N/A", "Not Available", ""],
            dtype=str,
        )
    except pd.errors.ParserError as exc:
        print(
            f"Could not parse cohort file at {data_path}: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    except OSError as exc:
        print(
            f"Could not read cohort file at {data_path}: {exc}",
            file=sys.stderr,
        )
        raise SystemExit(1)


def clean_stage_label(value: str | float | None) -> Optional[str]:
    if pd.isna(value):
        return None
    text = str(value).strip().upper()
    if not text:
        return None
    text = text.replace("STAGE", "", 1).strip()
    if not text or text in {"NA", "N/A"}:
        return None
    return f"Stage {text}"


def normalize_tnm_code(value: str | float | None) -> Optional[str]:
    if pd.isna(value):
        return None
    text = str(value).strip().upper()
    if not text or text in {"NA", "N/A"}:
        return None
    return text


def extract_status_text(value: str | float | None) -> Optional[str]:
    if pd.isna(value):
        return None
    text = str(value)
    if ":" in text:
        text = text.split(":", 1)[1]
    text = text.strip()
    return text.title() if text else None


def is_event(value: str | float) -> bool:
    if pd.isna(value):
        return False
    text = str(value).upper()
    return "DECEASED" in text or "DEAD" in text


def ordered_counts(series: pd.Series, order: Iterable[str]) -> pd.Series:
    counts = series.dropna().value_counts()
    if counts.empty:
        return counts
    ordered = [item for item in order if item in counts.index]
    remaining = [item for item in counts.index if item not in order]
    return counts.reindex(ordered + remaining)


def determine_category_order(values: pd.Series, base_order: Iterable[str]) -> list[str]:
    """
    Preserve the preferred AJCC ordering but append any additional categories
    observed in the data (e.g., N3A/N3B) so casting to Categorical does not drop them.
    """
    base = list(base_order)
    observed = values.dropna().unique().tolist()
    if not observed:
        return base
    ordered = [item for item in base if item in observed]
    extras = [item for item in observed if item not in base]
    return ordered + extras


def validate_required_columns(
    df: pd.DataFrame,
    required: Iterable[str],
    reverse_column_map: dict[str, str],
) -> None:
    missing = [column for column in required if column not in df.columns]
    if not missing:
        return
    source_labels = ", ".join(
        reverse_column_map.get(column, column) for column in missing
    )
    normalized = ", ".join(missing)
    raise KeyError(
        "Cohort TSV is missing required fields after harmonisation. "
        f"Normalized columns not found: {normalized}. "
        f"Expected source headers: {source_labels}. "
        "Please verify the TSV matches the schema described in README.md."
    )


def preprocess(df: pd.DataFrame, config: Config) -> pd.DataFrame:
    column_map = config["column_mapping"]
    reverse_column_map = config["reverse_column_mapping"]
    required_columns = config["required_columns"]

    df = df.rename(columns=column_map)
    validate_required_columns(df, required_columns, reverse_column_map)

    df["ajcc_stage"] = df["ajcc_stage"].apply(clean_stage_label)
    for key in ("t_stage", "n_stage", "m_stage"):
        if key in df.columns:
            df[key] = df[key].apply(normalize_tnm_code)
        else:
            df[key] = np.nan

    os_months = df.get("overall_survival_months")
    if os_months is not None:
        df["overall_survival_months"] = pd.to_numeric(os_months, errors="coerce")
    else:
        df["overall_survival_months"] = pd.Series(
            np.nan, index=df.index, dtype="float64"
        )
    raw_status = df.get("overall_survival_status")
    if raw_status is not None:
        df["survival_event"] = raw_status.apply(is_event)
        df["overall_survival_status"] = raw_status.apply(extract_status_text)
    else:
        df["survival_event"] = False
        df["overall_survival_status"] = pd.Series(
            np.nan, index=df.index, dtype="object"
        )

    try:
        df = PreprocessedCohortSchema.validate(df)
    except pa.errors.SchemaError as exc:
        print(f"Data validation failed: {exc}", file=sys.stderr)
        raise ValueError(f"Cohort data failed schema validation: {exc}") from exc

    return df


def plot_stage_panels(
    t_counts: pd.Series,
    n_counts: pd.Series,
    m_counts: pd.Series,
    stage_counts: pd.Series,
    output_path: Path,
    figure_dpi: int,
) -> bool:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    axes[0, 0].bar(t_counts.index, t_counts.values, color="steelblue", alpha=0.8)
    axes[0, 0].set_title("T Stage Distribution", fontsize=12, fontweight="bold")
    axes[0, 0].set_ylabel("Patients", fontsize=10)
    axes[0, 0].grid(axis="y", alpha=0.3)

    axes[0, 1].bar(n_counts.index, n_counts.values, color="coral", alpha=0.8)
    axes[0, 1].set_title("N Stage Distribution", fontsize=12, fontweight="bold")
    axes[0, 1].set_ylabel("Patients", fontsize=10)
    axes[0, 1].grid(axis="y", alpha=0.3)

    axes[1, 0].bar(m_counts.index, m_counts.values, color="seagreen", alpha=0.8)
    axes[1, 0].set_title("M Stage Distribution", fontsize=12, fontweight="bold")
    axes[1, 0].set_ylabel("Patients", fontsize=10)
    axes[1, 0].grid(axis="y", alpha=0.3)

    axes[1, 1].bar(stage_counts.index, stage_counts.values, color="purple", alpha=0.8)
    axes[1, 1].set_title("AJCC Overall Stage", fontsize=12, fontweight="bold")
    axes[1, 1].set_ylabel("Patients", fontsize=10)
    axes[1, 1].tick_params(axis="x", rotation=45)
    axes[1, 1].grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_path, dpi=figure_dpi, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_tn_heatmap(
    df: pd.DataFrame,
    output_path: Path,
    t_stage_order: Iterable[str],
    n_stage_order: Iterable[str],
    figure_dpi: int,
) -> bool:
    tn_subset = df.dropna(subset=["t_stage", "n_stage"]).copy()
    if tn_subset.empty:
        return False

    t_categories = determine_category_order(tn_subset["t_stage"], t_stage_order)
    n_categories = determine_category_order(tn_subset["n_stage"], n_stage_order)
    tn_subset["t_stage"] = pd.Categorical(
        tn_subset["t_stage"], categories=t_categories, ordered=True
    )
    tn_subset["n_stage"] = pd.Categorical(
        tn_subset["n_stage"], categories=n_categories, ordered=True
    )

    heatmap_data = tn_subset.pivot_table(
        index="n_stage",
        columns="t_stage",
        aggfunc="size",
        fill_value=0,
        observed=False,
    )
    heatmap_data = heatmap_data.reindex(
        index=[stage for stage in n_categories if stage in heatmap_data.index],
        columns=[stage for stage in t_categories if stage in heatmap_data.columns],
        fill_value=0,
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=".0f",
        cmap="magma",
        cbar_kws={"label": "Patient count"},
        ax=ax,
    )
    ax.set_title("T vs N Stage Burden (AJCC)", fontsize=13, fontweight="bold")
    ax.set_xlabel("T Stage", fontsize=11)
    ax.set_ylabel("N Stage", fontsize=11)
    plt.tight_layout()
    fig.savefig(output_path, dpi=figure_dpi, bbox_inches="tight")
    plt.close(fig)
    return True


def build_survival_summary(
    df: pd.DataFrame,
    stage_order: Iterable[str],
) -> pd.DataFrame:
    survival = df.dropna(subset=["overall_survival_months", "ajcc_stage"]).copy()
    if survival.empty:
        return pd.DataFrame()

    kmf = KaplanMeierFitter()
    results: list[dict[str, float | int | str]] = []

    for stage in survival["ajcc_stage"].unique():
        stage_df = survival[survival["ajcc_stage"] == stage]
        kmf.fit(
            durations=stage_df["overall_survival_months"],
            event_observed=stage_df["survival_event"],
        )

        median_os = kmf.median_survival_time_
        ci = median_survival_times(kmf.confidence_interval_)
        median_os_value = (
            float(median_os) if median_os is not None and np.isfinite(median_os) else np.nan
        )
        ci_lower = ci.iloc[0, 0] if not ci.empty else np.nan
        ci_upper = ci.iloc[0, 1] if not ci.empty else np.nan
        ci_lower_value = float(ci_lower) if np.isfinite(ci_lower) else np.nan
        ci_upper_value = float(ci_upper) if np.isfinite(ci_upper) else np.nan

        results.append(
            {
                "ajcc_stage": stage,
                "median_os": median_os_value,
                "median_os_ci_lower": ci_lower_value,
                "median_os_ci_upper": ci_upper_value,
                "patient_count": len(stage_df),
                "event_rate": stage_df["survival_event"].mean(),
            }
        )

    summary = pd.DataFrame(results).set_index("ajcc_stage")
    ordered_index = [
        stage for stage in stage_order if stage in summary.index
    ] + [stage for stage in summary.index if stage not in stage_order]
    summary = summary.reindex(ordered_index)
    summary["event_rate_pct"] = summary["event_rate"] * 100
    return summary


def compute_logrank_statistics(
    df: pd.DataFrame, min_group_size: int
) -> Optional[dict]:
    """
    Compute omnibus and pairwise log-rank statistics across AJCC stages.
    Filters out underpowered stages to avoid unstable estimates. Pairwise p-values
    are adjusted using Benjamini-Hochberg FDR correction.
    """
    km_df = df.dropna(subset=["overall_survival_months", "ajcc_stage"]).copy()
    if km_df.empty:
        return None

    stage_counts = km_df["ajcc_stage"].value_counts()
    valid_stages = stage_counts[stage_counts >= min_group_size].index
    km_df = km_df[km_df["ajcc_stage"].isin(valid_stages)]

    if len(valid_stages) < 2:
        return None

    omnibus = multivariate_logrank_test(
        km_df["overall_survival_months"],
        km_df["ajcc_stage"],
        km_df["survival_event"],
    )

    pairwise_summary = pairwise_logrank_test(
        km_df["overall_survival_months"],
        km_df["ajcc_stage"],
        event_observed=km_df["survival_event"],
    ).summary

    pairwise_results: list[dict[str, Any]] = []
    for (stage_a, stage_b), row in pairwise_summary.iterrows():
        pairwise_results.append(
            {
                "stages": (stage_a, stage_b),
                "test_statistic": float(row["test_statistic"]),
                "p_value": float(row["p"]),
            }
        )
    raw_pvalues = [result["p_value"] for result in pairwise_results]

    if raw_pvalues:
        _, adjusted_pvalues, _, _ = multipletests(
            raw_pvalues, alpha=0.05, method="fdr_bh"
        )
        for result, adjusted_p in zip(pairwise_results, adjusted_pvalues):
            result["p_value_adjusted"] = float(adjusted_p)

    pairwise_results.sort(
        key=lambda item: cast(float, item.get("p_value_adjusted", item["p_value"]))
    )

    return {
        "test_statistic": omnibus.test_statistic,
        "p_value": omnibus.p_value,
        "degrees_freedom": len(valid_stages) - 1,
        "stages_compared": list(valid_stages),
        "pairwise_results": pairwise_results,
        "correction_method": "Benjamini-Hochberg FDR",
    }


def bootstrap_rmst_ci(
    durations: pd.Series,
    events: pd.Series,
    horizon: float,
    n_bootstrap: int = 1000,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """Compute RMST with bootstrap confidence intervals."""
    if horizon <= 0:
        raise ValueError("horizon must be positive for RMST computation.")
    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be a positive integer.")

    n = len(durations)
    if n == 0:
        return (np.nan, np.nan, np.nan)

    kmf = KaplanMeierFitter()
    bootstrap_rmsts: list[float] = []

    for _ in range(n_bootstrap):
        idx = np.random.choice(n, size=n, replace=True)
        kmf.fit(durations.iloc[idx], event_observed=events.iloc[idx])
        bootstrap_rmsts.append(restricted_mean_survival_time(kmf, t=horizon))

    kmf.fit(durations, event_observed=events)
    point_est = restricted_mean_survival_time(kmf, t=horizon)

    ci_lower = np.percentile(bootstrap_rmsts, 100 * alpha / 2)
    ci_upper = np.percentile(bootstrap_rmsts, 100 * (1 - alpha / 2))

    return point_est, ci_lower, ci_upper


def compute_rmst_by_stage(
    df: pd.DataFrame,
    horizon_months: float,
    stage_order: Iterable[str],
    min_group_size: int,
    n_bootstrap: int = 1000,
) -> pd.DataFrame:
    """
    Compute restricted mean survival time (RMST) for each AJCC stage that meets the
    minimum group size threshold. RMST integrates the Kaplan-Meier survival curve
    up to the supplied time horizon to provide a censoring-robust summary metric.
    Bootstrap resampling provides 95% confidence intervals to convey uncertainty.
    """
    if horizon_months <= 0:
        raise ValueError("rmst_months must be a positive value.")
    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be a positive integer.")

    km_df = df.dropna(subset=["overall_survival_months", "ajcc_stage"]).copy()
    if km_df.empty:
        return pd.DataFrame()

    ordered_stages = [
        stage for stage in stage_order if stage in km_df["ajcc_stage"].unique()
    ] + [
        stage
        for stage in km_df["ajcc_stage"].unique()
        if stage not in stage_order
    ]

    results: list[dict[str, float | int | str]] = []

    for stage in ordered_stages:
        stage_df = km_df[km_df["ajcc_stage"] == stage]
        if len(stage_df) < min_group_size:
            continue

        rmst_value, ci_lower, ci_upper = bootstrap_rmst_ci(
            stage_df["overall_survival_months"],
            stage_df["survival_event"],
            horizon=horizon_months,
            n_bootstrap=n_bootstrap,
        )
        results.append(
            {
                "ajcc_stage": stage,
                "rmst_months": float(rmst_value),
                "rmst_ci_lower": float(ci_lower),
                "rmst_ci_upper": float(ci_upper),
                "patient_count": int(len(stage_df)),
                "events": int(stage_df["survival_event"].sum()),
            }
        )

    return pd.DataFrame(results)


def plot_survival_by_stage(
    summary: pd.DataFrame,
    output_path: Path,
    min_group_size: int,
    figure_dpi: int,
) -> bool:
    if summary.empty:
        return False

    excluded_stages = summary[summary["patient_count"] < min_group_size].index.tolist()
    summary = summary[summary["patient_count"] >= min_group_size].copy()
    if summary.empty:
        if excluded_stages:
            print(
                f"Skipped {output_path.name}: no stages met min_group_size="
                f"{min_group_size} (excluded: {', '.join(excluded_stages)})"
            )
        return False

    if excluded_stages:
        print(
            f"Excluded stages from survival by stage plot "
            f"(min_group_size={min_group_size}): {', '.join(excluded_stages)}"
        )

    fig, ax = plt.subplots(figsize=(12, 7))

    new_labels = [f"{stage}\n(n={int(count)})" for stage, count in zip(summary.index, summary["patient_count"])]

    bars = ax.bar(new_labels, summary["median_os"], color="slateblue", alpha=0.85)

    ax.set_ylabel("Median OS (months)", fontsize=11)
    ax.set_title("Overall Survival by AJCC Stage (TCGA STAD)", fontsize=13, fontweight="bold")
    ax.set_xlabel("Stage", fontsize=11)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    ax.grid(axis="y", alpha=0.3)

    for bar, rate in zip(bars, summary["event_rate_pct"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{rate:.0f}% events",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    plt.tight_layout()
    fig.savefig(output_path, dpi=figure_dpi, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_kaplan_meier(
    df: pd.DataFrame,
    output_path: Path,
    stage_order: Iterable[str],
    min_group_size: int,
    figure_dpi: int,
) -> bool:
    km_df = df.dropna(subset=["overall_survival_months", "ajcc_stage"]).copy()
    if km_df.empty:
        return False

    ordered_stages = [
        stage for stage in stage_order if stage in km_df["ajcc_stage"].unique()
    ] + [
        stage for stage in km_df["ajcc_stage"].unique() if stage not in stage_order
    ]

    kmf = KaplanMeierFitter()
    plotted = False

    fig, ax = plt.subplots(figsize=(12, 8))
    for stage in ordered_stages:
        stage_df = km_df[km_df["ajcc_stage"] == stage]
        if len(stage_df) < min_group_size:
            continue
        kmf.fit(
            durations=stage_df["overall_survival_months"],
            event_observed=stage_df["survival_event"],
            label=f"{stage} (n={len(stage_df)})",
        )
        kmf.plot_survival_function(ax=ax, ci_show=True)
        plotted = True

    if not plotted:
        plt.close(fig)
        return False

    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Overall survival probability", fontsize=11)
    ax.set_xlabel("Months since diagnosis", fontsize=11)
    ax.set_title(
        "Kaplan–Meier Survival by AJCC Stage (TCGA STAD)",
        fontsize=13,
        fontweight="bold",
    )
    ax.grid(alpha=0.3)
    ax.legend(title="Stage", fontsize=9)
    plt.tight_layout()
    fig.savefig(output_path, dpi=figure_dpi, bbox_inches="tight")
    plt.close(fig)
    return True


def print_summary(
    stage_counts: pd.Series,
    survival_summary: pd.DataFrame,
    logrank_results: Optional[dict] = None,
    rmst_summary: Optional[pd.DataFrame] = None,
    rmst_horizon: Optional[float] = None,
) -> str:
    lines: list[str] = [
        "Gastric Cancer Staging Visualization",
        "=" * 60,
    ]

    if stage_counts.empty:
        lines.append("No AJCC staging data available in the supplied cohort.")
        text = "\n".join(lines)
        print(text)
        return text

    total = int(stage_counts.sum())
    lines.append(f"Total patients with AJCC staging: {total}")
    for stage, count in stage_counts.items():
        percent = (count / total * 100) if total else 0
        lines.append(f"{stage:>10}: {int(count)} patients ({percent:.1f}%)")

    if survival_summary.empty:
        text = "\n".join(lines)
        print(text)
        return text

    lines.append("")
    lines.append("Overall Survival Highlights:")
    lines.append("-" * 60)
    for stage, row in survival_summary.iterrows():
        ci_lower = row.get("median_os_ci_lower", np.nan)
        ci_upper = row.get("median_os_ci_upper", np.nan)
        median_os = row.get("median_os", np.nan)

        median_str = f"{median_os:.1f}" if pd.notna(median_os) else "NR"
        if pd.notna(ci_lower) and pd.notna(ci_upper):
            ci_str = f" (95% CI: {ci_lower:.1f}–{ci_upper:.1f})"
        else:
            ci_str = " (95% CI: NR)"

        lines.append(
            f"{stage:>10} | Median OS: {median_str} mo{ci_str} | "
            f"Events: {row['event_rate_pct']:.0f}% | n={int(row['patient_count'])}"
        )

    if logrank_results:
        lines.append("")
        lines.append("Log-Rank Test (Omnibus):")
        lines.append(f"  Chi-square: {logrank_results['test_statistic']:.2f}")
        lines.append(f"  p-value: {logrank_results['p_value']:.4f}")
        lines.append(f"  df: {logrank_results['degrees_freedom']}")

        pairwise_results = logrank_results.get("pairwise_results") or []
        if pairwise_results:
            lines.append("  Pairwise comparisons (p-values):")
            lines.append(
                f"  Note: Adjusted using {logrank_results.get('correction_method', 'N/A')}"
            )
            max_pairs = 6
            for pair_result in pairwise_results[:max_pairs]:
                stage_a, stage_b = pair_result["stages"]
                raw_p = pair_result["p_value"]
                adj_p = pair_result.get("p_value_adjusted", raw_p)
                lines.append(
                    f"    {stage_a} vs {stage_b}: p={raw_p:.4f} "
                    f"(adj. p={adj_p:.4f}, chi^2={pair_result['test_statistic']:.2f})"
                )
            remaining = len(pairwise_results) - max_pairs
            if remaining > 0:
                lines.append(f"    ... {remaining} additional comparisons not shown")

    if rmst_summary is not None and not rmst_summary.empty and rmst_horizon:
        lines.append("")
        lines.append(
            f"Restricted Mean Survival Time (t={rmst_horizon:.1f} months):"
        )
        lines.append("-" * 60)
        for _, row in rmst_summary.iterrows():
            rmst_ci_lower = row.get("rmst_ci_lower", np.nan)
            rmst_ci_upper = row.get("rmst_ci_upper", np.nan)
            rmst_value = row.get("rmst_months", np.nan)
            rmst_str = f"{rmst_value:.1f}" if pd.notna(rmst_value) else "NR"
            if pd.notna(rmst_ci_lower) and pd.notna(rmst_ci_upper):
                ci_str = f" (95% CI: {rmst_ci_lower:.1f}–{rmst_ci_upper:.1f})"
            else:
                ci_str = " (95% CI: NR)"
            lines.append(
                f"{row['ajcc_stage']:>10} | RMST: {rmst_str} mo{ci_str} | "
                f"Events: {int(row['events'])} | n={int(row['patient_count'])}"
            )

    text = "\n".join(lines)
    print(text)
    return text


def main() -> None:
    args = parse_args()
    config: Config = load_config(args.config)
    file_settings = config["file_settings"]
    stage_ordering = config["stage_ordering"]
    visualization = config["visualization"]

    data_path = (
        resolve_path(args.data, relative_to_base=False)
        if args.data
        else file_settings["default_data_path"]
    )
    figures_dir = (
        resolve_path(args.output_dir, relative_to_base=False)
        if args.output_dir
        else file_settings["default_output_dir"]
    )
    km_min_group = (
        args.km_min_group
        if args.km_min_group is not None
        else visualization["km_min_group_size"]
    )
    figure_dpi = file_settings["figure_dpi"]
    figure_filenames = visualization["filenames"]

    df_raw = load_cohort(data_path)
    try:
        df = preprocess(df_raw, config)
    except (ValueError, KeyError) as exc:
        print(f"Data validation failed: {exc}", file=sys.stderr)
        raise SystemExit(1)

    figures_dir.mkdir(parents=True, exist_ok=True)

    stage_counts = ordered_counts(df["ajcc_stage"], stage_ordering["ajcc_stages"])
    t_counts = ordered_counts(df["t_stage"], stage_ordering["t_stages"])
    n_counts = ordered_counts(df["n_stage"], stage_ordering["n_stages"])
    m_counts = ordered_counts(df["m_stage"], stage_ordering["m_stages"])

    generated_files: list[Path] = []
    skipped_outputs: list[str] = []

    stage_distribution_path = figures_dir / figure_filenames["stage_distribution"]
    plot_stage_panels(
        t_counts=t_counts,
        n_counts=n_counts,
        m_counts=m_counts,
        stage_counts=stage_counts,
        output_path=stage_distribution_path,
        figure_dpi=figure_dpi,
    )
    generated_files.append(stage_distribution_path)

    heatmap_path = figures_dir / figure_filenames["tn_heatmap"]
    if plot_tn_heatmap(
        df,
        heatmap_path,
        t_stage_order=stage_ordering["t_stages"],
        n_stage_order=stage_ordering["n_stages"],
        figure_dpi=figure_dpi,
    ):
        generated_files.append(heatmap_path)
    else:
        skipped_outputs.append(f"{heatmap_path} (insufficient paired T/N annotations)")

    survival_summary = build_survival_summary(df, stage_ordering["ajcc_stages"])
    survival_path = figures_dir / figure_filenames["survival_by_stage"]
    if plot_survival_by_stage(
        survival_summary,
        survival_path,
        min_group_size=km_min_group,
        figure_dpi=figure_dpi,
    ):
        generated_files.append(survival_path)
    else:
        skipped_outputs.append(
            f"{survival_path} (insufficient AJCC stage + survival duration data or "
            f"no stages met min_group_size={km_min_group})"
        )

    km_path = figures_dir / figure_filenames["km_by_stage"]
    if plot_kaplan_meier(
        df,
        km_path,
        stage_order=stage_ordering["ajcc_stages"],
        min_group_size=km_min_group,
        figure_dpi=figure_dpi,
    ):
        generated_files.append(km_path)
    else:
        skipped_outputs.append(
            f"{km_path} (no stages met km_min_group={km_min_group})"
        )

    logrank_results = compute_logrank_statistics(df, min_group_size=km_min_group)
    rmst_summary: Optional[pd.DataFrame] = None
    if args.rmst_months is not None:
        if args.rmst_months <= 0:
            print("rmst-months must be a positive value.", file=sys.stderr)
            raise SystemExit(1)
        rmst_summary = compute_rmst_by_stage(
            df,
            horizon_months=args.rmst_months,
            stage_order=stage_ordering["ajcc_stages"],
            min_group_size=km_min_group,
        )
        if not rmst_summary.empty:
            rmst_output_path = figures_dir / figure_filenames["rmst_by_stage"]
            rmst_summary.to_csv(rmst_output_path, index=False)
            generated_files.append(rmst_output_path)
        else:
            skipped_outputs.append(
                f"{figures_dir / figure_filenames['rmst_by_stage']} (insufficient AJCC stage + survival duration data "
                f"or no stages met min_group_size={km_min_group} for RMST)"
            )

    summary_text = print_summary(
        stage_counts, survival_summary, logrank_results, rmst_summary, args.rmst_months
    )
    stats_output_path = figures_dir / figure_filenames["statistical_summary"]
    stats_output_path.write_text(summary_text + "\n", encoding="utf-8")
    generated_files.append(stats_output_path)
    print("\nGenerated files:")
    if generated_files:
        for path in generated_files:
            print(f"  - {path}")
    else:
        print("  - None (no figures generated)")
    if skipped_outputs:
        print("\nSkipped outputs:")
        for message in skipped_outputs:
            print(f"  - {message}")


if __name__ == "__main__":
    main()
