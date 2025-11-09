"""
Gastric cancer staging visualizer backed by TCGA STAD clinical data.

The script ingests the publicly available TCGA PanCanAtlas 2018 gastric cancer
clinical TSV, harmonises AJCC TNM annotations, and outputs figures that can
support KLASS-standardised workflow discussions.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_style("whitegrid")
plt.rcParams["figure.figsize"] = (14, 8)

BASE_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_PATH = BASE_DIR / "data" / "tcga_2018_clinical_data.tsv"

COLUMN_MAP = {
    "Patient ID": "patient_id",
    "Sample ID": "sample_id",
    "Neoplasm Disease Stage American Joint Committee on Cancer Code": "ajcc_stage",
    "American Joint Committee on Cancer Tumor Stage Code": "t_stage",
    "Neoplasm Disease Lymph Node Stage American Joint Committee on Cancer Code": "n_stage",
    "American Joint Committee on Cancer Metastasis Stage Code": "m_stage",
    "Overall Survival (Months)": "overall_survival_months",
    "Overall Survival Status": "overall_survival_status",
}

REVERSE_COLUMN_MAP = {value: key for key, value in COLUMN_MAP.items()}
REQUIRED_COLUMNS = (
    "patient_id",
    "ajcc_stage",
    "t_stage",
    "n_stage",
    "m_stage",
    "overall_survival_months",
    "overall_survival_status",
)

STAGE_ORDER = [
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
]

T_STAGE_ORDER = [
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
]

N_STAGE_ORDER = ["N0", "N1", "N2", "N3", "N3A", "N3B", "NX"]
M_STAGE_ORDER = ["M0", "M1", "MX"]

FIG_STAGE_DISTRIBUTION = "tnm_staging_distribution.png"
FIG_TN_HEATMAP = "tn_heatmap.png"
FIG_SURVIVAL = "os_by_stage.png"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate AJCC staging visualisations from TCGA STAD cohort."
    )
    parser.add_argument(
        "--data",
        type=Path,
        default=DEFAULT_DATA_PATH,
        help="Path to the TCGA clinical TSV (default: %(default)s).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=BASE_DIR,
        help="Directory to store generated figures (default: project root).",
    )
    return parser.parse_args()


def load_cohort(data_path: Path) -> pd.DataFrame:
    if not data_path.exists():
        raise FileNotFoundError(
            f"Could not locate cohort file at {data_path}. Please supply --data."
        )
    return pd.read_csv(
        data_path,
        sep="\t",
        na_values=["NA", "N/A", "Not Available", ""],
        dtype=str,
    )


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


def validate_required_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    missing = [column for column in required if column not in df.columns]
    if not missing:
        return
    source_labels = ", ".join(REVERSE_COLUMN_MAP.get(column, column) for column in missing)
    normalized = ", ".join(missing)
    raise ValueError(
        "Cohort TSV is missing required fields after harmonisation. "
        f"Normalized columns not found: {normalized}. "
        f"Expected source headers: {source_labels}. "
        "Please verify the TSV matches the schema described in README.md."
    )


def preprocess(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns=COLUMN_MAP)
    validate_required_columns(df, REQUIRED_COLUMNS)

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
    return df


def plot_stage_panels(
    t_counts: pd.Series,
    n_counts: pd.Series,
    m_counts: pd.Series,
    stage_counts: pd.Series,
    output_path: Path,
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
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_tn_heatmap(df: pd.DataFrame, output_path: Path) -> bool:
    tn_subset = df.dropna(subset=["t_stage", "n_stage"]).copy()
    if tn_subset.empty:
        return False

    t_categories = determine_category_order(tn_subset["t_stage"], T_STAGE_ORDER)
    n_categories = determine_category_order(tn_subset["n_stage"], N_STAGE_ORDER)
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
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def build_survival_summary(df: pd.DataFrame) -> pd.DataFrame:
    survival = df.dropna(subset=["overall_survival_months", "ajcc_stage"]).copy()
    if survival.empty:
        return pd.DataFrame()

    summary = (
        survival.groupby("ajcc_stage")
        .agg(
            median_os=("overall_survival_months", "median"),
            patient_count=("patient_id", "count"),
            event_rate=("survival_event", "mean"),
        )
        .dropna()
    )

    ordered_index = [
        stage for stage in STAGE_ORDER if stage in summary.index
    ] + [stage for stage in summary.index if stage not in STAGE_ORDER]
    summary = summary.reindex(ordered_index)
    summary["event_rate_pct"] = summary["event_rate"] * 100
    return summary


def plot_survival_by_stage(summary: pd.DataFrame, output_path: Path) -> bool:
    if summary.empty:
        return False

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(summary.index, summary["median_os"], color="slateblue", alpha=0.85)
    ax.set_ylabel("Median OS (months)", fontsize=11)
    ax.set_title("Overall Survival by AJCC Stage (TCGA STAD)", fontsize=13, fontweight="bold")
    ax.set_xlabel("Stage", fontsize=11)
    ax.tick_params(axis="x", rotation=45)
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
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def print_summary(stage_counts: pd.Series, survival_summary: pd.DataFrame) -> None:
    if stage_counts.empty:
        print("Gastric Cancer Staging Visualization")
        print("=" * 60)
        print("No AJCC staging data available in the supplied cohort.")
        return

    total = int(stage_counts.sum())
    print("Gastric Cancer Staging Visualization")
    print("=" * 60)
    print(f"Total patients with AJCC staging: {total}")
    for stage, count in stage_counts.items():
        percent = (count / total * 100) if total else 0
        print(f"{stage:>10}: {int(count)} patients ({percent:.1f}%)")

    if survival_summary.empty:
        return

    print("\nOverall Survival Highlights:")
    print("-" * 60)
    for stage, row in survival_summary.iterrows():
        print(
            f"{stage:>10} | Median OS: {row['median_os']:.1f} mo | "
            f"Events: {row['event_rate_pct']:.0f}% | n={int(row['patient_count'])}"
        )


def main() -> None:
    args = parse_args()
    df_raw = load_cohort(args.data)
    df = preprocess(df_raw)

    figures_dir = args.output_dir
    figures_dir.mkdir(parents=True, exist_ok=True)

    stage_counts = ordered_counts(df["ajcc_stage"], STAGE_ORDER)
    t_counts = ordered_counts(df["t_stage"], T_STAGE_ORDER)
    n_counts = ordered_counts(df["n_stage"], N_STAGE_ORDER)
    m_counts = ordered_counts(df["m_stage"], M_STAGE_ORDER)

    generated_files: list[Path] = []
    skipped_outputs: list[str] = []

    stage_distribution_path = figures_dir / FIG_STAGE_DISTRIBUTION
    plot_stage_panels(
        t_counts=t_counts,
        n_counts=n_counts,
        m_counts=m_counts,
        stage_counts=stage_counts,
        output_path=stage_distribution_path,
    )
    generated_files.append(stage_distribution_path)

    heatmap_path = figures_dir / FIG_TN_HEATMAP
    if plot_tn_heatmap(df, heatmap_path):
        generated_files.append(heatmap_path)
    else:
        skipped_outputs.append(f"{heatmap_path} (insufficient paired T/N annotations)")

    survival_summary = build_survival_summary(df)
    survival_path = figures_dir / FIG_SURVIVAL
    if plot_survival_by_stage(survival_summary, survival_path):
        generated_files.append(survival_path)
    else:
        skipped_outputs.append(
            f"{survival_path} (insufficient AJCC stage + survival duration data)"
        )

    print_summary(stage_counts, survival_summary)
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
