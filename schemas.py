"""Pandera schemas for TCGA clinical data validation."""
from typing import Optional

import pandera.pandas as pa
from pandera.typing import Series


VALID_AJCC_STAGES = [
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

VALID_T_STAGES = [
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

VALID_N_STAGES = ["N0", "N1", "N2", "N3", "N3A", "N3B", "NX"]
VALID_M_STAGES = ["M0", "M1", "MX"]


class PreprocessedCohortSchema(pa.DataFrameModel):
    """Schema for cohort data after preprocessing."""

    patient_id: Series[str] = pa.Field(unique=True, nullable=False)

    ajcc_stage: Optional[Series[str]] = pa.Field(
        nullable=True,
        isin=VALID_AJCC_STAGES,
    )

    t_stage: Optional[Series[str]] = pa.Field(
        nullable=True,
        isin=VALID_T_STAGES,
    )

    n_stage: Optional[Series[str]] = pa.Field(
        nullable=True,
        isin=VALID_N_STAGES,
    )

    m_stage: Optional[Series[str]] = pa.Field(
        nullable=True,
        isin=VALID_M_STAGES,
    )

    overall_survival_months: Optional[Series[float]] = pa.Field(
        nullable=True,
        ge=0.0,  # Non-negative survival time
        le=360.0,  # Reasonable upper bound (30 years)
    )

    survival_event: Series[bool] = pa.Field(nullable=False)

    class Config:
        coerce = True  # Attempt type coercion
        strict = False  # Allow extra columns
