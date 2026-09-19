"""Explicit mortality layout, covariate coding and hip metric reshaping."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from .questionnaires import assert_unique

PAX_VARIANTS = ["primary", "valid_days_3", "valid_days_5", "calibration_1_or_2", "nonwear_90"]

def parse_mortality(path: Path) -> pd.DataFrame:
    """Parse official NHANES fixed-width positions from the official NCHS public-use read-in specification."""
    names = [
        "SEQN",
        "ELIGSTAT",
        "MORTSTAT",
        "UCOD_LEADING",
        "DIABETES_DEATH_FLAG",
        "HYPERTEN_DEATH_FLAG",
        "PERMTH_INT",
        "PERMTH_EXM",
    ]
    colspecs = [(0, 6), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21), (42, 45), (45, 48)]
    out = pd.read_fwf(
        path,
        colspecs=colspecs,
        names=names,
        dtype="Int64",
        na_values=["", "."],
        keep_default_na=True,
    )
    assert_unique(out, "SEQN", "mortality file")
    return out


def derive_covariates(demo: pd.DataFrame, smq: pd.DataFrame, bmx: pd.DataFrame,
                      mcq: pd.DataFrame, diq: pd.DataFrame) -> pd.DataFrame:
    for label, frame in [("DEMO", demo), ("SMQ", smq), ("BMX", bmx), ("MCQ", mcq), ("DIQ", diq)]:
        assert_unique(frame, "SEQN", label)
    out = demo[
        [
            "SEQN",
            "RIDSTATR",
            "RIDAGEYR",
            "RIAGENDR",
            "RIDRETH1",
            "DMDEDUC2",
            "INDFMPIR",
            "WTMEC2YR",
            "SDMVSTRA",
            "SDMVPSU",
        ]
    ].copy()
    out = out.rename(columns={"RIDAGEYR": "age", "INDFMPIR": "PIR"})
    out["sex"] = out["RIAGENDR"].map({1.0: "male", 2.0: "female"})
    out["race_ethnicity4"] = out["RIDRETH1"].map(
        {
            1.0: "Hispanic",
            2.0: "Hispanic",
            3.0: "Non-Hispanic White",
            4.0: "Non-Hispanic Black",
            5.0: "Other",
        }
    )
    out["education3"] = out["DMDEDUC2"].map(
        {
            1.0: "under high school",
            2.0: "under high school",
            3.0: "high school/GED",
            4.0: "more than high school",
            5.0: "more than high school",
        }
    )
    out["education_binary"] = out["DMDEDUC2"].map(
        {
            1.0: "high school or less",
            2.0: "high school or less",
            3.0: "high school or less",
            4.0: "more than high school",
            5.0: "more than high school",
        }
    )

    out = out.merge(smq[["SEQN", "SMQ020", "SMQ040"]], on="SEQN", how="left", validate="one_to_one")
    out["smoking3"] = pd.Series(pd.NA, index=out.index, dtype="string")
    out.loc[out["SMQ020"] == 2, "smoking3"] = "never"
    out.loc[(out["SMQ020"] == 1) & (out["SMQ040"] == 3), "smoking3"] = "former"
    out.loc[(out["SMQ020"] == 1) & out["SMQ040"].isin([1, 2]), "smoking3"] = "current"
    out["smoking_current"] = out["smoking3"].map(
        {"never": "not current", "former": "not current", "current": "current"}
    )

    out = out.merge(bmx[["SEQN", "BMXBMI"]], on="SEQN", how="left", validate="one_to_one")
    out["BMI"] = out["BMXBMI"].where(out["BMXBMI"] > 0)
    out = out.merge(diq[["SEQN", "DIQ010"]], on="SEQN", how="left", validate="one_to_one")
    out["diabetes_history"] = out["DIQ010"].map({1.0: "yes", 2.0: "no", 3.0: "no"})

    cvd_items = ["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]
    out = out.merge(mcq[["SEQN", *cvd_items, "MCQ220"]], on="SEQN", how="left", validate="one_to_one")
    any_yes = out[cvd_items].eq(1).any(axis=1)
    all_no = out[cvd_items].eq(2).all(axis=1)
    out["CVD_history"] = pd.Series(pd.NA, index=out.index, dtype="string")
    out.loc[any_yes, "CVD_history"] = "yes"
    out.loc[~any_yes & all_no, "CVD_history"] = "no"
    out["cancer_history"] = out["MCQ220"].map({1.0: "yes", 2.0: "no"})
    return out


def load_pax(path: Path) -> pd.DataFrame:
    pax = pd.read_csv(path, float_precision="round_trip")
    required = {
        "SEQN",
        "variant",
        "eligible_participant",
        "log1p_counts_per_wear_minute",
    }
    if not required.issubset(pax.columns):
        raise ValueError("PAX participant output schema changed")
    if set(pax["variant"].unique()) != set(PAX_VARIANTS):
        raise ValueError("PAX variants differ from frozen set")
    if pax.duplicated(["SEQN", "variant"]).any():
        raise ValueError("PAX participant key is not unique")
    metric_cols = [c for c in pax.columns if c not in ["SEQN", "variant"]]
    wide = pax.pivot(index="SEQN", columns="variant", values=metric_cols)
    wide.columns = [f"PAX_{metric}__{variant}" for metric, variant in wide.columns]
    wide = wide.reset_index()
    eligible = "PAX_eligible_participant__primary"
    primary_metric = "PAX_log1p_counts_per_wear_minute__primary"
    wide["PAX_primary_eligible"] = wide[eligible].astype("boolean")
    wide["PAX_log_counts_primary"] = wide[primary_metric].where(wide["PAX_primary_eligible"])
    return wide


def weighted_median(values: pd.Series, weights: pd.Series) -> float:
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    order = np.argsort(x, kind="mergesort")
    xs, ws = x[order], w[order]
    idx = int(np.searchsorted(np.cumsum(ws), 0.5 * ws.sum(), side="left"))
    return float(xs[min(idx, len(xs) - 1)])


def assign_discordance_group(s: float, d: float, s_median: float, d_median: float) -> str:
    s_high = s > s_median
    d_high = d > d_median
    if not s_high and not d_high:
        return "concordant_low"
    if s_high and d_high:
        return "concordant_high"
    if s_high and not d_high:
        return "self_report_relatively_high"
    return "self_report_relatively_low"
