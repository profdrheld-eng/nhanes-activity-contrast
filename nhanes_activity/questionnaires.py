"""Questionnaire scoring and released wrist-day summaries.

These are the study-specific rules, including plausibility limits, rather than
universal NHANES validity rules. Negative gates are zero, unknown/unable values
remain missing, and contradictory hip leisure records fail explicitly.
"""
from __future__ import annotations
import math
import numpy as np
import pandas as pd

def assert_unique(frame: pd.DataFrame, key: str, label: str) -> None:
    if frame[key].isna().any() or frame[key].duplicated().any():
        raise ValueError(f"{label}: {key} must be nonmissing and unique")


def numeric_between(series: pd.Series, low: float, high: float) -> pd.Series:
    return series.notna() & series.between(low, high, inclusive="both")


def hip_self_report(paq: pd.DataFrame, paqiaf: pd.DataFrame, adult_ids: set[int]) -> pd.DataFrame:
    """Implement the study self-report score rules for 2003–2006."""
    assert_unique(paq, "SEQN", "PAQ")
    required_activity = {"SEQN", "PADLEVEL", "PADTIMES", "PADDURAT"}
    if not required_activity.issubset(paqiaf.columns):
        raise ValueError("PAQIAF missing required activity columns")

    activity = paqiaf.loc[paqiaf["SEQN"].isin(adult_ids)].copy()
    if activity["SEQN"].isna().any():
        raise ValueError("PAQIAF contains missing SEQN")
    valid_level = activity["PADLEVEL"].isin([1, 2])
    if not valid_level.all():
        raise ValueError("Adult PAQIAF contains invalid PADLEVEL")

    summaries: dict[int, dict[str, float | int | bool]] = {}
    for seqn, rows in activity.groupby("SEQN", sort=False):
        item: dict[str, float | int | bool] = {"activity_row_count": int(len(rows))}
        for level, prefix in [(1, "moderate"), (2, "vigorous")]:
            part = rows.loc[rows["PADLEVEL"] == level]
            valid = numeric_between(part["PADTIMES"], 1, 300) & numeric_between(
                part["PADDURAT"], 10, 600
            )
            item[f"{prefix}_row_count"] = int(len(part))
            item[f"{prefix}_rows_complete"] = bool(valid.all())
            item[f"{prefix}_minutes_sum"] = (
                float((part["PADTIMES"] * part["PADDURAT"]).sum()) if valid.all() else math.nan
            )
        summaries[int(seqn)] = item
    if summaries:
        leisure = pd.DataFrame.from_dict(summaries, orient="index")
        leisure.index.name = "SEQN"
        leisure = leisure.reset_index()
    else:
        leisure = pd.DataFrame(
            columns=[
                "SEQN",
                "activity_row_count",
                "moderate_row_count",
                "moderate_rows_complete",
                "moderate_minutes_sum",
                "vigorous_row_count",
                "vigorous_rows_complete",
                "vigorous_minutes_sum",
            ]
        )

    out = paq.loc[paq["SEQN"].isin(adult_ids)].copy()
    out = out.merge(leisure, on="SEQN", how="left", validate="one_to_one")
    for col in ["activity_row_count", "moderate_row_count", "vigorous_row_count"]:
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0).astype(int)
    for col in ["moderate_rows_complete", "vigorous_rows_complete"]:
        out[col] = out[col].where(out[col].notna(), True).astype(bool)

    # Fail on released gate/activity contradictions in the analytic adult population.
    moderate_mismatch = (
        ((out["PAD320"] == 1) & (out["moderate_row_count"] == 0))
        | (out["PAD320"].isin([2, 3]) & (out["moderate_row_count"] > 0))
    )
    vigorous_mismatch = (
        ((out["PAD200"] == 1) & (out["vigorous_row_count"] == 0))
        | (out["PAD200"].isin([2, 3]) & (out["vigorous_row_count"] > 0))
    )
    if moderate_mismatch.any() or vigorous_mismatch.any():
        raise ValueError("Adult leisure gate/activity-row contradiction")

    out["SR_transport_30d"] = np.nan
    gate_yes = out["PAD020"] == 1
    transport_valid = (
        numeric_between(out["PAQ050Q"], 1, 100)
        & out["PAQ050U"].isin([1, 2, 3])
        & numeric_between(out["PAD080"], 10, 600)
    )
    multipliers = out["PAQ050U"].map({1.0: 30.0, 2.0: 30.0 / 7.0, 3.0: 1.0})
    out.loc[gate_yes & transport_valid, "SR_transport_30d"] = (
        out["PAQ050Q"] * multipliers * out["PAD080"]
    )
    out.loc[out["PAD020"] == 2, "SR_transport_30d"] = 0.0

    out["SR_home_30d"] = np.nan
    home_valid = numeric_between(out["PAD120"], 1, 120) & numeric_between(
        out["PAD160"], 10, 600
    )
    out.loc[(out["PAQ100"] == 1) & home_valid, "SR_home_30d"] = (
        out["PAD120"] * out["PAD160"]
    )
    out.loc[out["PAQ100"] == 2, "SR_home_30d"] = 0.0

    out["SR_moderate_30d"] = np.nan
    moderate_complete = (
        (out["PAD320"] == 1)
        & (out["moderate_row_count"] > 0)
        & out["moderate_rows_complete"]
    )
    out.loc[moderate_complete, "SR_moderate_30d"] = out["moderate_minutes_sum"]
    out.loc[out["PAD320"] == 2, "SR_moderate_30d"] = 0.0

    out["SR_vigorous_30d"] = np.nan
    vigorous_complete = (
        (out["PAD200"] == 1)
        & (out["vigorous_row_count"] > 0)
        & out["vigorous_rows_complete"]
    )
    out.loc[vigorous_complete, "SR_vigorous_30d"] = out["vigorous_minutes_sum"]
    out.loc[out["PAD200"] == 2, "SR_vigorous_30d"] = 0.0

    domains = [
        "SR_transport_30d",
        "SR_home_30d",
        "SR_moderate_30d",
        "SR_vigorous_30d",
    ]
    out["SR_complete"] = out[domains].notna().all(axis=1)
    out["SR_day"] = np.nan
    out.loc[out["SR_complete"], "SR_day"] = (
        out["SR_transport_30d"]
        + out["SR_home_30d"]
        + out["SR_moderate_30d"]
        + 2.0 * out["SR_vigorous_30d"]
    ) / 30.0
    if (out.loc[out["SR_complete"], "SR_day"] < 0).any():
        raise ValueError("Negative self-report score")
    return out[
        [
            "SEQN",
            *domains,
            "SR_complete",
            "SR_day",
            "activity_row_count",
            "moderate_row_count",
            "vigorous_row_count",
        ]
    ]


def wrist_self_report(paq: pd.DataFrame) -> pd.DataFrame:
    """Construct the frozen GPAQ activity-equivalent score."""
    assert_unique(paq, "SEQN", "PAQ")
    domains = [
        ("vigorous_work", "PAQ605", "PAQ610", "PAD615", 2.0),
        ("moderate_work", "PAQ620", "PAQ625", "PAD630", 1.0),
        ("transport", "PAQ635", "PAQ640", "PAD645", 1.0),
        ("vigorous_recreation", "PAQ650", "PAQ655", "PAD660", 2.0),
        ("moderate_recreation", "PAQ665", "PAQ670", "PAD675", 1.0),
    ]
    required = {"SEQN"}
    for _, gate, days, duration, _ in domains:
        required.update((gate, days, duration))
    if not required.issubset(paq.columns):
        raise ValueError(f"PAQ missing variables: {sorted(required - set(paq.columns))}")

    out = paq[["SEQN"]].copy()
    score_columns = []
    for label, gate, days, duration, intensity in domains:
        score = pd.Series(np.nan, index=paq.index, dtype=float)
        valid = numeric_between(paq[days], 1, 7) & numeric_between(
            paq[duration], 10, 1440
        )
        score.loc[(paq[gate] == 1) & valid] = (
            pd.to_numeric(paq.loc[(paq[gate] == 1) & valid, days])
            * pd.to_numeric(paq.loc[(paq[gate] == 1) & valid, duration])
            * intensity
        )
        score.loc[paq[gate] == 2] = 0.0
        column = f"SR_{label}_weekly_equivalent_minutes"
        out[column] = score
        score_columns.append(column)

    out["SR_complete"] = out[score_columns].notna().all(axis=1)
    out["SR_day"] = np.nan
    out.loc[out["SR_complete"], "SR_day"] = (
        out.loc[out["SR_complete"], score_columns].sum(axis=1) / 7.0
    )
    if (out.loc[out["SR_complete"], "SR_day"] < 0).any():
        raise ValueError("Negative self-report score")
    return out


def aggregate_wrist(paxday: pd.DataFrame) -> pd.DataFrame:
    """Aggregate released daily MIMS summaries using the frozen support rules."""
    required = {"SEQN", "PAXDAYD", "PAXWWMD", "PAXVMD", "PAXMTSD"}
    if not required.issubset(paxday.columns):
        raise ValueError(f"PAXDAY missing variables: {sorted(required - set(paxday.columns))}")
    if paxday.duplicated(["SEQN", "PAXDAYD"]).any():
        raise ValueError("PAXDAY contains duplicate participant-day records")

    d = paxday[list(required)].copy()
    for column in ("PAXWWMD", "PAXVMD", "PAXMTSD"):
        d[column] = pd.to_numeric(d[column], errors="coerce")
    d["valid_day"] = (
        d["PAXWWMD"].ge(600)
        & d["PAXVMD"].gt(0)
        & d["PAXMTSD"].ge(0)
    )

    rows = []
    for seqn, participant in d.groupby("SEQN", sort=False):
        valid = participant.loc[participant["valid_day"]]
        valid_days = len(valid)
        total_mims = float(valid["PAXMTSD"].sum()) if valid_days else math.nan
        total_valid_minutes = float(valid["PAXVMD"].sum()) if valid_days else math.nan
        mean_daily = total_mims / valid_days if valid_days else math.nan
        per_valid_minute = (
            total_mims / total_valid_minutes
            if valid_days and total_valid_minutes > 0
            else math.nan
        )
        rows.append(
            {
                "SEQN": seqn,
                "PAX_days_released": len(participant),
                "PAX_valid_days": valid_days,
                "PAX_total_mims_valid_days": total_mims,
                "PAX_total_valid_minutes": total_valid_minutes,
                "PAX_mean_daily_mims": mean_daily,
                "PAX_mims_per_valid_minute": per_valid_minute,
                "PAX_primary_eligible": valid_days >= 4,
                "PAX_one_day_eligible": valid_days >= 1,
            }
        )
    return pd.DataFrame(rows)
