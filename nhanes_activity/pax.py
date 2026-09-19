"""Hip accelerometer processing for NHANES 2003–2006.

The minute rules preserve the study's historical implementation. They operate
within seven 1,440-position monitor days. Missing positions are unknown, never
imputed zero. Calibration, nonwear and minimum valid-day sensitivities are
explicit. Participant-level structural contradictions cause quarantine.
"""
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Set, Tuple
import json
import math
import numpy as np
import pandas as pd

VARIANTS = (
    ("primary", {1}, 60, 4),
    ("valid_days_3", {1}, 60, 3),
    ("valid_days_5", {1}, 60, 5),
    ("calibration_1_or_2", {1, 2}, 60, 4),
    ("nonwear_90", {1}, 90, 4),
)

class DataInvariantError(ValueError):
    """Input or derived data violate an explicitly checked invariant."""

def classify_nonwear(
    known: np.ndarray, intensity: np.ndarray, minimum: int
) -> np.ndarray:
    """Classify non-wear in a logical seven-day frame.

    Low-count runs of one or two minutes are admitted only when immediately
    bounded by zero-count minutes. Unknown positions and day boundaries split
    candidates because classification is performed one monitor day at a time.
    """

    known = np.asarray(known)
    intensity = np.asarray(intensity, dtype=float)
    if known.shape != (10080,) or intensity.shape != (10080,):
        raise ValueError("Expected two seven-day vectors with 10080 positions")
    if known.dtype != np.dtype(bool):
        raise ValueError("Known-minute mask must be boolean")
    if isinstance(minimum, (bool, np.bool_)) or not isinstance(minimum, (int, np.integer)) or minimum <= 0:
        raise ValueError("Nonwear minimum must be a positive integer")
    result = np.zeros(known.shape, dtype=bool)
    for day_index in range(7):
        start = day_index * 1440
        stop = start + 1440
        day_known = known[start:stop]
        day_intensity = intensity[start:stop]
        zero = day_known & (day_intensity == 0)
        low = day_known & (day_intensity > 0) & (day_intensity < 100)
        allowed_low = np.zeros(1440, dtype=bool)
        padded = np.concatenate(([False], low, [False]))
        changes = np.flatnonzero(padded[1:] != padded[:-1])
        for left, right in changes.reshape(-1, 2):
            if right - left <= 2 and left > 0 and right < 1440:
                if zero[left - 1] and zero[right]:
                    allowed_low[left:right] = True
        candidate = zero | allowed_low
        padded_candidate = np.concatenate(([False], candidate, [False]))
        bounds = np.flatnonzero(
            padded_candidate[1:] != padded_candidate[:-1]
        ).reshape(-1, 2)
        for left, right in bounds:
            if right - left >= minimum:
                result[start + left : start + right] = True
    return result


def _value_counts(values: np.ndarray) -> str:
    clean = ["missing" if pd.isna(value) else str(int(value)) for value in values]
    counts: Dict[str, int] = {}
    for value in clean:
        counts[value] = counts.get(value, 0) + 1
    return json.dumps(counts, sort_keys=True, separators=(",", ":"))


def participant_qc(rows: Dict[str, np.ndarray]) -> Tuple[List[str], Dict[str, object]]:
    reasons: List[str] = []
    paxn = rows["PAXN"]
    if np.isnan(paxn).any():
        reasons.append("missing_paxn")
    finite = paxn[np.isfinite(paxn)]
    if finite.size and np.any(finite != np.floor(finite)):
        reasons.append("noninteger_paxn")
    if finite.size and np.any((finite < 1) | (finite > 10080)):
        reasons.append("out_of_range_paxn")
    valid_key = (
        finite.size == paxn.size
        and np.all(finite == np.floor(finite))
        and np.all((finite >= 1) & (finite <= 10080))
    )
    duplicate_count = 0
    nonincreasing_count = 0
    if valid_key and len(paxn) > 1:
        differences = np.diff(paxn)
        duplicate_count = int(np.sum(differences == 0))
        nonincreasing_count = int(np.sum(differences < 0))
        if duplicate_count:
            reasons.append("duplicate_seqn_paxn")
        if nonincreasing_count:
            reasons.append("nonincreasing_paxn")

    for field, reason in (
        ("PAXSTAT", "varying_paxstat"),
        ("PAXCAL", "varying_paxcal"),
    ):
        unique = set("missing" if pd.isna(value) else float(value) for value in rows[field])
        if len(unique) > 1:
            reasons.append(reason)

    clock_contradictions = 0
    if valid_key and duplicate_count == 0 and nonincreasing_count == 0 and len(paxn) > 1:
        previous_day = rows["PAXDAY"][:-1]
        previous_hour = rows["PAXHOUR"][:-1]
        previous_minute = rows["PAXMINUT"][:-1]
        current_day = rows["PAXDAY"][1:]
        current_hour = rows["PAXHOUR"][1:]
        current_minute = rows["PAXMINUT"][1:]
        delta = np.diff(paxn).astype(np.int64)
        complete = (
            np.isfinite(previous_day)
            & np.isfinite(previous_hour)
            & np.isfinite(previous_minute)
            & np.isfinite(current_day)
            & np.isfinite(current_hour)
            & np.isfinite(current_minute)
        )
        total = previous_hour * 60 + previous_minute + delta
        day_advance = np.floor_divide(total, 1440)
        minute_of_day = np.mod(total, 1440)
        expected_day = np.mod(previous_day - 1 + day_advance, 7) + 1
        expected_hour = np.floor_divide(minute_of_day, 60)
        expected_minute = np.mod(minute_of_day, 60)
        contradiction = (
            ~complete
            | (current_day != expected_day)
            | (current_hour != expected_hour)
            | (current_minute != expected_minute)
        )
        clock_contradictions = int(np.sum(contradiction))
        if clock_contradictions:
            reasons.append("clock_contradiction")

    internal_gap_positions = 0
    if valid_key and duplicate_count == 0 and nonincreasing_count == 0:
        internal_gap_positions = int(np.maximum(np.diff(paxn) - 1, 0).sum())
    qc = {
        "observed_rows": int(len(paxn)),
        "duplicate_keys": duplicate_count,
        "nonincreasing_keys": nonincreasing_count,
        "internal_gap_positions": internal_gap_positions,
        "clock_contradictions": clock_contradictions,
        "paxstat_values": _value_counts(rows["PAXSTAT"]),
        "paxcal_values": _value_counts(rows["PAXCAL"]),
        "missing_intensity": int(np.isnan(rows["PAXINTEN"]).sum()),
        "boundary_intensity_32767": int(np.sum(rows["PAXINTEN"] == 32767)),
        "negative_intensity": int(np.sum(rows["PAXINTEN"] < 0)),
        "missing_step": int(np.isnan(rows["PAXSTEP"]).sum()),
        "boundary_step_32767": int(np.sum(rows["PAXSTEP"] == 32767)),
        "negative_step": int(np.sum(rows["PAXSTEP"] < 0)),
    }
    return sorted(set(reasons)), qc


def derive_variant(
    rows: Dict[str, np.ndarray],
    variant: str,
    calibration_codes: Set[int],
    nonwear_minimum: int,
    participant_day_threshold: int,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    paxn = rows["PAXN"].astype(np.int64)
    indices = paxn - 1
    observed = np.zeros(10080, dtype=bool)
    observed[indices] = True
    status = np.full(10080, np.nan)
    calibration = np.full(10080, np.nan)
    intensity = np.full(10080, np.nan)
    steps = np.full(10080, np.nan)
    status[indices] = rows["PAXSTAT"]
    calibration[indices] = rows["PAXCAL"]
    intensity[indices] = rows["PAXINTEN"]
    steps[indices] = rows["PAXSTEP"]

    known = (
        observed
        & (status == 1)
        & np.isin(calibration, list(calibration_codes))
        & np.isfinite(intensity)
        & (intensity >= 0)
        & (intensity < 32767)
    )
    nonwear = classify_nonwear(known, intensity, nonwear_minimum)
    wear = known & ~nonwear
    observed_unknown = observed & ~known
    absent = ~observed
    internal_gap = np.zeros(10080, dtype=bool)
    if len(indices):
        internal_gap[min(indices) : max(indices) + 1] = True
        internal_gap &= absent

    days: List[Dict[str, object]] = []
    for day_number in range(1, 8):
        block = slice((day_number - 1) * 1440, day_number * 1440)
        day_wear = wear[block]
        wear_minutes = int(day_wear.sum())
        nonwear_minutes = int(nonwear[block].sum())
        observed_unknown_minutes = int(observed_unknown[block].sum())
        absent_unknown_minutes = int(absent[block].sum())
        unknown_minutes = observed_unknown_minutes + absent_unknown_minutes
        if wear_minutes + nonwear_minutes + unknown_minutes != 1440:
            raise DataInvariantError("day accounting identity failed")
        valid_day = wear_minutes >= 600
        step_usable = (
            np.isfinite(steps[block])
            & (steps[block] >= 0)
            & (steps[block] < 32767)
        )
        complete = bool(valid_day and np.all(step_usable[day_wear]))
        daily_steps = float(steps[block][day_wear].sum()) if complete else None
        intensity_sum = float(intensity[block][day_wear].sum())
        day_record = {
            "variant": variant,
            "monitor_day": day_number,
            "observed_rows": int(observed[block].sum()),
            "wear_minutes": wear_minutes,
            "nonwear_minutes": nonwear_minutes,
            "observed_unknown_minutes": observed_unknown_minutes,
            "absent_unknown_minutes": absent_unknown_minutes,
            "unknown_minutes": unknown_minutes,
            "internal_gap_positions": int(internal_gap[block].sum()),
            "valid_day": valid_day,
            "wear_intensity_sum": intensity_sum,
            "counts_per_wear_minute": (
                intensity_sum / wear_minutes if wear_minutes else None
            ),
            "step_complete": complete if valid_day else None,
            "daily_steps": daily_steps,
            "sedentary_minutes": int(
                np.sum(day_wear & np.isfinite(intensity[block]) & (intensity[block] < 100))
            ),
            "mvpa_minutes": int(
                np.sum(day_wear & np.isfinite(intensity[block]) & (intensity[block] >= 2020))
            ),
        }
        days.append(day_record)

    valid_days = [day for day in days if day["valid_day"]]
    complete_days = [day for day in valid_days if day["step_complete"]]
    incomplete_days = [day for day in valid_days if not day["step_complete"]]
    eligible = len(valid_days) >= participant_day_threshold
    total_wear = int(sum(day["wear_minutes"] for day in valid_days))
    total_intensity = float(sum(day["wear_intensity_sum"] for day in valid_days))
    counts_per_wear = total_intensity / total_wear if eligible and total_wear else None
    step_metric = (
        float(np.mean([day["daily_steps"] for day in complete_days]))
        if eligible and len(complete_days) >= participant_day_threshold
        else None
    )
    person = {
        "variant": variant,
        "participant_day_threshold": participant_day_threshold,
        "nonwear_minimum": nonwear_minimum,
        "calibration_codes": ",".join(map(str, sorted(calibration_codes))),
        "intensity_valid_days": len(valid_days),
        "step_complete_days": len(complete_days),
        "step_incomplete_days": len(incomplete_days),
        "eligible_participant": eligible,
        "total_valid_day_wear_minutes": total_wear if eligible else None,
        "counts_per_wear_minute": counts_per_wear,
        "log1p_counts_per_wear_minute": (
            math.log1p(counts_per_wear) if counts_per_wear is not None else None
        ),
        "mean_daily_steps": step_metric,
        "step_metric_denominator_days": (
            len(complete_days) if step_metric is not None else None
        ),
        "mean_sedentary_minutes_per_valid_day": (
            float(np.mean([day["sedentary_minutes"] for day in valid_days]))
            if eligible
            else None
        ),
        "mean_mvpa_minutes_per_valid_day": (
            float(np.mean([day["mvpa_minutes"] for day in valid_days]))
            if eligible
            else None
        ),
        "sedentary_share_of_wear": (
            sum(day["sedentary_minutes"] for day in valid_days) / total_wear
            if eligible and total_wear
            else None
        ),
        "mvpa_share_of_wear": (
            sum(day["mvpa_minutes"] for day in valid_days) / total_wear
            if eligible and total_wear
            else None
        ),
    }
    if person["intensity_valid_days"] != (
        person["step_complete_days"] + person["step_incomplete_days"]
    ):
        raise DataInvariantError("step day accounting failed")
    if eligible and total_wear < participant_day_threshold * 600:
        raise DataInvariantError("eligible participant wear minimum failed")
    return days, person


def participants_from_chunks(chunks):
    """Yield complete participants, rejecting nonintegral or unordered keys.

    Carry the last participant across chunk boundaries. Validate order before
    grouping, so an interleaved participant cannot be silently reassembled.
    """
    carry = None
    last_seen = None
    for chunk in chunks:
        if chunk.empty:
            continue
        keys = pd.to_numeric(chunk['SEQN'], errors='raise').to_numpy(dtype=float)
        if (not np.isfinite(keys).all() or (keys <= 0).any()
                or (keys != np.floor(keys)).any() or (np.diff(keys) < 0).any()):
            raise DataInvariantError('SEQN must be finite positive integers in sorted order')
        if last_seen is not None and keys[0] < last_seen:
            raise DataInvariantError('SEQN order decreases across chunk boundary')
        last_seen = keys[-1]
        chunk = chunk.copy()
        chunk['SEQN'] = keys.astype(np.int64)
        if carry is not None:
            chunk = pd.concat([carry, chunk], ignore_index=True)
        last = chunk['SEQN'].iloc[-1]
        carry = chunk.loc[chunk['SEQN'].eq(last)].copy()
        complete = chunk.loc[~chunk['SEQN'].eq(last)]
        for seqn, participant in complete.groupby('SEQN', sort=False, observed=True):
            yield int(seqn), participant
    if carry is not None:
        yield int(carry['SEQN'].iloc[0]), carry


def iter_xport_participants(path: Path, cycle: str, chunksize: int):
    """Read with ReadStat to preserve exact integer zeros in SAS XPORT files."""
    import pyreadstat
    columns = ['SEQN','PAXSTAT','PAXCAL','PAXDAY','PAXN','PAXHOUR','PAXMINUT','PAXINTEN']
    if cycle == '2005-2006':
        columns.append('PAXSTEP')
    elif cycle != '2003-2004':
        raise ValueError('Hip processing supports 2003-2004 and 2005-2006 only')
    if chunksize <= 0:
        raise ValueError('chunksize must be positive')
    chunks = pyreadstat.read_file_in_chunks(pyreadstat.read_xport, str(path),
                                            chunksize=chunksize, usecols=columns)
    def frames():
        for frame, _ in chunks:
            if list(frame.columns) != columns:
                raise DataInvariantError('PAX source columns do not match expected schema')
            yield frame
    for seqn, frame in participants_from_chunks(frames()):
        arrays = {c: pd.to_numeric(frame[c],errors='raise').to_numpy(dtype=float)
                  for c in columns if c != 'SEQN'}
        if 'PAXSTEP' not in arrays:
            arrays['PAXSTEP'] = np.full(len(frame),np.nan)
        yield seqn, arrays


def process_xport(path: Path, cycle: str, destination: Path, chunksize=5_000_000):
    """Process a complete raw XPORT into a new, local-only output directory.

    Reconcile source rows and participants against the documented NHANES
    release. A failed run leaves no success marker. No existing output folder
    is overwritten; reruns use a new destination.
    """
    import time
    expected = {'2003-2004': (72_250_027, 7176), '2005-2006': (74_874_095,7455)}
    if cycle not in expected:
        raise ValueError('Unsupported cycle')
    destination.mkdir(parents=True,exist_ok=False)
    start = time.monotonic()
    days, people, quality = [], [], []
    total_rows = 0
    for count,(seqn,rows) in enumerate(iter_xport_participants(path,cycle,chunksize),1):
        reasons,qc = participant_qc(rows)
        total_rows += qc['observed_rows']
        quality.append({'SEQN':seqn,'quarantined':bool(reasons),
                        'quarantine_reasons':';'.join(reasons),**qc})
        if not reasons:
            for name,calibration,minimum,threshold in VARIANTS:
                day,person = derive_variant(rows,name,calibration,minimum,threshold)
                days.extend({'SEQN':seqn,**entry} for entry in day)
                people.append({'SEQN':seqn,**person})
        if count % 500 == 0:
            print(f'{cycle}: {count} participants, {total_rows} rows',flush=True)
    if (total_rows,len(quality)) != expected[cycle]:
        raise DataInvariantError('Raw row or participant count differs from expected release')
    included = sum(not row['quarantined'] for row in quality)
    if len(days) != included*7*len(VARIANTS) or len(people) != included*len(VARIANTS):
        raise DataInvariantError('Output row accounting failed')
    for filename,records in [('pax_day_metrics',days),('pax_participant_metrics',people),('pax_qc_flags',quality)]:
        pd.DataFrame(records).to_csv(destination/(filename+'.csv.gz'),index=False,
            compression={'method':'gzip','mtime':0})
    summary={'cycle':cycle,'raw_rows':total_rows,'participants':len(quality),
             'quarantined':len(quality)-included,'elapsed_seconds':time.monotonic()-start}
    (destination/'complete.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
