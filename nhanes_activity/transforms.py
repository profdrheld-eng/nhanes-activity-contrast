"""Cycle-specific, examination-weighted transformations used in the paper.

Midranks pool the weights of exactly equal input values. Percentiles are
bounded at 0.0001 and 0.9999 before the inverse normal transformation.
Standardization uses a population-weighted SD, not a sample correction or SE.
The caller must select the paired reference domain before calling these
functions; missing observations are rejected rather than silently excluded.
"""
from statistics import NormalDist
import numpy as np


def _vectors(values, weights):
    x = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    if x.ndim != 1 or w.ndim != 1 or x.shape != w.shape or x.size == 0:
        raise ValueError('Values and weights must be nonempty matching vectors')
    if not np.isfinite(x).all() or not np.isfinite(w).all() or (w <= 0).any():
        raise ValueError('Values must be finite and weights finite and positive')
    if not np.isfinite(w.sum()):
        raise ValueError('Sum of weights must be finite')
    return x, w


def weighted_midrank_normal(values, weights):
    """Return inverse-normal weighted midranks in the original row order."""
    x, w = _vectors(values, weights)
    order = np.argsort(x, kind='mergesort')
    xs, ws = x[order], w[order]
    result_sorted = np.empty(len(x), dtype=float)
    total = ws.sum()
    cumulative = 0.0
    start = 0
    while start < len(xs):
        end = start + 1
        while end < len(xs) and xs[end] == xs[start]:
            end += 1
        tie_weight = ws[start:end].sum()
        percentile = (cumulative + 0.5 * tie_weight) / total
        percentile = min(0.9999, max(0.0001, float(percentile)))
        result_sorted[start:end] = NormalDist().inv_cdf(percentile)
        cumulative += tie_weight
        start = end
    result = np.empty(len(x), dtype=float)
    result[order] = result_sorted
    return result


def weighted_standardize(values, weights):
    """Return standardized values, weighted mean, and population-weighted SD."""
    x, w = _vectors(values, weights)
    mean = float(np.sum(w * x) / np.sum(w))
    sd = float(np.sqrt(np.sum(w * (x - mean) ** 2) / np.sum(w)))
    if not np.isfinite(mean) or not np.isfinite(sd) or sd <= 0:
        raise ValueError('Weighted mean must be finite and SD positive and finite')
    return (x - mean) / sd, mean, sd
