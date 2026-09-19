"""Descriptive reporting helpers; no imputation or model-specific re-selection."""
import numpy as np
import pandas as pd


def missing_counts(frame, variables):
    return [{'variable':v,'denominator':len(frame),'observed':int(frame[v].notna().sum()),
             'missing':int(frame[v].isna().sum())} for v in variables]


def quintile_cuts(values, weights):
    x=np.asarray(values,float);w=np.asarray(weights,float)
    if x.ndim!=1 or x.shape!=w.shape or not len(x) or not np.isfinite(x).all() or not np.isfinite(w).all() or not (w>0).all():
        raise ValueError('Quantiles require finite values and positive finite weights')
    order=np.argsort(x,kind='stable');x=x[order];w=w[order]
    return x[np.searchsorted(np.cumsum(w)/sum(w),[.2,.4,.6,.8],side='left')]


def assign_quintiles(values, cuts):
    return np.searchsorted(cuts,np.asarray(values,float),side='left')+1


def followup_summary(frame):
    t=frame.follow_up_years.to_numpy(float);death=frame.MORTSTAT.to_numpy(float)
    if not len(t) or not np.isfinite(t).all() or not (t>0).all() or not np.isin(death,[0,1]).all():
        raise ValueError('Follow-up summary requires positive times and known vital status')
    return {'n':len(t),'deaths':int(death.sum()),'person_years':float(t.sum()),
            'mean_years':float(t.mean()),'sd_years':float(t.std(ddof=1)),
            'minimum_years':float(t.min()),'maximum_years':float(t.max())}


def composite(frame, variables):
    x=frame[variables];out=pd.Series(pd.NA,index=frame.index,dtype='string')
    out.loc[x.eq('no').all(axis=1)]='no'
    out.loc[x.eq('yes').any(axis=1)]='yes'
    return out
