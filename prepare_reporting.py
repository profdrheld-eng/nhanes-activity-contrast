#!/usr/bin/env python3
"""Export nonsequential completeness, unweighted follow-up and exact quintile boundaries.

Uses the existing frozen model input; does not modify scores, raw data or samples.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from nhanes_activity.reporting import missing_counts, quintile_cuts, assign_quintiles, followup_summary, composite
from nhanes_activity.transforms import weighted_midrank_normal, weighted_standardize

CYCLES=['2003-2004','2005-2006','2011-2012','2013-2014']
VARIABLES=['age','sex','race_ethnicity4','education3','PIR','smoking3','BMI',
           'CVD_history','cancer_history','diabetes_history','disease','cvd_cancer',
           'ELIGSTAT','MORTSTAT','PERMTH_EXM','SR_day','device_raw','device_log','valid_days',
           'S','D','Delta','Level','SR_equal','Delta_equal','Level_equal',
           'Delta_zlog','Level_zlog','Delta_zraw','Level_zraw','quintile']
VARIANT_VARIABLES={**{f'PAX_log1p_counts_per_wear_minute__{v}':'hip' for v in ['valid_days_3','valid_days_5','calibration_1_or_2','nonwear_90']}, **{v:'wrist' for v in ['Delta_intensity','Level_intensity','Delta_one_day','Level_one_day']}}
DOMAINS=['adult_domain','reference_domain','primary_complete_case','age_domain']


def build(work):
    source=work/'models/revision_cohort.csv.gz';out=work/'reporting'
    if out.exists():raise FileExistsError('Refusing to overwrite reporting stage')
    d=pd.read_csv(source,float_precision='round_trip',low_memory=False)
    if set(d.cycle)!=set(CYCLES) or d.duplicated(['cycle','SEQN']).any():raise ValueError('Expected four cycles and unique person keys')
    for v in ['adult_domain','reference_domain','primary_complete_case','SR_complete','PAX_primary_eligible','mortality_followup_valid']:
        if not d[v].astype(str).str.lower().isin(['true','false','1','0']).all():raise ValueError('Invalid domain flag: '+v)
        d[v]=d[v].astype(str).str.lower().isin(['true','1'])
    d['disease']=composite(d,['CVD_history','cancer_history','diabetes_history'])
    d['cvd_cancer']=composite(d,['CVD_history','cancer_history'])
    d['age_domain']=d.primary_complete_case & np.isfinite(d.BMI) & np.where(d.era.eq('hip'),d.disease.notna(),d.cvd_cancer.notna())
    completeness=[];gates=[];bounds=[];qcounts=[];times=[]
    for cycle in CYCLES:
        x=d.loc[d.cycle.eq(cycle)]
        for domain in DOMAINS:
            sample=x.loc[x[domain]]
            applicable=VARIABLES+[v for v,era in VARIANT_VARIABLES.items() if era==('hip' if cycle in CYCLES[:2] else 'wrist')]
            for row in missing_counts(sample,applicable):
                completeness.append({'cycle':cycle,'domain':domain,**row})
            # False flags are criterion failures, not absent values. Counts overlap.
            criteria={'mortality linkage eligible':sample.ELIGSTAT.eq(1),
                      'known vital status':sample.MORTSTAT.isin([0,1]),
                      'positive examination follow-up':sample.PERMTH_EXM.gt(0),
                      'complete self-report':sample.SR_complete,
                      'at least four valid device days':sample.valid_days.ge(4),
                      'eligible device score':sample.PAX_primary_eligible & sample.device_log.notna(),
                      'complete principal covariates':sample[['age','sex','race_ethnicity4','education3','PIR','smoking3']].notna().all(axis=1)}
            for label,valid in criteria.items():gates.append({'cycle':cycle,'domain':domain,'criterion':label,'denominator':len(sample),'failed':int((~valid).sum())})
        ref=x.loc[x.reference_domain];w=ref.WTMEC2YR.to_numpy()
        # Match the exact arithmetic in prepare_models.py, not rounded display values.
        s=weighted_midrank_normal(np.log1p(ref.SR_day),w)
        dev=weighted_midrank_normal(ref.device_log,w)
        delta=weighted_standardize(s-dev,w)[0]
        if not np.allclose(delta,ref.Delta,rtol=0,atol=1e-12):raise ValueError('Contrast reconstruction differs')
        cuts=quintile_cuts(delta,w);assigned=assign_quintiles(delta,cuts)
        if not np.array_equal(assigned,ref.quintile.to_numpy()):raise ValueError('Quintile membership differs from original models')
        bounds.append({'cycle':cycle,'reference_n':len(ref),**dict(zip(['p20','p40','p60','p80'],cuts))})
        for domain in ['reference_domain','primary_complete_case']:
            a=x.loc[x[domain]]
            for quintile in range(1,6):qcounts.append({'cycle':cycle,'domain':domain,'quintile':quintile,'n':int(a.quintile.eq(quintile).sum()),'denominator':len(a)})
    for name,mask in [(c,d.cycle.eq(c)) for c in CYCLES]+[('hip',d.era.eq('hip')),('wrist',d.era.eq('wrist'))]:
        times.append({'group':name,**followup_summary(d.loc[mask & d.primary_complete_case])})
    out.mkdir()
    for name,rows in [('missingness',completeness),('criterion_failures',gates),('quintile_cutpoints',bounds),('quintile_counts',qcounts),('followup',times)]:
        pd.DataFrame(rows).to_csv(out/(name+'.csv'),index=False)
    (out/'descriptive_manifest.json').write_text(json.dumps({'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'input':'models/revision_cohort.csv.gz','variables':VARIABLES,'variant_variables':VARIANT_VARIABLES,'domains':DOMAINS,'counts':'unweighted; nonsequential and potentially overlapping','cutpoints':'cycle reference sample; weighted empirical inverse CDF; ties at boundary assigned to lower quintile'},indent=2)+'\n')
    print(f'Exported {len(completeness)} completeness cells, {len(gates)} criterion cells, four exact cutpoint sets and six follow-up summaries')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work-dir',type=Path,required=True);a=p.parse_args()
    work=a.work_dir.resolve();root=Path(__file__).resolve().parent
    if work==root or root in work.parents:p.error('Work must be outside repository')
    build(work)
