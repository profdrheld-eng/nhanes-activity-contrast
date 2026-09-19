#!/usr/bin/env python3
"""Build model inputs and descriptive exports from newly constructed cohorts.

Must be run after all four cycles. The current score and analysis definitions
are retained; raw questionnaire rescoring verifies complete-case inputs.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from nhanes_activity.questionnaires import hip_self_report, wrist_self_report
from nhanes_activity.transforms import weighted_midrank_normal, weighted_standardize

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-dir',type=Path,required=True)
parser.add_argument('--work-dir',type=Path,required=True)
args=parser.parse_args()
WORK=args.work_dir.resolve(); DATA=args.data_dir.resolve()
ROOT=Path(__file__).resolve().parent
if WORK==ROOT or ROOT in WORK.parents:
    parser.error('Work directory must be outside the code repository')
if WORK==DATA or DATA in WORK.parents or WORK in DATA.parents:
    parser.error('Work and raw-data directories must be disjoint')
OUT=WORK/'models'
OUT.mkdir(parents=True,exist_ok=False)
hip=SimpleNamespace(derive_self_report=hip_self_report,weighted_midrank_normal=weighted_midrank_normal,weighted_standardize=weighted_standardize)
wrist=SimpleNamespace(derive_self_report=wrist_self_report)

def flag(x):
    return x.astype(str).str.lower().isin(['true','1','1.0'])

def mean(x,w):
    return np.average(np.asarray(x,float),weights=np.asarray(w,float))

def sd(x,w):
    return np.sqrt(mean((np.asarray(x)-mean(x,w))**2,w))

def corr(x,y,w):
    return mean((x-mean(x,w))*(y-mean(y,w)),w)/sd(x,w)/sd(y,w)

def quant(x,w,p):
    order=np.argsort(x,kind='stable'); x=np.asarray(x)[order]; w=np.asarray(w)[order]
    return x[np.minimum(np.searchsorted(np.cumsum(w)/sum(w),p),len(x)-1)]

flows=[]; distributions=[]; correlations=[]; audits=[]; included=[]; modules=[]; frames=[]
cohorts={'2003-2004':WORK/'2003-2004/cohort/analytic_cohort.csv.gz',
         '2005-2006':WORK/'2005-2006/cohort/analytic_cohort.csv.gz',
         '2011-2012':WORK/'2011-2012/cohort/analytic_cohort.csv.gz',
         '2013-2014':WORK/'2013-2014/cohort/analytic_cohort.csv.gz'}
for cycle,suffix in [('2003-2004','C'),('2005-2006','D'),('2011-2012','G'),('2013-2014','H')]:
    raw=DATA/cycle
    demo=pd.read_sas(raw/f'DEMO_{suffix}.XPT',format='xport')
    paq=pd.read_sas(raw/f'PAQ_{suffix}.XPT',format='xport')
    d=pd.read_csv(cohorts[cycle],float_precision='round_trip',low_memory=False)
    if 'cycle' in d: d=d.loc[d.cycle==cycle].copy()
    d['cycle']=cycle; d['era']='hip' if suffix in 'CD' else 'wrist'
    for col in ['adult_domain','reference_domain','primary_complete_case','SR_complete','PAX_primary_eligible','mortality_followup_valid']:
        d[col]=flag(d[col])
    adult_ids=set(d.loc[d.adult_domain,'SEQN'])
    if suffix in 'CD':
        paqiaf=pd.read_sas(raw/f'PAQIAF_{suffix}.XPT',format='xport')
        sr=hip.derive_self_report(paq,paqiaf,adult_ids)
        equal=sr[['SR_transport_30d','SR_home_30d','SR_moderate_30d','SR_vigorous_30d']].sum(axis=1,min_count=4)/30
        d['SR_equal']=d.SEQN.map(pd.Series(equal.to_numpy(),index=sr.SEQN))
        d['device_raw']=d['PAX_counts_per_wear_minute__primary']
        d['device_log']=d['PAX_log_counts_primary']
        d['valid_days']=d['PAX_intensity_valid_days__primary']
        required=['PAD020','PAQ050Q','PAQ050U','PAD080','PAQ100','PAD120','PAD160','PAD320','PAD200']
    else:
        sr=wrist.derive_self_report(paq)
        cols=[c for c in sr if c.startswith('SR_') and c.endswith('_weekly_equivalent_minutes')]
        plain=sr[cols].copy()
        for col in plain:
            if 'vigorous' in col: plain[col]/=2
        equal=plain.sum(axis=1,min_count=5)/7
        d['SR_equal']=d.SEQN.map(pd.Series(equal.to_numpy(),index=sr.SEQN))
        d['device_raw']=d['PAX_mean_daily_mims'];d['device_log']=d['PAX_log_mean_daily_mims']
        d['valid_days']=d['PAX_valid_days']
        required=['PAQ605','PAQ610','PAD615','PAQ620','PAQ625','PAD630','PAQ635','PAQ640','PAD645','PAQ650','PAQ655','PAD660','PAQ665','PAQ670','PAD675']
    check=d[['SEQN','SR_day','SR_complete']].merge(sr[['SEQN','SR_day','SR_complete']],on='SEQN',how='left',suffixes=('_old','_new'))
    mask=check.SR_complete_old
    error=np.max(np.abs(check.loc[mask,'SR_day_old']-check.loc[mask,'SR_day_new']))
    # Historical deterministic CSVs use 12 significant digits; allow their
    # serialization rounding (largest observed discrepancy is ~5e-9 minutes).
    assert error<1e-7, (cycle, error, int(check.loc[mask,'SR_day_new'].isna().sum()))
    assert (flag(check.SR_complete_new)==check.SR_complete_old).all()
    audits.append({'cycle':cycle,'check':'raw_self_report_reproduction','max_abs_error':error,'passed':True})
    d['SR_day']=d.SEQN.map(sr.set_index('SEQN').SR_day)
    for col in required:
        modules.append({'cycle':cycle,'module':'PAQ','variable':col,'present':col in paq,'nonmissing':int(paq[col].notna().sum())})
    if suffix in 'CD':
        for col in ['PADACTIV','PADLEVEL','PADTIMES','PADDURAT']:
            modules.append({'cycle':cycle,'module':'PAQIAF','variable':col,'present':col in paqiaf,'nonmissing':int(paqiaf[col].notna().sum())})
    # Recover the full examined-adult denominator directly from DEMO.
    base=demo.loc[(demo.RIDAGEYR>=20)&(demo.RIDSTATR==2)].copy()
    cur=base.copy();prev=len(cur)
    flows.append({'cycle':cycle,'stage':0,'criterion':'examined adults age >=20','remaining':prev,'excluded':0})
    cur=cur.loc[cur.WTMEC2YR>0];flows.append({'cycle':cycle,'stage':1,'criterion':'positive examination weight','remaining':len(cur),'excluded':prev-len(cur)});prev=len(cur)
    cur=cur.loc[cur.SDMVSTRA.notna()&cur.SDMVPSU.notna()];flows.append({'cycle':cycle,'stage':2,'criterion':'complete survey identifiers','remaining':len(cur),'excluded':prev-len(cur)});prev=len(cur)
    cur=d.loc[d.SEQN.isin(cur.SEQN)].copy();assert len(cur)==prev
    steps=[('mortality linkage eligible',cur.ELIGSTAT.eq(1)),
           ('valid vital status',cur.MORTSTAT.isin([0,1])),
           ('positive nonmissing examination follow-up',cur.PERMTH_EXM.gt(0)),
           ('complete self-report',cur.SR_complete),
           ('at least four valid device days',cur.valid_days.ge(4)),
           ('valid device score and eligibility',cur.PAX_primary_eligible & cur.device_log.notna()),
           ('complete model covariates',cur.primary_complete_case)]
    for stage,(label,keep) in enumerate(steps,3):
        cur=cur.loc[keep.reindex(cur.index)];n=len(cur)
        flows.append({'cycle':cycle,'stage':stage,'criterion':label,'remaining':n,'excluded':prev-n});prev=n
    assert set(cur.SEQN)==set(d.loc[d.primary_complete_case,'SEQN'])
    flows.append({'cycle':cycle,'stage':10,'criterion':'deaths in final sample','remaining':int(cur.MORTSTAT.sum()),'excluded':np.nan})
    flows.append({'cycle':cycle,'stage':11,'criterion':'paired descriptive reference sample','remaining':int(d.reference_domain.sum()),'excluded':np.nan})
    ref=d.reference_domain; idx=d.index[ref]; w=d.loc[ref,'WTMEC2YR'].to_numpy()
    s=hip.weighted_midrank_normal(np.log1p(d.loc[ref,'SR_day']),w)
    dev=hip.weighted_midrank_normal(d.loc[ref,'device_log'],w)
    delta,mu,sigma=hip.weighted_standardize(s-dev,w)
    assert np.allclose(delta,d.loc[ref,'Delta'],atol=1e-9), (cycle,float(np.max(np.abs(delta-d.loc[ref,'Delta']))))
    for kind,sv,dv in [('equal',hip.weighted_midrank_normal(np.log1p(d.loc[ref,'SR_equal']),w),dev),
                       ('zlog',hip.weighted_standardize(np.log1p(d.loc[ref,'SR_day']),w)[0],hip.weighted_standardize(d.loc[ref,'device_log'],w)[0]),
                       ('zraw',hip.weighted_standardize(d.loc[ref,'SR_day'].to_numpy(),w)[0],hip.weighted_standardize(d.loc[ref,'device_raw'].to_numpy(),w)[0])]:
        d.loc[idx,f'Delta_{kind}']=hip.weighted_standardize(sv-dv,w)[0]
        d.loc[idx,f'Level_{kind}']=hip.weighted_standardize((sv+dv)/2,w)[0]
    cuts=quant(delta,w,[.2,.4,.6,.8]);d.loc[idx,'quintile']=np.searchsorted(cuts,delta,side='left')+1
    for sample,keep in [('reference',d.reference_domain),('complete_case',d.primary_complete_case)]:
        x=d.loc[keep];ww=x.WTMEC2YR.to_numpy();ones=np.ones(len(x))
        for name,a,b in [('raw_SR_device',x.SR_day,x.device_raw),('log_SR_device',np.log1p(x.SR_day),x.device_log),('normal_rank_SR_device',x.S,x.D),('average_difference',x.Level,x.Delta)]:
            correlations.append({'cycle':cycle,'sample':sample,'pair':name,'n':len(x),'weighted_r':corr(a,b,ww),'unweighted_r':corr(a,b,ones)})
    for weighted,ww in [('weighted',w),('unweighted',np.ones(len(w)))]:
        raw_delta=s-dev;qs=quant(raw_delta,ww,[.01,.05,.1,.25,.5,.75,.9,.95,.99])
        row={'cycle':cycle,'weighting':weighted,'n':len(w),'mean':mean(raw_delta,ww),'sd':sd(raw_delta,ww),'min':min(raw_delta),'max':max(raw_delta)}
        row.update(dict(zip(['p01','p05','p10','p25','median','p75','p90','p95','p99'],qs)));distributions.append(row)
    audits.append({'cycle':cycle,'check':'rank_contrast_reproduction','max_abs_error':float(np.max(np.abs(delta-d.loc[ref,'Delta']))),'passed':True})
    # Available-variable review covers files already available, without new downloads.
    for name in ['MCQ','DIQ','PAQ','SMQ','BMX']:
        mod=pd.read_sas(raw/f'{name}_{suffix}.XPT',format='xport')
        for c in mod:
            if c!='SEQN': modules.append({'cycle':cycle,'module':name,'variable':c,'present':True,'nonmissing':int(mod[c].notna().sum())})
    # Descriptive selection comparison in the examined, positive-weight adult domain.
    for label,keep in [('included',d.adult_domain&d.primary_complete_case),('excluded',d.adult_domain&~d.primary_complete_case)]:
        x=d.loc[keep]
        for var in ['age','BMI','PIR','SR_day']:
            a=x[var];valid=a.notna();ww=x.loc[valid,'WTMEC2YR']
            included.append({'cycle':cycle,'group':label,'variable':var,'n_group':len(x),'n_observed':int(valid.sum()),'mean_or_proportion':mean(a[valid],ww),'sd':sd(a[valid],ww)})
        for var,value in [('sex','female'),('smoking3','current'),('CVD_history','yes'),('cancer_history','yes'),('diabetes_history','yes')]:
            valid=x[var].notna();ww=x.loc[valid,'WTMEC2YR']
            included.append({'cycle':cycle,'group':label,'variable':var+'='+value,'n_group':len(x),'n_observed':int(valid.sum()),'mean_or_proportion':mean(x.loc[valid,var].eq(value).astype(float),ww),'sd':np.nan})
    frames.append(d)

all_data=pd.concat(frames,ignore_index=True)
all_data.to_csv(OUT/'revision_cohort.csv.gz',index=False,compression={'method':'gzip','mtime':0})
for name,rows in [('participant_flow',flows),('contrast_distribution',distributions),('correlations',correlations),('source_checks',audits),('inclusion_comparison',included),('module_variables',modules)]:
    pd.DataFrame(rows).drop_duplicates().to_csv(OUT/f'{name}.csv',index=False)
inputs=list(cohorts.values())
records=[{'path':str(p.relative_to(WORK)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in dict.fromkeys(inputs)]
(OUT/'input_hashes.json').write_text(json.dumps(records,indent=2)+'\n')
print(pd.DataFrame(flows).pivot(index='criterion',columns='cycle',values='remaining').to_string())
print('All original score and cohort checks passed')
