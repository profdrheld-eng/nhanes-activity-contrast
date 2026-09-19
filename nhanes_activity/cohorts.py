"""Construct each cycle from public-use modules and locally processed hip data."""
from __future__ import annotations
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .questionnaires import hip_self_report, wrist_self_report, aggregate_wrist, assert_unique
from .cohort_support import derive_covariates, parse_mortality, load_pax, weighted_median, assign_discordance_group
from .transforms import weighted_midrank_normal, weighted_standardize

CYCLES = {'2003-2004':'C','2005-2006':'D','2011-2012':'G','2013-2014':'H'}
CORE = ['age','sex','race_ethnicity4','education3','PIR','smoking3']
HEALTH = ['BMI','diabetes_history','CVD_history','cancer_history']


def add_scores(cohort, domain, device_column, suffix=''):
    index = cohort.index[domain]
    weights = cohort.loc[index,'WTMEC2YR']
    s = weighted_midrank_normal(np.log1p(cohort.loc[index,'SR_day']),weights)
    d = weighted_midrank_normal(cohort.loc[index,device_column],weights)
    delta,mean_delta,sd_delta = weighted_standardize(s-d,weights)
    level,mean_level,sd_level = weighted_standardize((s+d)/2,weights)
    for name,values in [('S',s),('D',d),('Delta_raw',s-d),('Level_raw',(s+d)/2),('Delta',delta),('Level',level)]:
        cohort.loc[index,name+suffix] = values
    return {'variant':suffix or 'primary','n':len(index),'delta_raw_mean':mean_delta,
            'delta_raw_sd':sd_delta,'level_raw_mean':mean_level,'level_raw_sd':sd_level}


def build_cycle(data_dir: Path, work_dir: Path, cycle: str):
    suffix = CYCLES[cycle]
    raw = data_dir/cycle
    def read(module):
        return pd.read_sas(raw/f'{module}_{suffix}.XPT',format='xport',encoding='utf-8')
    demo,paq,smq,bmx,mcq,diq = (read(m) for m in ['DEMO','PAQ','SMQ','BMX','MCQ','DIQ'])
    assert_unique(demo,'SEQN','DEMO')
    cohort = derive_covariates(demo,smq,bmx,mcq,diq)
    cohort['cycle'] = cycle
    cohort['survey_base'] = (cohort.RIDSTATR.eq(2) & cohort.WTMEC2YR.gt(0)
        & cohort.SDMVSTRA.notna() & cohort.SDMVPSU.notna())
    cohort['adult_domain'] = cohort.survey_base & cohort.age.ge(20)
    cohort = cohort.loc[cohort.survey_base].copy()
    is_hip = suffix in 'CD'
    if is_hip:
        adult_ids = set(cohort.loc[cohort.adult_domain,'SEQN'].astype(int))
        sr = hip_self_report(paq,read('PAQIAF'),adult_ids)
        device = load_pax(work_dir/cycle/'pax/pax_participant_metrics.csv.gz')
    else:
        sr = wrist_self_report(paq)
        device = aggregate_wrist(read('PAXDAY'))
    for frame in [sr,device,parse_mortality(raw/f'NHANES_{cycle.replace("-","_")}_MORT_2019_PUBLIC.dat')]:
        cohort = cohort.merge(frame,on='SEQN',how='left',validate='one_to_one')
    for name in ['SR_complete','PAX_primary_eligible'] + ([] if is_hip else ['PAX_one_day_eligible']):
        cohort[name] = cohort[name].astype('boolean').fillna(False).astype(bool)
    cohort['mortality_followup_valid'] = (cohort.ELIGSTAT.eq(1) & cohort.MORTSTAT.isin([0,1]) & cohort.PERMTH_EXM.gt(0)).fillna(False)
    cohort['follow_up_years'] = cohort.PERMTH_EXM/12.0
    common = cohort.adult_domain & cohort.mortality_followup_valid & cohort.SR_complete
    if is_hip:
        metric = 'PAX_log_counts_primary'
    else:
        cohort['PAX_log_mean_daily_mims'] = np.log1p(cohort.PAX_mean_daily_mims)
        cohort['PAX_log_mims_per_valid_minute'] = np.log1p(cohort.PAX_mims_per_valid_minute)
        metric = 'PAX_log_mean_daily_mims'
    cohort['reference_domain'] = common & cohort.PAX_primary_eligible & cohort[metric].notna()
    contrasts = [add_scores(cohort,cohort.reference_domain,metric)]
    cohort['primary_complete_case'] = cohort.reference_domain & cohort[CORE].notna().all(axis=1)
    if is_hip:
        cohort['model_c_complete_case'] = cohort.primary_complete_case & cohort[HEALTH].notna().all(axis=1)
        ref = cohort.reference_domain
        sm = weighted_median(cohort.loc[ref,'S'],cohort.loc[ref,'WTMEC2YR'])
        dm = weighted_median(cohort.loc[ref,'D'],cohort.loc[ref,'WTMEC2YR'])
        cohort.loc[ref,'discordance_group'] = [assign_discordance_group(s,d,sm,dm)
            for s,d in zip(cohort.loc[ref,'S'],cohort.loc[ref,'D'])]
    else:
        for variant,eligible,metric in [('intensity','PAX_primary_eligible','PAX_log_mims_per_valid_minute'),
                                      ('one_day','PAX_one_day_eligible','PAX_log_mean_daily_mims')]:
            domain = common & cohort[eligible] & cohort[metric].notna()
            cohort['reference_domain_'+variant] = domain
            contrasts.append(add_scores(cohort,domain,metric,'_'+variant))
            cohort[variant+'_complete_case'] = domain & cohort[CORE].notna().all(axis=1)
    return cohort,contrasts


def write_cycle(data_dir: Path, work_dir: Path, cycle: str):
    out = work_dir/cycle/'cohort'
    if out.exists():
        raise FileExistsError(f'Refusing to overwrite completed or partial cohort: {out}')
    cohort,contrasts = build_cycle(data_dir,work_dir,cycle)
    out.mkdir(parents=True,exist_ok=False)
    cohort.to_csv(out/'analytic_cohort.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    selected = cohort.loc[cohort.primary_complete_case]
    summary = {'cycle':cycle,'survey_base_n':len(cohort),'reference_n':int(cohort.reference_domain.sum()),
        'complete_case_n':len(selected),'deaths':int(selected.MORTSTAT.sum()),'contrasts':contrasts}
    (out/'complete.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
