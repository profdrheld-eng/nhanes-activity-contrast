#!/usr/bin/env python3
"""Format STROBE reporting additions from machine-readable descriptive and model outputs."""
import argparse
import csv
import json
from pathlib import Path
import pandas as pd
import numpy as np
from prepare_reporting import CYCLES, VARIABLES, VARIANT_VARIABLES, DOMAINS

LABELS={'race_ethnicity4':'Race and ethnicity','education3':'Education','PIR':'Poverty-income ratio',
 'smoking3':'Smoking','disease':'Any CVD, cancer or diabetes','cvd_cancer':'CVD or cancer',
 'ELIGSTAT':'Linkage eligibility code','MORTSTAT':'Vital status','PERMTH_EXM':'Examination follow-up months',
 'SR_day':'Self-report activity index','device_raw':'Primary device measure','device_log':'Log primary device measure',
 'valid_days':'Number of valid device days','S':'Self-report normal rank','D':'Device normal rank',
 'Delta':'Standardized contrast','Level':'Standardized average','SR_equal':'Equal-weight self-report index',
 'Delta_equal':'Contrast, equal vigorous weighting','Level_equal':'Average, equal vigorous weighting',
 'Delta_zlog':'Contrast, z scores of log values','Level_zlog':'Average, z scores of log values',
 'Delta_zraw':'Contrast, z scores of raw values','Level_zraw':'Average, z scores of raw values',
 'quintile':'Contrast quintile','CVD_history':'CVD history','cancer_history':'Cancer history',
 'diabetes_history':'Diabetes history','BMI':'BMI','age':'Age','sex':'Sex',
 'Delta_intensity':'Contrast, wrist intensity','Level_intensity':'Average, wrist intensity',
 'Delta_one_day':'Contrast, wrist one-day rule','Level_one_day':'Average, wrist one-day rule'}
for v in VARIANT_VARIABLES:
 if v.startswith('PAX_'):LABELS[v]='Hip log device, '+v.split('__')[1].replace('_',' ')
GROUPS=CYCLES+['hip','wrist']
def group_label(x):return {'hip':'2003–2006 pooled','wrist':'2011–2014 pooled'}.get(x,x)
def interval(r):return f'{r.estimate:.3f} [{r.ci_low:.3f}, {r.ci_high:.3f}]'
def pvalue(x):return '< 0.001' if x<.001 else f'{x:.3f}'

def build(work):
 r=work/'reporting';tables={};notes={};titles={}
 m=pd.read_csv(r/'missingness.csv');g=pd.read_csv(r/'criterion_failures.csv')
 if m.duplicated(['cycle','domain','variable']).any() or not (m.observed+m.missing==m.denominator).all():raise ValueError('Invalid completeness cells')
 labels=['A. Examined adult survey domain','B. Paired-activity reference sample','C. Principal mortality sample','D. Age-analysis sample']
 rows=[['Variable or domain','2003–2004','2005–2006','2011–2012','2013–2014']]
 for domain,label in zip(DOMAINS,labels):
  x=m.loc[m.domain.eq(domain)];ns=[int(x.loc[x.cycle.eq(c),'denominator'].iloc[0]) for c in CYCLES]
  rows.append([label]+[f'N = {n:,}' for n in ns])
  for variable in VARIABLES+list(VARIANT_VARIABLES):
   values=[]
   for c in CYCLES:
    a=x.loc[x.cycle.eq(c)&x.variable.eq(variable)]
    if len(a)==0:
     if variable not in VARIANT_VARIABLES:raise ValueError('Missing completeness cell')
     values.append('Not applicable');continue
    if len(a)!=1:raise ValueError('Duplicate completeness cell')
    v=a.iloc[0];values.append(f'{int(v.missing):,} ({100*v.missing/v.denominator:.1f}%)')
   rows.append([LABELS[variable]]+values)
 tables['Table S32']=rows;titles['Table S32']='Nonsequential missingness by survey cycle and analysis domain'
 notes['Table S32']='Values are unweighted missing counts (percent of the stated domain N), evaluated separately for each variable; counts overlap and must not be summed. The examined adult survey domain requires age ≥20 years, examination, positive examination weight and complete survey identifiers; these design fields have no missing values by construction. Analysis domains are nested within this domain. Missing constructed rank scores outside their eligible reference domain mean undefined scores, not necessarily questionnaire nonresponse. Device metrics may be observed despite failing the required number of valid days; criterion failures are reported separately in Table S33. Not applicable identifies a measurement variant not used in that period. BMI, body mass index; CVD, cardiovascular disease; z score, standardized value. Any-disease composites are positive if any component is positive, negative only if all are negative, and missing otherwise. No values were imputed.'
 rows=[['Criterion not met','2003–2004','2005–2006','2011–2012','2013–2014']]
 for domain,label in zip(DOMAINS,labels):
  x=g.loc[g.domain.eq(domain)];rows.append([label]+[f'N = {int(x.loc[x.cycle.eq(c),"denominator"].iloc[0]):,}' for c in CYCLES])
  for criterion in x.criterion.unique():
   values=[]
   for c in CYCLES:
    a=x.loc[x.cycle.eq(c)&x.criterion.eq(criterion)];assert len(a)==1;v=a.iloc[0]
    values.append(f'{int(v.failed):,} ({100*v.failed/v.denominator:.1f}%)')
   rows.append([criterion]+values)
 tables['Table S33']=rows;titles['Table S33']='Overlapping eligibility and completeness failures'
 notes['Table S33']='Values are unweighted numbers failing each criterion (percent of the stated domain N). Failures are assessed independently and can overlap; they differ from the sequential exclusions in Table S2. Missing or invalid values fail the relevant criterion. Principal covariates are age, sex, race and ethnicity, education, poverty-income ratio and smoking. A false eligibility flag is a criterion failure, not a missing flag.'
 f=pd.read_csv(r/'followup.csv');rows=[['Sample','n / deaths','Person-years','Mean ± SD, years','Range, years']]
 for name in GROUPS:
  x=f.loc[f.group.eq(name)];assert len(x)==1;v=x.iloc[0]
  rows.append([group_label(name),f'{int(v.n):,} / {int(v.deaths):,}',f'{v.person_years:,.2f}',f'{v.mean_years:.2f} ± {v.sd_years:.2f}',f'{v.minimum_years:.2f}–{v.maximum_years:.2f}'])
 tables['Table S34']=rows;titles['Table S34']='Unweighted mortality follow-up in the principal analysis samples'
 notes['Table S34']='All summaries are unweighted and use exactly the principal mortality complete-case samples. Person-years sum examination follow-up months divided by 12; arithmetic mean and sample standard deviation (SD) describe the observed follow-up. Pooled totals sum the constituent cycles. These means differ from the survey-weighted means in Table 1. Public-use follow-up times may be perturbed for confidentiality.'
 b=pd.read_csv(r/'quintile_cutpoints.csv');rows=[['Survey cycle','Reference n','20th percentile','40th percentile','60th percentile','80th percentile']]
 for c in CYCLES:
  x=b.loc[b.cycle.eq(c)];assert len(x)==1;v=x.iloc[0];rows.append([c,f'{int(v.reference_n):,}']+[f'{v[k]:.6f}' for k in ['p20','p40','p60','p80']])
 tables['Table S35']=rows;titles['Table S35']='Cycle-specific weighted cutpoints for contrast quintiles'
 notes['Table S35']='Cutpoints use the standardized signed contrast in each paired-activity reference sample and the two-year examination weight WTMEC2YR. The empirical weighted quantile is the smallest observed value whose cumulative weight share reaches the target. Q1 contains values ≤c20, Q2 values >c20 to ≤c40, Q3 values >c40 to ≤c60, Q4 values >c60 to ≤c80, and Q5 values >c80; Q3 is the model reference. Ties at a boundary enter the lower quintile. Displayed cutpoints are rounded to six decimals; full-precision cutpoints and verified memberships are retained in the reporting CSV exports. Categories were not recalculated in complete-case samples.'
 a=pd.read_csv(r/'mortality_adjustment.csv');a=a.loc[a.inference.eq('full_design_t')]
 if len(a)!=18 or a.duplicated(['group','adjustment']).any():raise ValueError('Incomplete adjustment models')
 if not np.isfinite(a[['estimate','ci_low','ci_high','p','n','deaths']]).all().all():raise ValueError('Nonfinite estimates')
 rows=[['Sample','Adjustment','n / deaths','HR [95% CI]','p']]
 for group in GROUPS:
  x=a.loc[a.group.eq(group)];assert x.n.nunique()==1 and x.deaths.nunique()==1
  for adjustment,label in [('unadjusted','Unadjusted'),('level_only','Average component only'),('principal','Principal adjustment')]:
   v=x.loc[x.adjustment.eq(adjustment)].iloc[0];rows.append([group_label(group),label,f'{int(v.n):,} / {int(v.deaths):,}',interval(v),pvalue(v.p)])
 tables['Table S36']=rows;titles['Table S36']='Mortality associations before and after adjustment in identical samples'
 notes['Table S36']='All models use identical principal complete-case participants and deaths within each sample. HR, hazard ratio per one within-cycle standard deviation higher signed contrast; CI, confidence interval. Unadjusted models contain the contrast only; average-component models additionally contain the standardized average activity component; principal models additionally adjust for age, sex, race and ethnicity, education, poverty-income ratio and smoking. All models retain survey weights, strata and primary sampling units; pooled models retain cycle-specific baseline hazards, including the unadjusted model. Efron ties and design-based t inference were used (degrees of freedom: 15, 15, 17 and 15 for successive cycles; 30 and 32 for pooled periods). Confidence intervals and p values are not multiplicity-adjusted. Unadjusted and conditional coefficients represent different associations; changes with adjustment do not establish confounding or mediation. These reporting comparisons were added post hoc.'
 return tables,titles,notes

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work-dir',type=Path,required=True);args=p.parse_args()
 tables,titles,notes=build(args.work_dir);out=args.work_dir/'displays';out.mkdir(exist_ok=True)
 target=out/'reporting_tables.json'
 if target.exists() or any((out/(k.lower().replace(' ','_')+'.csv')).exists() for k in tables):raise FileExistsError('Refusing to overwrite reporting tables')
 target.write_text(json.dumps(tables,ensure_ascii=False,indent=2)+'\n')
 (out/'reporting_table_notes.json').write_text(json.dumps({'titles':titles,'notes':notes},ensure_ascii=False,indent=2)+'\n')
 for name,rows in tables.items():
  with (out/(name.lower().replace(' ','_')+'.csv')).open('x',newline='') as handle:csv.writer(handle).writerows(rows)
 print('Exported Tables S32–S36')
