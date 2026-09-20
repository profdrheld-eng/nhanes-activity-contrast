#!/usr/bin/env python3
"""Export both S40 panels from completed age-comparison models, without refitting."""
import argparse
import csv
import json
import math
from pathlib import Path
from export_reporting_diagnostics import read, index, number

MODELS={
    'historical_under80':'Original domains; ages 20–79',
    'harmonized_under80':'Harmonized domains; ages 20–79',
    'harmonized_cap80':'Harmonized domains; common 80+ category',
}
INFERENCES={'full_design':'Full design','residual':'Residual','asymptotic':'Asymptotic'}
TITLE='Additional exploratory comparisons of age relationships between measurement periods'

NOTE='Note. The outcome was the standardized signed self-report–device rank contrast. Panel A jointly tests the linear and quadratic age-by-period interactions in fully period-interacted covariate models with cycle intercepts. The asymptotic test uses chi-square=2×F with two degrees of freedom. Panel B presents age-related differences from age 50, not absolute predicted outcome levels or hazard ratios. Full-design confidence intervals and p values use 62 denominator degrees of freedom; residual alternatives use 39, and asymptotic alternatives use normal/chi-square references. All confidence intervals are pointwise. Panel A uses the same denominator conventions. Original domains reproduce the earlier period-specific age eligibility rules. Harmonized domains require BMI and known combined diabetes/cardiovascular-disease/cancer status in both periods. The common 80+ sensitivity caps recorded age at 80 before fitting. All scores retain their original within-cycle standardization. All tests are exploratory and unadjusted for multiplicity. BMI, body mass index; CI, confidence interval; df, degrees of freedom; SD, standard deviation.'

def probability(row):
    p=number(row,'p')
    if not 0<=p<=1:raise ValueError('Invalid p value')
    return f'{p:.5e}'


def build_tables(work):
    path=work/'age_comparison'
    marker=path/'COMPLETE.txt'
    if not marker.exists() or marker.read_text().strip()!='3 age-comparison models complete':
        raise ValueError('Age-comparison models are incomplete')
    joint=index(read(path/'joint_tests.csv'),['model','inference'])
    support=index(read(path/'support.csv'),['model','era'])
    contrasts=index(read(path/'contrasts.csv'),['model','inference','type'])
    pairs={(m,i) for m in MODELS for i in INFERENCES}
    if set(joint)!=pairs or set(support)!={(m,e) for m in MODELS for e in ['hip','wrist']}:
        raise ValueError('Incomplete or unexpected model coverage')
    if set(contrasts)!={(m,i,t) for m,i in pairs for t in ['hip','wrist','wrist_minus_hip']}:
        raise ValueError('Incomplete or unexpected contrasts')
    for row in support.values():
        for key in ['n','design_df']:
            x=number(row,key)
            if x<=0 or x!=int(x):raise ValueError('Invalid model support')
    for (model,inf),row in joint.items():
        df=float(row['denominator_df'])
        if inf=='asymptotic':
            if df!=math.inf:raise ValueError('Invalid asymptotic df')
        elif not math.isfinite(df) or df<=0 or df!=int(df):raise ValueError('Invalid finite df')
        if inf=='full_design' and df!=sum(number(support[(model,e)],'design_df') for e in ['hip','wrist']):
            raise ValueError('Full design df mismatch')
        if number(row,'numerator_df')!=2 or number(row,'F')<0:raise ValueError('Invalid joint test')
        if number(row,'F')!=number(joint[(model,'full_design')],'F'):raise ValueError('Joint statistic mismatch')
        probability(row)
        for kind in ['hip','wrist','wrist_minus_hip']:
            c=contrasts[(model,inf,kind)]
            if number(c,'age')!=75 or number(c,'reference')!=50 or float(c['denominator_df'])!=df:
                raise ValueError('Unexpected contrast or inference')
            if not number(c,'ci_low')<=number(c,'estimate')<=number(c,'ci_high') or number(c,'se')<=0:
                raise ValueError('Invalid interval or standard error')
            probability(c)
        delta=number(contrasts[(model,inf,'wrist')],'estimate')-number(contrasts[(model,inf,'hip')],'estimate')
        if abs(delta-number(contrasts[(model,inf,'wrist_minus_hip')],'estimate'))>1e-10:
            raise ValueError('Direct contrast does not equal wrist minus hip')
    a=[['Analysis','Hip n / Wrist n','F (2 df)','Full-design p','Residual p','Asymptotic p']]
    b=[['Analysis / inference','Hip difference, SD (95% CI)','Wrist difference, SD (95% CI)','Wrist minus hip, SD (95% CI)','Direct p']]
    for model,label in MODELS.items():
        counts=[int(number(support[(model,e)],'n')) for e in ['hip','wrist']]
        a.append([label,f'{counts[0]} / {counts[1]}',f"{number(joint[(model,'full_design')],'F'):.3f}"]+
                 [probability(joint[(model,i)]) for i in INFERENCES])
        for inf,inf_label in INFERENCES.items():
            rows=[contrasts[(model,inf,t)] for t in ['hip','wrist','wrist_minus_hip']]
            values=[f"{number(c,'estimate'):.3f} ({number(c,'ci_low'):.3f} to {number(c,'ci_high'):.3f})" for c in rows]
            b.append([label+' / '+inf_label]+values+[probability(rows[2])])
    return {'A':a,'B':b}


def export(work):
    tables=build_tables(work);out=work/'displays'
    names=['table_s40_a.csv','table_s40_b.csv','age_comparison_tables.json']
    if any((out/n).exists() for n in names):raise FileExistsError('Refusing to overwrite S40 exports')
    out.mkdir(exist_ok=True)
    for panel,name in zip(['A','B'],names):
        with (out/name).open('x',newline='') as f:csv.writer(f).writerows(tables[panel])
    with (out/names[2]).open('x') as f:
        json.dump({'Table S40':{'title':TITLE,'note':NOTE,'panels':tables}},f,ensure_ascii=False,indent=2);f.write('\n')
    print('Exported Table S40 (two panels; three models; nine contrast rows)')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    export(parser.parse_args().work_dir)
