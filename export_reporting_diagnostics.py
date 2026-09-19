#!/usr/bin/env python3
"""Export Table S37 from completed reporting diagnostics, without refitting models."""
import argparse
import csv
import json
import math
from pathlib import Path

GROUPS=['2003-2004','2005-2006','2011-2012','2013-2014','hip','wrist']
ADJUSTMENTS={'unadjusted':'Unadjusted','level_only':'Average component only','principal':'Principal adjustment'}
MODELS=[group+'__'+adjustment for group in GROUPS for adjustment in ADJUSTMENTS]
TITLE='Time-dependence diagnostics for the adjustment comparisons'
NOTE=('Contrast and global p values are conventional Schoenfeld-residual diagnostics, not survey-design-based tests. '
      'The late/early hazard ratio (HR) ratio compares the contrast association after five years with that during the first five years. '
      'Its 95% confidence interval (CI) and step p value use full-survey design-based t inference. '
      'All 18 models correspond to the same-case adjustment comparisons in Table S36. '
      'These post-hoc diagnostics do not replace the original models; a nonsignificant test does not establish proportionality.')


def read(path):
    with path.open(newline='') as handle:return list(csv.DictReader(handle))


def index(rows, keys):
    result={}
    for row in rows:
        key=tuple(row[k] for k in keys)
        if key in result:raise ValueError('Duplicate diagnostic key: '+str(key))
        result[key]=row
    return result


def number(row, name):
    value=float(row[name])
    if not math.isfinite(value):raise ValueError('Nonfinite '+name)
    return value


def probability(row):
    p=number(row,'p')
    if not 0<=p<=1:raise ValueError('Invalid p value')
    return '<0.001' if p<.001 else f'{p:.3f}'


def build_table(work):
    path=work/'reporting_diagnostics'
    marker=path/'COMPLETE.txt'
    if not marker.exists() or marker.read_text().strip()!='18 model diagnostics complete':
        raise ValueError('Reporting diagnostics are incomplete')
    ph=index(read(path/'schoenfeld.csv'),['model','term'])
    steps=index(read(path/'step_time.csv'),['model','contrast'])
    support=index(read(path/'support.csv'),['model'])
    expected={(m,c) for m in MODELS for c in ['early','late','late_to_early']}
    if set(steps)!=expected or set(support)!={(m,) for m in MODELS}:
        raise ValueError('Incomplete or unexpected model coverage')
    if {m for m,t in ph}!=set(MODELS) or any((m,t) not in ph for m in MODELS for t in ['Delta','GLOBAL']):
        raise ValueError('Incomplete Schoenfeld diagnostics')
    for (model,term),row in ph.items():
        if model!=row['group']+'__'+row['adjustment']:raise ValueError('Inconsistent model label')
        probability(row)
    for (model,contrast),row in steps.items():
        if model!=row['group']+'__'+row['adjustment'] or row['inference']!='full_design_t' or number(row,'cut_years')!=5:
            raise ValueError('Unexpected model, inference or time cut')
        probability(row)
        if not 0<number(row,'ci_low')<=number(row,'estimate')<=number(row,'ci_high'):
            raise ValueError('Invalid HR or confidence interval')
        ref=support[(model,)]
        for key in ['n','deaths','design_df']:
            value=number(row,key)
            if value!=number(ref,key) or value!=int(value):raise ValueError('Inconsistent model support')
        if not 0<number(row,'deaths')<=number(row,'n') or number(row,'design_df')<=0:
            raise ValueError('Invalid model support')
    table=[['Sample / adjustment','Contrast p','Global p','Late/early HR ratio [95% CI]','Step p']]
    for model in MODELS:
        row=steps[(model,'late_to_early')];group=row['group']
        sample={'hip':'2003–2006 pooled','wrist':'2011–2014 pooled'}.get(group,group)
        table.append([sample+' / '+ADJUSTMENTS[row['adjustment']],probability(ph[(model,'Delta')]),
            probability(ph[(model,'GLOBAL')]),f"{float(row['estimate']):.3f} [{float(row['ci_low']):.3f}, {float(row['ci_high']):.3f}]",probability(row)])
    return table


def export(work):
    table=build_table(work);out=work/'displays'
    names=['table_s37.csv','reporting_diagnostic_tables.json','reporting_diagnostic_table_notes.json']
    if any((out/name).exists() for name in names):raise FileExistsError('Refusing to overwrite S37 exports')
    out.mkdir(exist_ok=True)
    with (out/names[0]).open('x',newline='') as handle:csv.writer(handle).writerows(table)
    with (out/names[1]).open('x') as handle:json.dump({'Table S37':table},handle,ensure_ascii=False,indent=2);handle.write('\n')
    with (out/names[2]).open('x') as handle:
        json.dump({'titles':{'Table S37':TITLE},'notes':{'Table S37':NOTE}},handle,ensure_ascii=False,indent=2);handle.write('\n')
    print('Exported Table S37 (18 models)')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--work-dir',type=Path,required=True)
    export(parser.parse_args().work_dir)
