#!/usr/bin/env python3
"""Generate current main-text tables from new full-design model outputs."""
from pathlib import Path
import argparse,csv,json
import pandas as pd
def tables(work):
    here = work / 'models'
    m = pd.read_csv(here / 'mortality_results.csv')
    m = m[(m.inference == 'full_design_t') & (m.contrast == 'Common contrast')]
    def r(name):
        subset = m[m.model == name]
        assert len(subset) == 1, name
        return subset.iloc[0]
    def val(row): return f'{row.estimate:.3f} [{row.ci_low:.3f}-{row.ci_high:.3f}]'
    def p(row): return '< 0.001' if row.p < 0.001 else f'{row.p:.3f}'
    t2 = [['Analysis', 'n / deaths', 'Common HR', '95% CI', 'p', 'df']]
    for name, label in [('2003-2004','2003-2004'),('2005-2006','2005-2006'),('hip','Pooled 2003-2006'),('2011-2012','2011-2012'),('2013-2014','2013-2014'),('wrist','Pooled 2011-2014')]:
        row = r(f'{name}__common__primary')
        t2.append([label, f'{int(row.n):,} / {int(row.deaths):,}', f'{row.estimate:.3f}', f'{row.ci_low:.3f}-{row.ci_high:.3f}', p(row), str(int(row.design_df))])
    t3 = [['Sensitivity analysis', 'Period', 'n / deaths', 'Common HR [95% CI]', 'p']]
    for code, label in [('exclude_deaths_3y','Exclude deaths within 3 years'),('exclude_deaths_5y','Exclude deaths within 5 years'),('censor_5y','Censor follow-up at 5 years'),('equal','Vigorous activity weighted once'),('zlog','Conventional z scores of log values'),('zraw','Conventional z scores of raw values')]:
        for name, period in [('hip','2003-2006'),('wrist','2011-2014')]:
            row = r(f'{name}__common__{code}')
            t3.append([label, period, f'{int(row.n):,} / {int(row.deaths):,}', val(row), p(row)])
    age=pd.read_csv(here / 'age_results.csv')
    age=age[(age.inference=='full_design_t') & age.model.str.endswith('__primary')]
    tests=pd.read_csv(here / 'joint_tests.csv')
    t4=[['Analysis','Age 80 versus 50, SD [95% CI]','p','Global age p','df']]
    for name,label in [('2003-2004','2003-2004'),('2005-2006','2005-2006'),('hip','Pooled 2003-2006'),('2011-2012','2011-2012'),('2013-2014','2013-2014'),('wrist','Pooled 2011-2014')]:
        row=age[age.model==f'{name}__age__primary'].iloc[0]
        test=tests[(tests.model==row.model)&(tests.test=='Global age')].iloc[0]
        t4.append([label,val(row),p(row),'< 0.001' if test.p_full<.001 else f'{test.p_full:.3f}',str(int(row.design_df))])
    reg = pd.read_csv(here/'regime_nuisance_results.csv',float_precision='round_trip')
    def regime(horizon):
        row = reg[(reg.model==f'all_nuisance_era__{horizon}') & (reg.inference=='full_design_t')]
        assert len(row)==1
        return row.iloc[0]
    row=regime('Inf')
    t2=[['Panel A. Period- and cycle-specific associations']]+t2+[
        ['Panel B. Direct measurement-regime comparison'],
        ['Comparison','n / deaths','HR ratio','95% CI','p','df'],
        ['Wrist versus hip',f'{int(row.n):,} / {int(row.deaths):,}',f'{row.estimate:.3f}',
         f'{row.ci_low:.3f}-{row.ci_high:.3f}',p(row),str(int(row.design_df))]]
    row=regime('5')
    t3=[['Panel A. Period-specific sensitivity analyses']]+t3+[
        ['Panel B. Direct measurement-regime sensitivity comparison'],
        ['Sensitivity analysis','Comparison','n / deaths','HR ratio [95% CI]','p'],
        ['Follow-up censored at 5 years','Wrist versus hip',f'{int(row.n):,} / {int(row.deaths):,}',val(row),p(row)]]
    return {'Table 2':t2,'Table 3':t3,'Table 4':t4}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    work=parser.parse_args().work_dir.resolve()
    out=work/'displays'
    out.mkdir(exist_ok=True)
    target=out/'main_tables.json'
    if target.exists():raise FileExistsError('Refusing to overwrite main tables')
    result=tables(work)
    with (out/'table_1_two_panel_characteristics.csv').open() as handle:
        result['Table 1']=list(csv.reader(handle))[1:]
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    for name,rows in result.items():
        with (out/(name.lower().replace(' ','_')+'.csv')).open('x',newline='') as handle:
            csv.writer(handle).writerows(rows)
