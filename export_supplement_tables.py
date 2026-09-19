#!/usr/bin/env python3
"""Generate numerical supplement tables from newly fitted model outputs.

Method-specification tables S1/S4/S5 are handled separately from numeric results.
Presentation labels require reconciliation to the current working document.
"""
from pathlib import Path
import argparse,csv,json,re
import pandas as pd
import numpy as np
CY=['2003-2004','2005-2006','2011-2012','2013-2014']
def read(name): return pd.read_csv(OUT/(name+'.csv'))

def p(v): return 'not estimable' if pd.isna(v) else '< 0.001' if v < .001 else f'{v:.3f}'

def n(v): return '' if pd.isna(v) else f'{int(v):,}'

def f(v): return '' if pd.isna(v) else f'{v:.3f}'

def presentation_label(value):
    """Expand legacy display labels without changing source values or terms."""
    value = {'female_vs_male_slope': 'Female/Male slope ratio',
             'Regime heterogeneity': 'Regime difference',
             'Cycle heterogeneity': 'Cycle difference',
             'quintile joint': 'All quintiles'}.get(str(value), str(value))
    value = re.sub(r'\b(?:women|female)\b', 'Female', str(value), flags=re.I)
    value = re.sub(r'\b(?:men|male)\b', 'Male', value, flags=re.I)
    value = value.replace('High school or GED', 'High school or equivalent')
    value = value.replace('bmi extension', 'BMI extension')
    return value

def variable(v, cycle=None):
    if v == "SR_day" and cycle in ("2003-2004", "2005-2006"):
        return "Self-report index, units/day"
    return {'age':'Age, years','BMI':'Body mass index, kg/m²','PIR':'Poverty-income ratio',
        'SR_day':'Self-report, minutes/day','sex=female':'Female sex',
        'smoking3=current':'Current smoking','CVD_history=yes':'Cardiovascular disease',
        'cancer_history=yes':'Cancer','diabetes_history=yes':'Diabetes'}.get(v,str(v).replace('_',' '))

def term(v):
    labels={'Delta':'Activity contrast','Level':'Common activity component',
        'PIR':'Poverty-income ratio','age':'Age','sex':'Sex','sexfemale':'Female versus Male',
        'Delta:sex':'Contrast by sex','Delta:sexfemale':'Contrast by Female sex',
        'GLOBAL':'All model terms','education3':'Education','race4':'Race and ethnicity',
        'smoking3':'Smoking status','race4Hispanic':'Hispanic versus non-Hispanic White',
        'race4Non-Hispanic Black':'Non-Hispanic Black versus non-Hispanic White',
        'race4Other':'Other race/ethnicity versus non-Hispanic White',
        'education3high school/GED':'High school/equivalent versus more than high school',
        'education3under high school':'Less than high school versus more than high school',
        'smoking3former':'Former versus never smoking','smoking3current':'Current versus never smoking',
        'late_delta':'Contrast change after 5 years','late_sex':'Sex coefficient change after 5 years',
        'late_interaction':'Contrast-by-sex change after 5 years',
        'tt(Delta)':'Contrast by log follow-up time','tt(sex_female)':'Female sex by log follow-up time'}
    assert v in labels, v
    return labels[v]

def model(v):
    swaps = {'hip':'2003-2006','wrist':'2011-2014','common':'common slope','interaction':'sex interaction',
        'primary':'primary definition','equal':'vigorous weighted once','zlog':'z scores of log values',
        'zraw':'z scores of raw values','exclude_deaths_2y':'exclude deaths within 2 years',
        'exclude_deaths_3y':'exclude deaths within 3 years','exclude_deaths_5y':'exclude deaths within 5 years',
        'censor_5y':'censor at 5 years','censor_8y':'censor at 8 years','landmark_3y':'3-year landmark',
        'landmark_5y':'5-year landmark','matched_core':'matched core model','matched_BMI':'matched plus BMI',
        'matched_BMI_disease':'matched plus BMI and disease','exclude_CVD_cancer':'exclude CVD or cancer',
        'quintiles':'quintiles','spline':'restricted cubic spline','age':'age model',
        'all_nuisance_era':'all covariate slopes vary by regime','sex_era':'sex slope varies by regime',
        'Inf':'full follow-up','5':'5-year censoring','8':'8-year censoring',
        'exclude_deaths_first_2y':'exclude deaths within 2 years',
        'exclude_baseline_cvd_or_cancer':'exclude baseline cardiovascular disease or cancer',
        'pax_valid_days_3':'at least 3 valid wear days',
        'pax_valid_days_5':'at least 5 valid wear days',
        'pax_calibration_1_or_2':'device calibration codes 1 or 2',
        'pax_nonwear_90':'90-minute nonwear rule',
        'pax_response_weighted':'response-weighted model'}
    return '; '.join(swaps.get(x,x.replace('_',' ')) for x in str(v).split('__'))

def reference_table(frame):
    """One row per coefficient, with all three references side by side."""
    headers=['Analysis','Contrast','n / deaths','Estimate','Full-design 95% CI; p','Residual 95% CI; p','Normal 95% CI; p']
    rows=[headers]
    for (name,contrast), group in frame.groupby(['model','contrast'],sort=False):
        refs={r.inference:r for r in group.itertuples()}
        assert set(refs)=={'full_design_t','residual_t','normal'},(name,contrast)
        main=refs['full_design_t']
        for r in refs.values():
            assert np.isclose(r.estimate,main.estimate) and np.isclose(r.se,main.se)
        intervals=[f'{refs[k].ci_low:.3f} to {refs[k].ci_high:.3f}; {p(refs[k].p)}' for k in ['full_design_t','residual_t','normal']]
        rows.append([model(name),contrast,f'{n(main.n)} / {n(main.deaths)}',f(main.estimate),*intervals])
    return rows

def old_tables():
    im=pd.read_csv(LEGACY/'mortality_modifiers.csv')
    im=im[im.convention=='full_design_t']
    s1=[['Modifier','Scale','n / deaths','Interaction HR ratio [95% CI]','p','BH q','Design df']]
    for r in im.itertuples():
        s1.append([r.modifier,r.scale,f'{n(r.n)} / {n(r.deaths)}',f'{r.estimate:.3f} [{r.ci_low:.3f}-{r.ci_high:.3f}]',p(r.p),p(r.q_bh),n(r.design_df)])
    om=pd.read_csv(LEGACY/'participant_omnibus.csv')
    om=om[om.convention=='full_design_F']
    s2=[['Characteristic','Model','n','F (df1, df2)','p','BH q']]
    for r in om.itertuples():
        characteristic = {'PIR':'Poverty-income ratio', 'BMI':'Body mass index'}.get(r.characteristic, r.characteristic.replace('_',' '))
        s2.append([characteristic,r.model.replace('_',' '),n(r.n),f'{r.F:.3f} ({int(r.numerator_df)}, {int(r.denominator_df)})',p(r.p),p(r.q_bh)])
    co=pd.read_csv(LEGACY/'participant_correlates.csv')
    co=co[co.convention=='full_design_t']
    # Original Table S3 contains core contrasts plus the designated extension only.
    selected=co[((co.model=='core') & ~co.term.isin(['Level','age10_c50','I(age10_c50^2)'])) |
        ((co.model=='bmi_extension') & co.term.str.contains('BMI')) |
        ((co.model=='smoking_extension') & co.term.str.contains('smoking')) |
        ((co.model=='disease_extension') & co.term.str.contains('disease'))]
    s3=[['Contrast','n','Difference, SD','95% CI','p','Model','Design df']]
    for r in selected.itertuples():
        s3.append([r.label.replace(' versus male',' versus Male'),n(r.n),f(r.estimate),f'{r.ci_low:.3f} to {r.ci_high:.3f}',p(r.p),r.model.replace('_',' '),n(r.design_df)])
    assert len(s3)==12,selected[['model','term']]
    m=read('mortality_results')
    primary=m[m.model.str.endswith('__primary')]
    s4=[['Analysis','Contrast','Reference','df','HR [95% CI]','p']]
    labels={'full_design_t':'Full-design t','residual_t':'Residual t','normal':'Normal'}
    for r in primary.itertuples():
        s4.append([model(r.model),r.contrast,labels[r.inference], 'Infinity' if np.isinf(r.denominator_df) else n(r.denominator_df),f'{r.estimate:.3f} [{r.ci_low:.3f}-{r.ci_high:.3f}]',p(r.p)])
    return {188:s1,192:s2,197:s3,205:s4}

def additions():
    blocks = {}

    def table(number, rows):
        mapping = {6: 2, 8: 6, '8a': 7, '8b': 8, 16: 3}
        number = mapping.get(number, number)
        blocks[f'Table S{number}'] = rows
    flow = read('participant_flow')
    rows = [['Criterion', *CY]]
    for stage, g in flow.groupby('stage', sort=True):
        cells = []
        for cy in CY:
            r = g[g.cycle == cy].iloc[0]
            cells.append(f'{n(r.remaining)}' + (f' (excluded {n(r.excluded)})' if r.excluded > 0 else ''))
        rows.append([g.iloc[0].criterion, *cells])
    table(6, rows)
    dist = read('contrast_distribution')
    rows = [['Cycle and weighting', 'n', 'Mean ± SD', 'Median', 'Range', '5th-95th percentile']]
    for r in dist.itertuples():
        rows.append([r.cycle + '; ' + r.weighting, n(r.n), f'{r.mean:.3f} ± {r.sd:.3f}', f(r.median), f'{r.min:.3f} to {r.max:.3f}', f'{r.p05:.3f} to {r.p95:.3f}'])
    table(8, rows)
    ties = read('rank_ties_and_bounds')
    zero = read('contrast_zero_reference')
    rows = [['Cycle', 'Measure', 'n', 'Unique values', 'Participants in tied groups', 'Clipped low / high']]
    for r in ties.itertuples():
        rows.append([r.cycle, 'Self-report' if r.measure == 'SR_day' else 'Device', n(r.n), n(r.unique_values), n(r.persons_in_tied_groups), f'{n(r.clipped_low)} / {n(r.clipped_high)}'])
    table('8a', rows)
    rows = [['Cycle', 'Mean raw difference', 'Raw difference SD', 'Centered contrast at equal ranks']]
    for r in zero.itertuples():
        rows.append([r.cycle, f(r.mean), f(r.sd), f(r.equal_rank_location_standardized)])
    table('8b', rows)
    co = read('correlations')
    rows = [['Cycle', 'Sample', 'Pair', 'n', 'Weighted r', 'Unweighted r']]
    for r in co.itertuples():
        rows.append([r.cycle, r.sample, r.pair, n(r.n), f(r.weighted_r), f(r.unweighted_r)])
    table(9, rows)
    sup = read('model_support')
    prim = sup[sup.model.str.contains('__(?:common|interaction)__primary$', regex=True)]
    rows = [['Model', 'n / deaths', 'Strata', 'PSUs', 'Design df', 'Parameters', 'Residual df']]
    for r in prim.itertuples():
        rows.append([model(r.model), f'{n(r.n)} / {n(r.deaths)}', n(r.strata), n(r.PSUs), n(r.design_df), n(r.parameters), n(r.residual_df)])
    assert len(prim) == 12
    table(10, rows)
    m = read('mortality_results')
    table(11, reference_table(m[~m.model.str.endswith('__primary')]))
    table(12, reference_table(read('additional_sensitivity_results')))
    reg = read('regime_nuisance_results')
    table(13, reference_table(reg))
    jt = pd.concat([read('joint_tests'), read('regime_nuisance_tests')], ignore_index=True)
    rows = [['Analysis', 'Test', 'F (numerator df)', 'Design / residual df', 'Full p', 'Residual p', 'Normal p']]
    for r in jt.itertuples():
        rows.append([model(r.model), r.test, f'{r.F:.3f} ({int(r.numerator_df)})', f'{n(r.design_df)} / {n(r.residual_df)}', p(r.p_full), p(r.p_residual), p(r.p_normal)])
    table(14, rows)
    legacy = read('legacy_sensitivities_unified')
    rows = [['2005-2006 analysis', 'Contrast', 'n / deaths', 'HR', 'Full-design CI; p', 'Residual CI; p', 'Normal CI; p']]
    for r in legacy.itertuples():
        rows.append([model(r.model_id), r.sex, f'{n(r.n)} / {n(r.deaths)}', f(r.hazard_ratio), f'{r.ci_low_unified:.3f} to {r.ci_high_unified:.3f}; {p(r.p_unified)}', f'{r.ci_low_residual_df:.3f} to {r.ci_high_residual_df:.3f}; {p(r.p_value_residual_df)}', f'{r.native_ci_low:.3f} to {r.native_ci_high:.3f}; {p(r.native_p_value)}'])
    table(15, rows)
    inc = read('inclusion_comparison')
    rows = [['Cycle', 'Group', 'Characteristic', 'Group n', 'Observed n', 'Weighted mean ± SD or proportion']]
    for r in inc.itertuples():
        rows.append([r.cycle, r.group, variable(r.variable, r.cycle), n(r.n_group), n(r.n_observed), f'{r.mean_or_proportion:.3f} ± {r.sd:.3f}' if pd.notna(r.sd) else f(r.mean_or_proportion)])
    table(16, rows)
    ph = read('ph_diagnostics')
    rows = [['Model', 'Term', 'Chi-square', 'df', 'p']]
    for r in ph.itertuples():
        rows.append([model(r.model), term(r.term), f(r.chisq), n(r.df), p(r.p)])
    table(17, rows)
    tv = read('time_varying')
    rows = [['Model and method', 'Term', 'n / deaths', 'Ratio [95% CI]', 'p', 'Reference df']]
    for r in tv.itertuples():
        rows.append([model(r.model) + '; ' + r.method, term(r.term), f'{n(r.n)} / {n(r.deaths)}', f'{r.ratio:.3f} [{r.ci_low:.3f}-{r.ci_high:.3f}]', p(r.p), 'Normal' if pd.isna(r.denominator_df) or np.isinf(r.denominator_df) else n(r.denominator_df)])
    table(18, rows)
    ds = read('diagnostic_support')
    rows = [['Model', 'n / deaths', 'Standardized matrix condition number', 'Legacy Gram condition number', 'Maximum absolute weighted dfbeta']]
    for r in ds.itertuples():
        rows.append([model(r.model), f'{n(r.n)} / {n(r.deaths)}', f(r.weighted_matrix_condition_number), f(r.legacy_gram_condition_number), f(r.max_abs_weighted_dfbeta)])
    table(19, rows)
    influence = read('influence')
    rows = [['Model', 'Coefficient', 'Maximum absolute weighted dfbeta']]
    for r in influence.itertuples():
        rows.append([model(r.model), term(r.term), f'{r.max_abs_weighted_dfbeta:.6f}'])
    table('19a', rows)
    sp = read('spline_support')
    rows = [['Period', 'Knots (10%, 50%, 90%)', 'Display range (1%-99%)', 'Tail bounds (5%-95%)', 'Low-tail n / deaths', 'High-tail n / deaths']]
    for r in sp.itertuples():
        rows.append([model(r.era), f'{r.k10:.3f}, {r.k50:.3f}, {r.k90:.3f}', f'{r.p01:.3f} to {r.p99:.3f}', f'{r.p05:.3f} to {r.p95:.3f}', f'{n(r.low_tail_n)} / {n(r.low_tail_deaths)}', f'{n(r.high_tail_n)} / {n(r.high_tail_deaths)}'])
    table(20, rows)
    table(21, reference_table(read('age_results')))
    table(22, reference_table(read('age_support_sensitivity')))
    overlay = pd.read_csv(OUT / 'age_sex_overlay.csv')
    age_n = {era: int(g.n.iloc[0]) for era, g in overlay.groupby('era')}
    assert all((g.n.nunique() == 1 for _, g in overlay.groupby('era')))
    rows = [list(overlay.columns)] + overlay.fillna('').astype(str).values.tolist()
    table(23, [['Period', 'Participants', 'Common model', 'Exploratory sex-specific model'], ['2003-2006', n(age_n['hip']), 'Linear and quadratic age plus core age covariates and cycle', 'Same domain and covariates, adding both age-by-sex interactions'], ['2011-2014', n(age_n['wrist']), 'Linear and quadratic age plus core age covariates and cycle', 'Same domain and covariates, adding both age-by-sex interactions']])
    for suffix, filename in [('a', 'mortality_modifiers'), ('b', 'participant_correlates'), ('c', 'participant_omnibus')]:
        x = pd.read_csv(LEGACY / (filename + '.csv'))
        if suffix == 'a':
            cols = ['modifier', 'scale', 'convention', 'estimate', 'ci_low', 'ci_high', 'p', 'q_bh']
        elif suffix == 'b':
            cols = ['model', 'label', 'convention', 'estimate', 'ci_low', 'ci_high', 'p']
        else:
            cols = ['model', 'characteristic', 'convention', 'F', 'numerator_df', 'denominator_df', 'p', 'q_bh']
        headers = {'modifier': 'Modifier', 'scale': 'Scale', 'convention': 'Reference', 'estimate': 'Estimate', 'ci_low': 'CI lower', 'ci_high': 'CI upper', 'p': 'p', 'q_bh': 'BH q', 'model': 'Model', 'label': 'Comparison', 'characteristic': 'Characteristic', 'F': 'F', 'numerator_df': 'df1', 'denominator_df': 'df2'}
        rows = [[headers[c] for c in cols]]
        for _, r in x.iterrows():
            rows.append([p(r[c]) if c in ['p', 'q_bh'] else f(r[c]) if isinstance(r[c], (int, float, np.number)) else str(r[c]).replace('_', ' ') for c in cols])
        table('24' + suffix, rows)
    return blocks

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    work=parser.parse_args().work_dir.resolve()
    OUT=work/'models';LEGACY=work/'exploratory'
    output=work/'displays/supplement_tables.json'
    if output.exists():raise FileExistsError('Refusing to overwrite supplement tables')
    result=additions()
    method_file=Path(__file__).resolve().parent/'config/method_tables.json'
    result.update(json.loads(method_file.read_text()))
    result.update({f'Table S{number}':rows for number,rows in zip(range(25,29),old_tables().values())})
    wrist=read('wrist_mortality_sensitivity_results')
    result['Table S29']=[['Analysis','Contrast','n / deaths','HR [95% CI]','p']]+[
        [r.analysis,r.contrast,f'{n(r.n)} / {n(r.deaths)}',f'{r.hazard_ratio:.3f} [{r.ci_low:.3f}-{r.ci_high:.3f}]',p(r.p_value)]
        for r in wrist.itertuples()]
    # Current presentation wording is distinct from the fitted variable names.
    headers={
        'Table S6':['Survey cycle and weighting','n','Mean ± SD','Median','Range','5th-95th percentile'],
        'Table S7':['Survey cycle','Activity Measure','n','Unique observed values','Participants in tied groups','Bounded in lower tail /bounded in upper tail'],
        'Table S8':['Survey cycle','Weighted mean of D','Weighted SD of D','Value of Δ at equal ranks'],
        'Table S9':['Survey cycle','Analytical sample','Variable pair','n','Survey weighted r','Unweighted r'],
        'Table S29':['Analysis','Contrast','n / deaths','HR [95% CI]','p value']}
    for name,header in headers.items():result[name][0]=header
    pairs={'raw_SR_device':'Raw self-report and device measures',
           'log_SR_device':'Log-transformed self-report and device measures',
           'normal_rank_SR_device':'Inverse-normal self-report and device scores',
           'average_difference':'Common activity component and signed contrast'}
    for row in result['Table S9'][1:]:row[2]=pairs[row[2]]
    analyses={'Pooled 2011-2014':'Primary','Exclude deaths in first 2 years':'Exclude early deaths',
              'Exclude baseline CVD or cancer':'Exclude CVD/cancer','MIMS per valid monitor minute':'MIMS/valid minute',
              'At least 1 valid monitor day':'>=1 valid day'}
    contrasts={'Relative contrast, Male':'Male','Relative contrast, Female':'Female',
               'Female-to-male slope ratio':'Female-to-Male ratio'}
    for row in result['Table S29'][1:]:row[0]=analyses[row[0]];row[1]=contrasts[row[1]]
    result={name:[[presentation_label(cell) for cell in row] for row in rows] for name,rows in result.items()}
    output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    for name,rows in result.items():
        with (output.parent/(name.lower().replace(' ','_')+'.csv')).open('x',newline='') as handle:
            csv.writer(handle).writerows(rows)
