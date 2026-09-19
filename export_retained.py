#!/usr/bin/env python3
"""Export retained supplement sensitivities from freshly fitted model results.

No historic coefficients are inputs. Repeated tables share identical fresh fits.
"""
import argparse
from pathlib import Path
import pandas as pd

HIP = {
    'exclude_deaths_first_2y': 'exclude_deaths_2y',
    'exclude_baseline_cvd_or_cancer': 'exclude_CVD_cancer',
    'pax_valid_days_3': 'valid_days_3',
    'pax_valid_days_5': 'valid_days_5',
    'pax_calibration_1_or_2': 'calibration_1_or_2',
    'pax_nonwear_90': 'nonwear_90',
    'pax_response_weighted': 'response_weighted',
}
WRIST = {
    'Pooled 2011-2014': 'primary',
    'Exclude deaths in first 2 years': 'exclude_deaths_2y',
    'Exclude baseline CVD or cancer': 'exclude_CVD_cancer',
    'MIMS per valid monitor minute': 'intensity',
    'At least 1 valid monitor day': 'one_day',
}
CONTRASTS = {
    'male': 'Male', 'female': 'Female',
    'female_vs_male_slope': 'Female/Male slope ratio',
}


def export(work):
    out = work / 'models'
    targets = [out/'legacy_sensitivities_unified.csv',
               out/'wrist_mortality_sensitivity_results.csv']
    if any(p.exists() for p in targets):
        raise FileExistsError('Refusing to overwrite retained exports')
    results = pd.concat([pd.read_csv(out/f'{name}.csv',float_precision='round_trip')
        for name in ('mortality_results','additional_sensitivity_results','retained_sensitivity_results')])
    if results.duplicated(['model','contrast','inference']).any():
        raise ValueError('Duplicate model/contrast/reference combination')
    results = results.set_index(['model','contrast','inference'])
    hip = []
    for label, model in HIP.items():
        for sex, contrast in CONTRASTS.items():
            fit_type = 'legacy_interaction' if model in ('valid_days_3','valid_days_5','calibration_1_or_2','nonwear_90') else 'interaction'
            rows = {conv: results.loc[(f'2005-2006__{fit_type}__{model}',contrast,conv)]
                    for conv in ('full_design_t','residual_t','normal')}
            full = rows['full_design_t']
            row = dict(model_id=label,sex=sex,coefficient=full.coefficient,
                       hazard_ratio=full.estimate,design_based_se=full.se,
                       n=full.n,deaths=full.deaths,design_df=full.design_df,
                       residual_df=full.residual_df)
            for conv, columns in {
                'full_design_t': ('ci_low_unified','ci_high_unified','p_unified'),
                'residual_t': ('ci_low_residual_df','ci_high_residual_df','p_value_residual_df'),
                'normal': ('native_ci_low','native_ci_high','native_p_value'),
            }.items():
                row.update(zip(columns,rows[conv][['ci_low','ci_high','p']]))
            hip.append(row)
    wrist = []
    labels = ['Relative contrast, Male','Relative contrast, Female','Female-to-male slope ratio']
    for label,model in WRIST.items():
        for contrast,display in zip(CONTRASTS.values(),labels):
            r = results.loc[(f'wrist__interaction__{model}',contrast,'full_design_t')]
            wrist.append(dict(analysis=label,contrast=display,coefficient=r.coefficient,
                hazard_ratio=r.estimate,design_based_se=r.se,statistic=r.coefficient/r.se,
                design_df=r.design_df,residual_df=r.residual_df,inference='full_design_t',
                denominator_df=r.denominator_df,ci_low=r.ci_low,ci_high=r.ci_high,
                p_value=r.p,n=r.n,deaths=r.deaths))
    for target,rows in zip(targets,(hip,wrist)):
        pd.DataFrame(rows).to_csv(target,index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,required=True)
    export(parser.parse_args().work_dir.resolve())
