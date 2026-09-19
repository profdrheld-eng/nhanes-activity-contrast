#!/usr/bin/env python3
"""Export same-sample, post-hoc hip-score sensitivity tables from fresh fits."""
import argparse
import csv
import json
import math
from pathlib import Path
import pandas as pd


def pvalue(value):
    return '< 0.001' if value < .001 else f'{value:.3f}'


def interval(row):
    return f'{row.estimate:.3f} [{row.ci_low:.3f}, {row.ci_high:.3f}]'


def build_tables(work):
    tables = {}
    for name, label in [('mortality', 'Table S30'), ('age', 'Table S31')]:
        frame = pd.read_csv(work/'models'/f'no_transport_{name}.csv')
        frame = frame.loc[frame.inference.eq('full_design_t')]
        contrasts = [('age', 'Age 80 versus 50')] if name == 'age' else [
            ('common', 'Common contrast'), ('interaction', 'Male'),
            ('interaction', 'Female'), ('interaction', 'Female/Male slope ratio')]
        expected_keys = {
            (f'{group}__{model}__{definition}', contrast)
            for group in ['2003-2004', '2005-2006', 'hip']
            for model, contrast in contrasts
            for definition in ['original', 'no_transport']}
        keys = list(zip(frame.model, frame.contrast))
        if len(keys) != len(set(keys)) or set(keys) != expected_keys:
            raise ValueError('Sensitivity export requires complete, unique expected model contrasts')
        for row in frame.itertuples():
            if not all(math.isfinite(v) for v in [row.estimate, row.ci_low, row.ci_high,
                                                  row.p, row.n, row.deaths]):
                raise ValueError('Sensitivity export contains non-finite results')
            if not (0 <= row.p <= 1 and row.ci_low <= row.estimate <= row.ci_high):
                raise ValueError('Sensitivity export contains invalid intervals or probabilities')
            if not (row.n > 0 and 0 <= row.deaths <= row.n
                    and row.n == int(row.n) and row.deaths == int(row.deaths)):
                raise ValueError('Sensitivity export contains invalid sample counts')
        rows = [['Analysis', 'Contrast', 'n / deaths' if name == 'mortality' else 'n',
                 'Original estimate [95% CI]', 'p', 'Without transport [95% CI]', 'p']]
        for original in frame.loc[frame.model.str.endswith('__original')].itertuples():
            expected = original.model.removesuffix('__original')+'__no_transport'
            matched = frame.loc[frame.model.eq(expected) & frame.contrast.eq(original.contrast)]
            if len(matched) != 1:
                raise ValueError('Each original estimate requires exactly one sensitivity estimate')
            other = next(matched.itertuples())
            if (original.n, original.deaths) != (other.n, other.deaths):
                raise ValueError('Sensitivity analysis changed the analysis sample or events')
            group = original.model.split('__')[0]
            group = '2003–2006 pooled' if group == 'hip' else group
            size = f'{original.n:,} / {original.deaths:,}' if name == 'mortality' else f'{original.n:,}'
            rows.append([group, original.contrast, size, interval(original), pvalue(original.p),
                         interval(other), pvalue(other.p)])
        tables[label] = rows
    return tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    out = args.work_dir/'displays'
    targets = [out/'no_transport_tables.json', out/'table_s30.csv', out/'table_s31.csv']
    if any(path.exists() for path in targets):
        raise FileExistsError('Refusing to overwrite completed or partial no-transport table outputs')
    tables = build_tables(args.work_dir)
    out.mkdir(parents=True, exist_ok=True)
    targets[0].write_text(json.dumps(tables, indent=2, ensure_ascii=False)+'\n')
    for label, rows in tables.items():
        with (out/(label.lower().replace(' ', '_')+'.csv')).open('x', newline='') as handle:
            csv.writer(handle).writerows(rows)


if __name__ == '__main__':
    main()
