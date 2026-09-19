"""Reject incomplete or invalid sensitivity exports before writing tables."""
from pathlib import Path
import unittest
from unittest.mock import patch
import pandas as pd
from export_no_transport import build_tables


def fixture(kind):
    contrasts = [('age', 'Age 80 versus 50')] if kind == 'age' else [
        ('common', 'Common contrast'), ('interaction', 'Male'),
        ('interaction', 'Female'), ('interaction', 'Female/Male slope ratio')]
    return pd.DataFrame([
        dict(model=f'{group}__{model}__{definition}', contrast=contrast,
             inference='full_design_t', n=100, deaths=10, estimate=1.2,
             ci_low=1.1, ci_high=1.3, p=.01)
        for group in ['2003-2004', '2005-2006', 'hip']
        for model, contrast in contrasts
        for definition in ['original', 'no_transport']])


class ExportSensitivityTests(unittest.TestCase):
    def build(self, mortality=None, age=None):
        with patch('export_no_transport.pd.read_csv', side_effect=[
                fixture('mortality') if mortality is None else mortality,
                fixture('age') if age is None else age]):
            return build_tables(Path('synthetic-work'))

    def test_complete_pair_coverage(self):
        tables = self.build()
        self.assertEqual([len(tables[k]) for k in ['Table S30', 'Table S31']], [13, 4])

    def test_missing_original_or_pair_rejected(self):
        frame = fixture('mortality')
        for partial in [frame.iloc[1:], frame.iloc[2:], frame.iloc[:0]]:
            with self.subTest(rows=len(partial)), self.assertRaises(ValueError):
                self.build(mortality=partial)
        with self.assertRaises(ValueError):
            self.build(age=fixture('age').iloc[2:])

    def test_duplicate_and_unexpected_keys_rejected(self):
        frame = fixture('mortality')
        for extra in [frame.iloc[[0]], frame.iloc[[0]].assign(model='unexpected__original')]:
            with self.subTest(extra=extra.model.iloc[0]), self.assertRaises(ValueError):
                self.build(mortality=pd.concat([frame, extra], ignore_index=True))

    def test_nonfinite_and_invalid_results_rejected(self):
        for column, value in [('estimate', float('inf')), ('p', float('nan')),
                              ('p', -0.1), ('p', 1.1), ('ci_low', 1.5),
                              ('n', 100.5), ('deaths', 101)]:
            frame = fixture('mortality')
            if column in ['n', 'deaths']:frame[column] = frame[column].astype(float)
            frame.loc[0, column] = value
            with self.subTest(column=column, value=value), self.assertRaises(ValueError):
                self.build(mortality=frame)

    def test_mismatched_pair_support_rejected(self):
        frame = fixture('mortality');frame.loc[1, 'n'] = 101
        with self.assertRaises(ValueError):
            self.build(mortality=frame)


if __name__ == '__main__':
    unittest.main()
