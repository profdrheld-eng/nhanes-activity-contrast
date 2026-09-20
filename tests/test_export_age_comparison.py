"""Synthetic coverage, consistency and output-guard tests for S40."""
import csv
import tempfile
import unittest
from pathlib import Path
from export_age_comparison import build_tables, export, MODELS, INFERENCES

class AgeExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.work=Path(self.tmp.name)
        self.out=self.work/'age_comparison';self.out.mkdir()
        self.joint=[];self.support=[];self.contrasts=[]
        for model in MODELS:
            for era in ['hip','wrist']:
                self.support.append(dict(model=model,era=era,n=100,design_df=31))
            for inf in INFERENCES:
                df={'full_design':62,'residual':39,'asymptotic':'Inf'}[inf]
                self.joint.append(dict(model=model,inference=inf,F=3,numerator_df=2,denominator_df=df,p=.05))
                for kind,estimate in [('hip',1),('wrist',.5),('wrist_minus_hip',-.5)]:
                    self.contrasts.append(dict(model=model,inference=inf,type=kind,age=75,reference=50,denominator_df=df,estimate=estimate,se=.1,ci_low=estimate-.2,ci_high=estimate+.2,p=.05))
        self.save()
    def save(self):
        for name,rows in [('joint_tests',self.joint),('support',self.support),('contrasts',self.contrasts)]:
            with (self.out/(name+'.csv')).open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
        (self.out/'COMPLETE.txt').write_text('3 age-comparison models complete\n')
    def test_all_panels_and_contrasts(self):
        tables=build_tables(self.work)
        self.assertEqual(len(tables['A']),4)
        self.assertEqual(len(tables['B']),10)
        self.assertEqual(tables['B'][1][3],'-0.500 (-0.700 to -0.300)')
    def test_missing_completion(self):
        (self.out/'COMPLETE.txt').unlink()
        with self.assertRaises(ValueError):build_tables(self.work)
    def test_duplicate_joint(self):
        self.joint.append(self.joint[0]);self.save()
        with self.assertRaises(ValueError):build_tables(self.work)
    def test_missing_contrast(self):
        self.contrasts.pop();self.save()
        with self.assertRaises(ValueError):build_tables(self.work)
    def test_wrong_age(self):
        self.contrasts[0]['age']=80;self.save()
        with self.assertRaises(ValueError):build_tables(self.work)
    def test_invalid_interval(self):
        self.contrasts[0]['ci_low']=2;self.save()
        with self.assertRaises(ValueError):build_tables(self.work)
    def test_nonfinite_probability(self):
        self.joint[0]['p']='NaN';self.save()
        with self.assertRaises(ValueError):build_tables(self.work)
    def test_difference_identity(self):
        self.contrasts[2]['estimate']=-.4;self.save()
        with self.assertRaises(ValueError):build_tables(self.work)
    def test_inference_df(self):
        self.joint[0]['denominator_df']=31;self.save()
        with self.assertRaises(ValueError):build_tables(self.work)
    def test_export_refuses_overwrite(self):
        export(self.work)
        before=(self.work/'displays/age_comparison_tables.json').read_bytes()
        with self.assertRaises(FileExistsError):export(self.work)
        self.assertEqual(before,(self.work/'displays/age_comparison_tables.json').read_bytes())

if __name__=='__main__':unittest.main()
