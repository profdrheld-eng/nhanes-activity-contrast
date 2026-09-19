"""Validate S37 assembly and reject incomplete diagnostic inputs."""
import csv
import tempfile
import unittest
from pathlib import Path
from export_reporting_diagnostics import build_table, export

GROUPS=['2003-2004','2005-2006','2011-2012','2013-2014','hip','wrist']
ADJUSTMENTS=['unadjusted','level_only','principal']

class DiagnosticExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.work=Path(self.tmp.name);self.path=self.work/'reporting_diagnostics';self.path.mkdir()
        (self.path/'COMPLETE.txt').write_text('18 model diagnostics complete\n')
        self.ph=[];self.steps=[];self.support=[]
        for group in GROUPS:
            for adjustment in ADJUSTMENTS:
                base=dict(model=group+'__'+adjustment,group=group,adjustment=adjustment)
                self.ph.extend([dict(base,term=term,p=.0004) for term in ['Delta','GLOBAL']])
                self.steps.extend([dict(base,contrast=c,inference='full_design_t',n=100,deaths=10,
                    estimate=1.2,ci_low=1.05,ci_high=1.5,p=.025,cut_years=5,design_df=15)
                    for c in ['early','late','late_to_early']])
                self.support.append(dict(model=base['model'],n=100,deaths=10,design_df=15))
        self.write()

    def write(self):
        for name,rows in [('schoenfeld',self.ph),('step_time',self.steps),('support',self.support)]:
            with (self.path/(name+'.csv')).open('w',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

    def test_complete_and_order_independent(self):
        table=build_table(self.work)
        self.assertEqual(len(table),19)
        self.assertEqual(table[1],['2003-2004 / Unadjusted','<0.001','<0.001','1.200 [1.050, 1.500]','0.025'])
        self.steps.reverse();self.ph.reverse();self.write()
        self.assertEqual(build_table(self.work),table)

    def test_missing_duplicate_unexpected_models(self):
        for mutate in [lambda: self.steps.pop(),lambda:self.steps.append(self.steps[0]),
                       lambda:self.steps[0].update(model='unexpected')]:
            original=[r.copy() for r in self.steps];mutate();self.write()
            with self.assertRaises(ValueError):build_table(self.work)
            self.steps=original

    def test_invalid_numbers_and_support(self):
        for key,value in [('p','nan'),('p',1.1),('ci_low',2),('estimate',0),('n',101),
                          ('cut_years',4),('inference','normal')]:
            row=self.steps[0].copy();self.steps[0][key]=value;self.write()
            with self.assertRaises(ValueError):build_table(self.work)
            self.steps[0]=row
        self.ph[0]['p']=-1;self.write()
        with self.assertRaises(ValueError):build_table(self.work)

    def test_completion_and_overwrite_guards(self):
        (self.path/'COMPLETE.txt').unlink()
        with self.assertRaises(ValueError):build_table(self.work)
        (self.path/'COMPLETE.txt').write_text('18 model diagnostics complete\n')
        export(self.work)
        target=self.work/'displays/table_s37.csv';before=target.read_bytes()
        with self.assertRaises(FileExistsError):export(self.work)
        self.assertEqual(target.read_bytes(),before)

if __name__=='__main__':unittest.main()
