from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from nhanes_activity.cohort_support import derive_covariates, parse_mortality

class CovariateTests(unittest.TestCase):
    def test_disease_unknown_and_borderline(self):
        keys=[1,2,3]
        demo=pd.DataFrame({'SEQN':keys,'RIDSTATR':[2]*3,'RIDAGEYR':[50]*3,'RIAGENDR':[1,2,1],
            'RIDRETH1':[3]*3,'DMDEDUC2':[3]*3,'INDFMPIR':[1.5]*3,'WTMEC2YR':[1]*3,'SDMVSTRA':[1]*3,'SDMVPSU':[1]*3})
        smq=pd.DataFrame({'SEQN':keys,'SMQ020':[2,1,1],'SMQ040':[np.nan,3,1]})
        bmx=pd.DataFrame({'SEQN':keys,'BMXBMI':[22,25,np.nan]})
        diq=pd.DataFrame({'SEQN':keys,'DIQ010':[3,1,9]})
        mcq=pd.DataFrame({'SEQN':keys,'MCQ220':[1,2,9],
            'MCQ160B':[1,2,2],'MCQ160C':[np.nan,2,np.nan],'MCQ160D':[np.nan,2,2],
            'MCQ160E':[np.nan,2,2],'MCQ160F':[np.nan,2,2]})
        d=derive_covariates(demo,smq,bmx,mcq,diq)
        self.assertEqual(d.diabetes_history.iloc[0],'no')
        self.assertTrue(pd.isna(d.diabetes_history.iloc[2]))
        self.assertEqual(d.CVD_history.iloc[0],'yes')
        self.assertEqual(d.CVD_history.iloc[1],'no')
        self.assertTrue(pd.isna(d.CVD_history.iloc[2]))
        self.assertEqual(d.cancer_history.iloc[0],'yes')
        self.assertTrue(pd.isna(d.cancer_history.iloc[2]))
        self.assertEqual(d.smoking3.tolist(),['never','former','current'])
        with self.assertRaises(ValueError):derive_covariates(pd.concat([demo,demo.iloc[:1]]),smq,bmx,mcq,diq)

    def test_fixed_width_mortality_fields(self):
        fields=[(0,6,'101'),(14,15,'1'),(15,16,'0'),(16,19,'.'),(19,20,'.'),(20,21,'.'),(42,45,'120'),(45,48,'119')]
        chars=[' ']*48
        for start,end,value in fields:chars[start:end]=value.rjust(end-start)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'mortality.dat';path.write_text(''.join(chars)+'\n')
            d=parse_mortality(path).iloc[0]
            self.assertEqual(d.SEQN,101);self.assertEqual(d.MORTSTAT,0)
            self.assertEqual(d.PERMTH_INT,120);self.assertEqual(d.PERMTH_EXM,119)
            self.assertTrue(pd.isna(d.UCOD_LEADING))

if __name__=='__main__':unittest.main()
