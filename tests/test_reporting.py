import unittest
import numpy as np
import pandas as pd
from nhanes_activity.reporting import missing_counts, quintile_cuts, assign_quintiles, followup_summary, composite

class ReportingTests(unittest.TestCase):
    def test_nonsequential_missingness_and_zero(self):
        d=pd.DataFrame({'a':[None,0,2],'b':[None,False,True]})
        rows=missing_counts(d,['a','b'])
        self.assertEqual([r['missing'] for r in rows],[1,1])
        self.assertEqual([r['denominator'] for r in rows],[3,3])
    def test_weighted_empirical_cutpoints_and_ties(self):
        x=np.array([0.,1.,2.,3.,4.]);w=np.array([1.,1.,6.,1.,1.])
        cuts=quintile_cuts(x,w)
        np.testing.assert_array_equal(cuts,[1.,2.,2.,2.])
        np.testing.assert_array_equal(assign_quintiles(x,cuts),[1,1,2,5,5])
    def test_reject_bad_weights(self):
        with self.assertRaises(ValueError):quintile_cuts([1,2],[0,0])
    def test_person_years_unweighted(self):
        d=pd.DataFrame({'follow_up_years':[2.,4.],'MORTSTAT':[1,0],'WTMEC2YR':[100.,1.]})
        x=followup_summary(d)
        self.assertEqual(x['person_years'],6.)
        self.assertEqual(x['mean_years'],3.)
        self.assertEqual(x['deaths'],1)
    def test_positive_composite_overrides_unknown(self):
        d=pd.DataFrame({'a':['yes','no',None,'no'],'b':[None,'no','no',None]})
        v=composite(d,['a','b'])
        self.assertEqual(v.iloc[0],'yes');self.assertEqual(v.iloc[1],'no')
        self.assertTrue(v.iloc[2:].isna().all())

if __name__=='__main__':unittest.main()
