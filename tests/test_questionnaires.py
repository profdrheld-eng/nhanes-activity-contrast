import unittest
import numpy as np
import pandas as pd
from nhanes_activity.questionnaires import hip_self_report, wrist_self_report, aggregate_wrist

class HipScoringTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame({'SEQN':[1],'PAD020':[2],'PAQ050Q':[np.nan],'PAQ050U':[np.nan],
          'PAD080':[np.nan],'PAQ100':[2],'PAD120':[np.nan],'PAD160':[np.nan],'PAD320':[2],'PAD200':[2]})
    def empty(self):return pd.DataFrame(columns=['SEQN','PADLEVEL','PADTIMES','PADDURAT'])

    def test_negative_gates_are_zero(self):
        d=hip_self_report(self.frame(),self.empty(),{1})
        self.assertTrue(d.SR_complete.iloc[0]);self.assertEqual(d.SR_day.iloc[0],0)

    def test_unable_is_missing(self):
        p=self.frame();p['PAD320']=3
        self.assertFalse(hip_self_report(p,self.empty(),{1}).SR_complete.iloc[0])

    def test_vigorous_weight_and_month_denominator(self):
        p=self.frame();p['PAD200']=1
        activity=pd.DataFrame({'SEQN':[1],'PADLEVEL':[2],'PADTIMES':[3],'PADDURAT':[20]})
        self.assertEqual(hip_self_report(p,activity,{1}).SR_day.iloc[0],4)

    def test_gate_row_conflict_fails(self):
        p=self.frame();p['PAD320']=1
        with self.assertRaises(ValueError):hip_self_report(p,self.empty(),{1})

class WristScoringTests(unittest.TestCase):
    def test_week_denominator_and_missing_gate(self):
        row={'SEQN':1}
        for gate,days,duration in [('PAQ605','PAQ610','PAD615'),('PAQ620','PAQ625','PAD630'),
           ('PAQ635','PAQ640','PAD645'),('PAQ650','PAQ655','PAD660'),('PAQ665','PAQ670','PAD675')]:
            row.update({gate:2,days:np.nan,duration:np.nan})
        row.update({'PAQ605':1,'PAQ610':2,'PAD615':35})
        self.assertEqual(wrist_self_report(pd.DataFrame([row])).SR_day.iloc[0],20)
        row['PAQ635']=7
        self.assertFalse(wrist_self_report(pd.DataFrame([row])).SR_complete.iloc[0])

    def test_device_boundary_and_mean(self):
        d=pd.DataFrame({'SEQN':[1]*5,'PAXDAYD':[1,2,3,4,5],'PAXWWMD':[600,700,800,900,599],
                        'PAXVMD':[1400]*5,'PAXMTSD':[100,200,300,400,500]})
        x=aggregate_wrist(d).iloc[0]
        self.assertEqual(x.PAX_valid_days,4);self.assertTrue(x.PAX_primary_eligible)
        self.assertEqual(x.PAX_mean_daily_mims,250)
        with self.assertRaises(ValueError):aggregate_wrist(pd.concat([d,d.iloc[:1]]))

if __name__=='__main__':unittest.main()
