import unittest
import numpy as np
import pandas as pd
from nhanes_activity.pax import classify_nonwear, participant_qc, derive_variant, participants_from_chunks

class NonwearTests(unittest.TestCase):
    def mask(self, positions, values=None, minimum=60):
        known=np.zeros(10080,dtype=bool);counts=np.full(10080,np.nan)
        known[positions]=True;counts[positions]=0 if values is None else values
        return classify_nonwear(known,counts,minimum)

    def test_duration_threshold(self):
        self.assertEqual(self.mask(slice(0,59)).sum(),0)
        self.assertEqual(self.mask(slice(0,60)).sum(),60)

    def test_short_bounded_interruptions(self):
        v=np.zeros(60);v[10:12]=50;v[30]=99
        self.assertEqual(self.mask(slice(0,60),v).sum(),60)
        v[10:13]=50
        self.assertEqual(self.mask(slice(0,60),v).sum(),0)

    def test_low_run_needs_zero_on_both_sides(self):
        v=np.zeros(60);v[0]=50
        self.assertEqual(self.mask(slice(0,60),v).sum(),0)

    def test_unknown_and_day_boundary_split(self):
        self.assertEqual(self.mask(np.r_[0:30,31:61]).sum(),0)
        self.assertEqual(self.mask(slice(1410,1470)).sum(),0)

    def test_100_counts_not_permitted(self):
        v=np.zeros(60);v[30]=100
        self.assertEqual(self.mask(slice(0,60),v).sum(),0)

    def test_shape_and_threshold_guard(self):
        for n in [0,1440,10079]:
            with self.assertRaises(ValueError):classify_nonwear(np.ones(n,bool),np.zeros(n),60)
        for minimum in [0,-1,1.5,True]:
            with self.assertRaises(ValueError):self.mask(slice(0,60),minimum=minimum)

class ParticipantTests(unittest.TestCase):
    def test_valid_day_boundary_and_calibration(self):
        n=1440*4;p=np.arange(n)
        rows={'PAXN':p+1,'PAXDAY':p//1440+1,'PAXHOUR':p%1440//60,'PAXMINUT':p%60,
              'PAXSTAT':np.ones(n),'PAXCAL':np.ones(n),'PAXINTEN':np.full(n,200.),'PAXSTEP':np.full(n,np.nan)}
        reasons,_=participant_qc(rows);self.assertEqual(reasons,[])
        days,person=derive_variant(rows,'primary',{1},60,4)
        self.assertTrue(person['eligible_participant']);self.assertEqual(person['intensity_valid_days'],4)
        self.assertEqual(person['counts_per_wear_minute'],200)
        self.assertIsNone(person['mean_daily_steps'])
        rows['PAXCAL'][:]=2
        _,person=derive_variant(rows,'primary',{1},60,4);self.assertFalse(person['eligible_participant'])
        _,person=derive_variant(rows,'calibration_1_or_2',{1,2},60,4);self.assertTrue(person['eligible_participant'])
        rows['PAXN'][1]=rows['PAXN'][0]
        self.assertIn('duplicate_seqn_paxn',participant_qc(rows)[0])

class ChunkTests(unittest.TestCase):
    def test_participant_across_chunk_boundary(self):
        chunks=[pd.DataFrame({'SEQN':[1,1],'value':[10,11]}),pd.DataFrame({'SEQN':[1,2],'value':[12,20]})]
        result=list(participants_from_chunks(chunks))
        self.assertEqual([key for key,_ in result],[1,2])
        self.assertEqual(result[0][1].value.tolist(),[10,11,12])

    def test_interleaved_and_noninteger_keys_fail(self):
        for keys in [[1,2,1],[1,1.5],[1,float('nan')]]:
            with self.assertRaises(ValueError):list(participants_from_chunks([pd.DataFrame({'SEQN':keys})]))

    def test_cross_chunk_order_fails(self):
        with self.assertRaises(ValueError):list(participants_from_chunks([pd.DataFrame({'SEQN':[1,2]}),pd.DataFrame({'SEQN':[1,3]})]))

if __name__=='__main__':unittest.main()
